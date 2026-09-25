"""LLM-loop tools for the OSP PBPK benchmark: inspect + try-a-model.

Two tools give Claude the closed loop:

  * ``osp_inspect``   (observe) - the current model (parameters, distribution
    method, processes), the literature priors, and the observed-data overview.
  * ``osp_try_model`` (act)     - apply an edit spec, run PK-Sim headless, and
    return the fit (GMFE, %-within-2-fold, and per-route BIAS so the agent knows
    which way to move each parameter), plus plausibility flags.

The task context (the base snapshot to edit, the observed data, the PK-Sim CLI)
is captured in the handlers, so these are registered per-task, not globally.
"""

from __future__ import annotations

import json
import os
from typing import Any

from ..engines.osp_cli import OSPCli
from ..engines import osp_score
from ..engines import osp_catalog
from ..engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
from .registry import Tool, ToolRegistry, ToolResult


def _expressed_molecules(snapshot_path: str) -> list[dict]:
    with open(snapshot_path, encoding="utf-8") as fh:
        data = json.load(fh)
    seen, out = set(), []
    for ep in data.get("ExpressionProfiles") or []:
        mol = ep.get("Molecule")
        if mol and mol not in seen:
            seen.add(mol)
            out.append({"molecule": mol, "type": ep.get("Type")})
    return out


_GIVEN_SOURCES = {"Publication", "In Vitro", "In vitro", "Database", "Internet"}


def _value_status(value_origin) -> str:
    """Classify a parameter's CURRENT value by its provenance, so the agent knows
    which starting numbers to trust. A measured/published value is a real input;
    a blanked benchmark placeholder is only a naive starting default the agent
    must DETERMINE (from the data or established physchem knowledge) - not trust.
    Everything else is a structural value present in the model as-is (e.g. a
    physical constant like molecular weight, which is correct and trustworthy).
    This is problem-setup information (which inputs are given vs unknown), not the
    answer: it never reveals the target value."""
    vo = value_origin or {}
    src = vo.get("Source")
    if src in _GIVEN_SOURCES:
        return "given"          # measured / published - trust it
    desc = str(vo.get("Description") or "").lower()
    if src == "Unknown" and ("blanked" in desc or "naive prior" in desc):
        return "placeholder"    # a benchmark unknown - determine it, do not trust
    # untagged / structural value present in the model as-is (constants, defaults)
    return "structural"


def _current_model(snapshot_path: str) -> dict[str, Any]:
    with open(snapshot_path, encoding="utf-8") as fh:
        _snap = json.load(fh)
    comp = (_snap.get("Compounds") or [{}])[0]
    _forms = _snap.get("Formulations") or []
    params = []
    seen = set()

    # compound-level params (everything except the Processes block)
    def walk(o):
        if isinstance(o, dict):
            nm = o.get("Name")
            if isinstance(nm, str) and isinstance(o.get("Value"), (int, float)) \
                    and nm not in seen:
                seen.add(nm)
                params.append({"name": nm, "value": o["Value"],
                               "unit": o.get("Unit", ""),
                               "value_status": _value_status(o.get("ValueOrigin"))})
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk({k: v for k, v in comp.items() if k != "Processes"})

    # process params: if the SAME parameter name appears on more than one process
    # (e.g. a per-enzyme CLspec on UGT1A9/UGT2B7/CYP), expose each with a QUALIFIED
    # name '<Name>@<Molecule>' so the agent can set/estimate them independently.
    from collections import defaultdict
    occ = defaultdict(list)
    for p in comp.get("Processes") or []:
        mol = p.get("Molecule") or p.get("InternalName")
        pnames = {par.get("Name") for par in p.get("Parameters") or []}
        # for a specific-clearance metabolization process, 'CLspec/[Enzyme]' is
        # THE fittable clearance; hide its structural/derived siblings so the
        # agent doesn't estimate 'Specific clearance' (a derived, usually-0 value)
        # or 'Enzyme concentration' by mistake.
        skip = set()
        if "CLspec/[Enzyme]" in pnames:
            skip = {"Specific clearance", "Enzyme concentration"}
        for par in p.get("Parameters") or []:
            nm = par.get("Name")
            if nm in skip:
                continue
            if isinstance(nm, str) and isinstance(par.get("Value"), (int, float)):
                occ[nm].append((mol, par["Value"], par.get("Unit", ""),
                                _value_status(par.get("ValueOrigin"))))
    for nm, os_ in occ.items():
        if len(os_) == 1:
            if nm not in seen:
                seen.add(nm)
                params.append({"name": nm, "value": os_[0][1], "unit": os_[0][2],
                               "value_status": os_[0][3]})
        else:                                   # collision -> qualify by molecule
            for mol, val, unit, status in os_:
                params.append({"name": f"{nm}@{mol}", "value": val, "unit": unit,
                               "on_process": mol, "value_status": status})

    # FORMULATION params (Weibull dissolution = oral absorption). Expose them as
    # editable so the agent can FIT absorption. With more than one formulation
    # (tablet vs capsule, different release), qualify by formulation name
    # ('Dissolution time (50% dissolved)@Halcion') so each can be fit independently.
    multi = len([f for f in _forms if f.get("Parameters")]) > 1
    for form in _forms:
        fname = form.get("Name")
        for par in form.get("Parameters") or []:
            nm = par.get("Name")
            if not (isinstance(nm, str) and isinstance(par.get("Value"), (int, float))):
                continue
            if "dissolution" not in nm.lower():        # only the fittable dissolution knobs
                continue
            key = f"{nm}@{fname}" if multi else nm
            if key in seen:
                continue
            seen.add(key)
            params.append({"name": key, "value": par["Value"], "unit": par.get("Unit", ""),
                           "on_formulation": fname,
                           "value_status": _value_status(par.get("ValueOrigin"))})

    return {
        "parameters": params,
        "calculation_methods": comp.get("CalculationMethods") or [],
        "processes": [{"molecule": p.get("Molecule"),
                       "internal": p.get("InternalName")}
                      for p in comp.get("Processes") or []],
    }


def _observed_overview(observed: list[dict]) -> dict[str, Any]:
    routes: dict[str, dict] = {}
    for o in observed:
        r = (o.get("route") or "NA")
        routes.setdefault(r, {"n_datasets": 0, "studies": set(), "doses": set()})
        routes[r]["n_datasets"] += 1
        if o.get("study"):
            routes[r]["studies"].add(o["study"])
        if o.get("dose"):
            routes[r]["doses"].add(str(o["dose"]))
    for r in routes.values():
        r["studies"] = sorted(r["studies"])
        r["doses"] = sorted(r["doses"])
    return {"n_datasets": len(observed), "by_route": routes}


def _reference_index(reflib: "dict | None") -> "dict | None":
    """Turn the assembled reference library into a bare INDEX for osp_inspect:
    which analogues exist and which enzymes/transporters each involves - but NOT
    their methods or parameter names. To actually use an analogue's structure the
    agent must open it with osp_read_reference, so it reads rather than skims."""
    if not reflib:
        return None
    idx = []
    for m in reflib.get("models") or []:
        mols = set()
        for p in m.get("processes") or []:
            if isinstance(p, dict) and p.get("molecule"):
                mols.add(p["molecule"])
        for e in m.get("estimated_parameters") or []:
            pn = e if isinstance(e, str) else (e.get("parameter") or "")
            if "@" in pn:
                mols.add(pn.split("@", 1)[1])
        idx.append({"compound": m.get("compound"), "involves": sorted(mols)})
    return {"mode": reflib.get("mode"),
            "note": (reflib.get("note", "") + " This is only an INDEX. To USE an "
                     "analogue you MUST open it: call osp_read_reference(compound) "
                     "to read its full structure (distribution/permeability method, "
                     "processes, fitted parameter names)."),
            "available": idx}


def _reference_markdown(m: dict) -> str:
    methods = m.get("calculation_methods") or []
    procs = m.get("processes") or []
    est = m.get("estimated_parameters") or []
    lines = [f"# Reference model: {m.get('compound')} (analogue, leave-one-out)",
             "", "## Calculation methods"]
    lines += [f"- {x}" for x in methods] or ["- (none listed)"]
    lines += ["", "## Processes (mechanisms)"]
    if procs:
        for p in procs:
            mol = p.get("molecule") if isinstance(p, dict) else None
            ty = p.get("type") if isinstance(p, dict) else None
            lines.append(f"- {mol or '?'} — {ty or '?'}")
    else:
        lines.append("- (none)")
    lines += ["", "## Parameters the modeler ESTIMATED (fitted)"]
    if est:
        for e in est:
            if isinstance(e, dict):
                v = e.get("value")
                lines.append(f"- {e.get('parameter')}" + (f" = {v}" if v is not None else ""))
            else:
                lines.append(f"- {e}")
    else:
        lines.append("- (none)")
    lines += ["", "Reuse this STRUCTURE by analogy where the chemistry matches; fit "
              "the numbers to THIS compound's own data (do not copy values blindly)."]
    return "\n".join(lines)


def _method_guidance(observed: list[dict]) -> dict[str, Any]:
    """Best-practice METHOD guidance derived only from the data shape (not the
    answer): the staged IV->PO fit (Kuepfer 2016) and whether saturable kinetics
    are even identifiable from the doses on hand. Guidance, not instructions."""
    def _r(o):
        s = str(o.get("route") or "").upper()
        return "IV" if "IV" in s else ("PO" if ("PO" in s or "ORAL" in s) else "other")
    routes = {_r(o) for o in observed}
    doses = {str(o.get("dose")) for o in observed if o.get("dose")}
    out: dict[str, Any] = {}
    if "IV" in routes and "PO" in routes:
        out["route_staging"] = (
            "Both IV and PO data are present. Fit in STAGES: (1) on the IV data "
            "alone, choose the distribution method and fit distribution + systemic "
            "clearance; (2) then on the PO data, HOLD those fixed and fit ONLY the "
            "absorption parameters (intestinal permeability, dissolution/solubility, "
            "lag/meal). Fitting all routes jointly lets absorption and clearance "
            "compensate and can hide a structural gap. (If an IV-relevant process is "
            "also an absorption process, e.g. a gut uptake transporter, fit IV+PO "
            "together instead.)")
    elif "PO" in routes and "IV" not in routes:
        out["route_staging"] = (
            "PO data only: absorption and clearance are not separately identifiable "
            "without IV data — fix clearance from a prior/IVIVE where you can, and "
            "treat the absolute clearance as uncertain.")
    n_doses = len(doses)
    out["dose_levels"] = n_doses
    if n_doses < 2:
        out["saturable_identifiability"] = (
            f"Only {n_doses} dose level in the data. Michaelis–Menten / saturable "
            "kinetics (Km, Vmax), saturable transport, or dose-nonlinear absorption "
            "are NOT identifiable without ≥2 dose levels — prefer linear clearance "
            "unless a prior fixes the saturable constants.")
    return out


def _structure_problems(structure, estimate, link_groups=None) -> list[str]:
    """Catch model-config patterns that make PK-Sim build/run produce no output,
    so the agent gets an actionable message instead of the opaque 'no result CSVs
    parsed'. Pure config checks - no engine needed."""
    structure = structure or {}
    estimate = estimate or {}
    problems: list[str] = []
    # (1) two+ added processes on the SAME molecule -> their per-molecule parameter
    #     keys collide (both e.g. kcat@CYP3A4 / Intrinsic clearance@CYP3A4), so
    #     PK-Sim cannot build them apart and the run yields nothing.
    counts: dict[str, int] = {}
    for p in structure.get("add_processes") or []:
        mol = (p or {}).get("molecule")
        if mol:
            counts[mol] = counts.get(mol, 0) + 1
    for mol, n in counts.items():
        if n > 1:
            problems.append(
                f"you added {n} processes on the same molecule '{mol}': their "
                f"parameters share one key (kcat@{mol} / Intrinsic clearance@{mol}), "
                "so PK-Sim cannot distinguish the parallel routes and the run "
                f"produces no output. Model the parallel {mol} routes as ONE lumped "
                f"pathway — a single process on {mol} with one total clearance. "
                "(A single plasma curve cannot split parallel routes on the same "
                f"enzyme anyway.) Do NOT move a second {mol} route onto a different "
                "molecule to dodge this — that models the wrong biology.")
    # (2) link_scale group with fewer than two DISTINCT members (e.g. the same key
    #     twice) - degenerate, and usually rides on the duplicate-molecule bug.
    for members in (link_groups or []):
        if len(set(members)) < 2:
            problems.append(
                f"a link_scale group has fewer than two distinct members ({members}); "
                "link_scale fits the shared magnitude of >=2 DIFFERENT parameters - "
                "pass distinct keys, or fit the single parameter directly with "
                "'estimate'.")
    # (3) fitting 'Permeability' while a permeability method that COMPUTES it is
    #     active -> Permeability is not a free parameter, so PK-Sim fails to run.
    perm = (structure.get("calculation_methods") or {}).get("permeability")
    if "Permeability" in estimate and perm and "Charge dependent" in str(perm):
        problems.append(
            f"'Permeability' is being fitted but the permeability method '{perm}' "
            "COMPUTES cellular permeability from physchem, so Permeability is not a "
            "settable/fittable parameter under it and PK-Sim fails to run. Switch "
            "permeability to 'PK-Sim Standard' to fit Permeability, or remove it "
            "from estimate.")
    return problems


def _given_physchem_value(base: str, lit: list) -> "float | None":
    """The GIVEN value of a DIRECT physicochemical property (currently lipophilicity) if the input
    provided it, else None. Only direct measurements are returned - NOT in-vitro kinetic inputs
    (CLint/Vmax/Km), which are legitimately refined. Used to FIX a given measurement rather than let
    the agent re-fit it (fitting a given logP frees a knob that makes the distribution method
    unidentifiable - the Tizanidine failure)."""
    b = (base or "").lower()
    if not ("lipophil" in b or b in ("logp", "logd")):
        return None
    for e in lit or []:
        n = (e.get("parameter") or "").lower()
        if "lipophil" in n or "logp" in n or "logd" in n:
            v = e.get("value")
            if isinstance(v, (int, float)):
                return float(v)
    return None


# Lipophilicity is an EFFECTIVE distribution parameter in PK-Sim, so its best-fit
# value commonly differs from the measured logP by up to ~1 log unit; the fit is
# constrained to a WINDOW this wide around the given logP rather than pinned at it.
_LIP_WINDOW = 1.0


def _is_lipophilicity(base: str) -> bool:
    b = (base or "").lower()
    return "lipophil" in b or b in ("logp", "logd")


def _fix_given_physchem(estimate: dict, fix: dict, lit: list):
    """Constrain estimate parameters whose value the input GAVE as a direct measurement.

    A true measurement (fraction unbound, pKa, solubility) is FIXED at its value.
    LIPOPHILICITY is the exception: in PK-Sim it is an EFFECTIVE distribution
    parameter that commonly sits up to ~1 log unit off the measured logP (the
    reference alfentanil model fits 1.85 against a measured 2.16), so pinning it
    at the measured value biases Vd and can make the model unfittable no matter
    which distribution method is swept. Instead we keep it ESTIMABLE within a
    bounded window (measured +/- 1 log unit) - narrow enough to keep the method
    identifiable, wide enough to let the effective value settle. Returns
    (estimate, fix, notes)."""
    estimate, fix, notes = dict(estimate), dict(fix or {}), []
    for name in list(estimate):
        base = name.split("@", 1)[0]
        gv = _given_physchem_value(base, lit)
        if gv is None:
            continue
        if _is_lipophilicity(base):
            lo, hi = gv - _LIP_WINDOW, gv + _LIP_WINDOW
            req = estimate.get(name)
            if isinstance(req, (list, tuple)) and len(req) == 2:
                # intersect with the agent's requested bounds, but never collapse
                clo, chi = max(lo, float(req[0])), min(hi, float(req[1]))
                if chi > clo:
                    lo, hi = clo, chi
            estimate[name] = [lo, hi]
            notes.append(f"{name}: measured logP {gv:g} given, but Lipophilicity is an "
                         f"EFFECTIVE distribution parameter - fitted within a bounded window "
                         f"[{lo:g}, {hi:g}] around the measured value (not pinned at it, which "
                         f"would bias Vd).")
        else:
            estimate.pop(name)
            fix[name] = gv
            notes.append(f"{name}: a measured value ({gv:g}) was given in the input - FIXED at it, "
                         f"not estimated (a true measurement is used, not re-fitted).")
    return estimate, fix, notes


def _defaulted_process_params(add_processes, provided_keys) -> list[dict]:
    """Added-process parameters left at their CATALOG DEFAULT - i.e. the agent
    added the process but neither gave the parameter a value in the process spec
    nor fitted/fixed it. These defaults are otherwise injected SILENTLY (GFR
    fraction 1.0 = full filtration, Km 1.0, Intrinsic clearance 0.1), so a default
    can quietly become the model value. Surfacing them lets the agent SEE what it
    left unset and decide to fit or set it. `provided_keys` = every parameter name
    the agent set or chose to fit (base name or '<param>@<molecule>')."""
    prov = {str(k).lower() for k in (provided_keys or [])}
    out = []
    for p in add_processes or []:
        if not isinstance(p, dict):
            continue
        spec = osp_catalog.PROCESS_TYPES.get(p.get("type")) or {}
        given = {str(k).lower() for k in (p.get("parameters") or {})}
        mol = str(p.get("molecule") or "")
        for prm in spec.get("parameters", []):
            base = str(prm["name"]).lower()
            qual = f"{base}@{mol.lower()}" if mol else base
            if base in given or base in prov or qual in prov:
                continue
            out.append({"molecule": mol or None, "parameter": prm["name"],
                        "default": prm.get("default")})
    return out


def _given_measurement(base: str, lit: list) -> bool:
    """True if the input's literature_physicochemical GIVES a measured value for this parameter
    (matched by kind), so it must be respected - never widened or freely re-estimated. General: keys
    off the parameter kind, not any specific drug."""
    b = (base or "").lower()
    keys = [b]
    if "lipophil" in b:
        keys += ["lipophil", "logp", "logd"]
    elif "unbound" in b:
        keys += ["unbound", "fraction unbound", "fu", "fup"]
    elif "gfr" in b:
        keys += ["gfr"]
    names = [(e.get("parameter") or "").lower() for e in (lit or [])]
    return any(any(k in n for k in keys) for n in names)


def register_osp_loop_tools(registry: ToolRegistry, config, ctx: dict) -> None:
    """ctx: {cli: OSPCli, snapshot_path, observed: [...], input: {...}}."""
    cli: OSPCli = ctx["cli"]
    snapshot_path: str = ctx["snapshot_path"]
    observed: list[dict] = ctx["observed"]
    inp: dict = ctx.get("input") or {}
    # SELF-EXTRACT: the agent extracts givens itself, so we must NOT hand it the
    # reference model's fix-vs-fit split (which params it fit vs took as given) -
    # that is a chosen method. The agent decides from its own recorded givens.
    self_extract: bool = bool(ctx.get("self_extract"))
    # FIT-VISION: hand the model a rendered semilog overlay of its own predicted
    # curves vs the observed data, so it reasons about SHAPE (missing distribution
    # phase, wrong terminal slope, Cmax/tmax offset) the scalar GMFE hides. On by
    # default; the copilot toggle can turn it off for an ablation.
    fit_vision: bool = bool(ctx.get("fit_vision", True))

    def _hide_status(model: dict) -> dict:
        return {**model, "parameters": [{k: v for k, v in p.items() if k != "value_status"}
                                        for p in model.get("parameters", [])]}

    def _fit_images(series: list, title: str) -> list:
        """Render the observed-vs-simulated overlay to a PNG image block, when
        fit-vision is on and there is something plottable. Never raises into a run."""
        if not fit_vision or not series:
            return []
        try:
            from ..engines import fit_plot
            im = fit_plot.render_fit_b64(series, title=title)
            return [im] if im else []
        except Exception:                            # noqa: BLE001 - vision is best-effort
            return []

    def _cap_evals(n: int) -> int:
        """Clamp a requested optimizer budget to config.max_evals_cap (if set), so a
        big multi-compound run stays affordable - each eval is a full PK-Sim rebuild."""
        cap = getattr(config, "max_evals_cap", None)
        return min(n, cap) if cap else n

    # -- observe -------------------------------------------------------- #
    def inspect(args: dict, session) -> ToolResult:
        gd = inp.get("given_data", {}) or {}
        bg = inp.get("background") or {}
        model = _current_model(snapshot_path)
        # up-front: which starting values are real inputs vs placeholders to
        # determine, so the agent does not trust a naive default as if measured.
        to_determine = [p["name"] for p in model["parameters"]
                        if p.get("value_status") == "placeholder"]
        given_in_model = [p["name"] for p in model["parameters"]
                          if p.get("value_status") == "given"]
        out = dict(
            objective=inp.get("objective"),
            background=bg.get("description"),
            known_biology=bg.get("literature_facts"),
            # mechanism-DISCOVERY mode: the enzyme identity is withheld and the
            # agent must pick the clearing molecule(s) out of this candidate pool.
            candidate_clearance_molecules=bg.get("candidate_clearance_molecules"),
            compound_identity=gd.get("compound_identity"),
            literature_physicochemical=gd.get("literature_physicochemical")
            or inp.get("literature_physicochemical"),
            unknowns_guidance=inp.get("unknowns_guidance"),
            # LIBRARY-ASSISTED mode: an INDEX of finished models of OTHER compounds
            # (leave-one-out); the agent must osp_read_reference to see any one's
            # structure. None in de-novo mode.
            reference_library=_reference_index(inp.get("reference_library")),
            evaluation_rubric=inp.get("evaluation_rubric"),
            current_model=model,
            observed_overview=_observed_overview(observed),
            method_guidance=_method_guidance(observed),   # staged IV->PO + saturable identifiability
            valid_partition_methods=PARTITION_METHODS,
            valid_permeability_methods=PERMEABILITY_METHODS,
        )
        if self_extract:
            # withhold the reference's fix-vs-fit split and the given/placeholder
            # labels; the agent decides from ITS recorded givens + the tiers.
            out["current_model"] = _hide_status(model)
            out["fix_vs_fit_note"] = (
                "You decide the split: FIX each parameter you recorded as a measured "
                "given (via osp_record_givens) at that value, and FIT the rest that "
                "are uncertain. The task does NOT tell you which the reference fit vs "
                "fixed - use the parameter tiers in osp_options (tier=estimate is "
                "normally fit; tier=measured_soft/constant is normally fixed) and your "
                "own givens.")
        else:
            out["parameters_to_determine"] = to_determine
            out["given_input_parameters"] = given_in_model
        return ToolResult.success(
            "task description: objective, known biology, literature priors, "
            "current model, and the observed clinical data", **out)

    # -- act (run + score) ---------------------------------------------- #
    def try_model(args: dict, session) -> ToolResult:
        edits = args.get("edits") or {}
        res = cli.build_and_run(snapshot_path, edits=edits)
        if not res["ok"]:
            return ToolResult.error(
                f"PK-Sim run failed: {res['message']}",
                edits_applied=res.get("edits_applied"))

        predicted, unmatched = osp_score.map_predictions(res["profiles"], observed)
        score = osp_score.score_fit(observed, predicted)
        applied = res.get("edits_applied") or {}
        param_list = [{"parameter": k, "value": v, "unit": ""}
                      for k, v in applied.get("parameters", {}).items()]
        flags = osp_score.plausibility(param_list)

        overall = score["overall"]
        # track history + best
        hist = session.get("osp_history") or []
        hist.append({"edits": edits, "gmfe": overall.get("gmfe")})
        session.put("osp_history", hist)
        best = session.get("osp_best_gmfe")
        if overall.get("gmfe") is not None and (best is None or overall["gmfe"] < best):
            session.put("osp_best_gmfe", overall["gmfe"])
            session.put("osp_best_edits", edits)

        worst = [{"dataset": d["dataset"], "route": d["route"],
                  "gmfe": d["gmfe"], "bias": d["bias"]}
                 for d in score["per_dataset"][:3]]
        pkp = score.get("pk_parameters") or {}
        overlay = osp_score.overlay_series(observed, predicted)
        images = _fit_images(overlay, f"observed vs simulated — GMFE {overall.get('gmfe')}")
        defaulted = _defaulted_process_params(
            edits.get("add_processes"),
            list((edits.get("parameters") or {}).keys()) + list((edits.get("fix") or {}).keys()))
        dflt_note = ""
        if defaulted:
            dflt_note = (" | ⚠ added-process parameter(s) left at their DEFAULT (not set "
                         "or fitted): "
                         + ", ".join(f"{d['parameter']}"
                                     + (f"@{d['molecule']}" if d['molecule'] else "")
                                     + f"={d['default']}" for d in defaulted)
                         + " — a default silently defines the model; set or fit these.")
        out = ToolResult.success(
            f"GMFE {overall.get('gmfe')} overall "
            f"(within2fold {overall.get('within_2fold_pct')}%); "
            f"best so far {session.get('osp_best_gmfe')}"
            + (" — a fit overlay is attached; read it for shape mismatches "
               "(distribution phase, terminal slope, Cmax/tmax)." if images else "")
            + dflt_note,
            gmfe_overall=overall.get("gmfe"),
            within_2fold_pct=overall.get("within_2fold_pct"),
            bias_overall=overall.get("bias"),
            by_route=score["by_route"],
            # PK-parameter fold-errors (Cmax/tmax/AUC/t1/2) so shape & exposure,
            # not just the pointwise GMFE, guide the next edit.
            pk_parameter_gmfe=pkp.get("gmfe_by_metric"),
            worst_datasets=worst,
            parameter_flags=flags,
            defaulted_process_params=defaulted or None,   # params silently at their default
            edits_applied=applied,
            not_found=applied.get("not_found"),
            n_matched=len(predicted), n_total=len(observed),
            best_gmfe_so_far=session.get("osp_best_gmfe"),
            iteration=len(hist),
            saw_fit_curve=bool(images),               # transparency marker for the UI
        )
        out.images = images
        return out

    # -- observe (authoritative action space) --------------------------- #
    def options(args: dict, session) -> ToolResult:
        model = _current_model(snapshot_path)
        expressed = _expressed_molecules(snapshot_path)
        # mechanism-DISCOVERY (hard) mode: the clearing molecule is withheld, so
        # the "molecules you may attach a mechanism to" menu must be the FIXED,
        # case-independent discovery panel - NOT the model's own expressed set,
        # which is exactly the answer and would leak it. Every hard task shows
        # the same broad panel; any candidate is attachable (materialized on
        # demand), so the true molecule sits among decoys.
        hard_pool = (inp.get("background") or {}).get("candidate_clearance_molecules")
        if hard_pool:
            expressed = osp_catalog.hard_candidate_pool(hard_pool)
        with open(snapshot_path, encoding="utf-8") as fh:
            comp0 = (json.load(fh).get("Compounds") or [{}])[0]
        is_small = comp0.get("IsSmallMolecule", True)
        molecule_type = "small molecule" if is_small else (
            "large molecule / protein (biologic): disposition is via size-limited "
            "distribution + FcRn recycling + target binding, NOT enzyme processes; "
            "tune the compound parameters (radius, FcRn Kd, target binding)")
        lit = (inp.get("given_data", {}) or {}).get("literature_physicochemical", [])
        editable = []
        for p in model["parameters"]:
            # a qualified per-process name ('CLspec/[Enzyme]@UGT1A9') is described
            # by its BASE parameter for tier/range/role.
            base = p["name"].split("@", 1)[0]
            cat = osp_catalog.describe_parameter(base)
            tier = osp_catalog.param_tier(base)
            entry = {**p, "description": cat.get("description"),
                     "plausible_range": cat.get("range"),
                     "role": cat.get("role", "unknown"), "tier": tier}
            if self_extract:
                entry.pop("value_status", None)      # don't leak the reference fix-vs-fit split
            if "@" in p["name"]:
                entry["per_process"] = True     # set/estimate with this exact qualified name
            # measured_soft params carry their measured range (a fit may not exceed it).
            if tier == "measured_soft":
                mr = osp_catalog.measured_range(p["name"], lit)
                if mr:
                    entry["measured_range"] = list(mr)
            editable.append(entry)
        # method options are just names (the sweep chooses the method deterministically; the
        # agent does not need per-method prose). The measured-quantity principle and the
        # process-catalog rules live in the task PROMPT, not re-sent on every call.
        return ToolResult.success(
            "authoritative action space: what you may edit and legal choices. tier=constant "
            "is never estimated; tier=measured_soft is fixed by default (estimate only within "
            "measured_range); tier=estimate is the identifiable set to fit.",
            molecule_type=molecule_type,
            editable_parameters=editable,
            calculation_methods={
                "partition": {
                    "current": next((m for m in model["calculation_methods"]
                                     if "partition" in m.lower()), None),
                    "options": list(PARTITION_METHODS)},
                "permeability": {
                    "current": next((m for m in model["calculation_methods"]
                                     if "permeability" in m.lower()), None),
                    "options": list(PERMEABILITY_METHODS)}},
            processes_present=model["processes"],
            expressed_molecules=expressed,
            addable_process_types=osp_catalog.addable_process_types(expressed),
            interaction_process_types=osp_catalog.interaction_process_types(),
            edit_spec_help={
                "parameters": "{name: value}", "processes": "{molecule: false} to disable",
                "calculation_methods": "{partition:.., permeability:..}",
                "add_processes": "[{type, molecule, parameters}] (enzyme needs an expressed enzyme)"},
        )

    registry.register(Tool(
        name="osp_options",
        description=(
            "OBSERVE the authoritative ACTION SPACE for this model: every editable "
            "compound parameter (with a description, plausible range, and whether "
            "it is normally measured or estimated), the legal distribution/"
            "permeability calculation methods (with descriptions), the processes "
            "currently present, the molecules the model EXPRESSES (enzymes/"
            "transporters you can attach a mechanism to), and the process types "
            "you may ADD. Read straight from the model - you do not need prior OSP "
            "knowledge. Call this before deciding structure/parameters."),
        input_schema={"type": "object", "properties": {}},
        handler=options, phase="observe"))

    registry.register(Tool(
        name="osp_inspect",
        description=(
            "OBSERVE the PBPK task: returns the current compound model "
            "(parameters with values/units, the distribution & permeability "
            "calculation methods, and the metabolizing/clearance processes), the "
            "literature physicochemical priors, and an overview of the observed "
            "clinical datasets (routes, studies, doses). Call this first."),
        input_schema={"type": "object", "properties": {}},
        handler=inspect, phase="observe"))

    context_reports: dict = ctx.get("context_reports") or {}   # compound -> report .md path

    def read_reference(args: dict, session) -> ToolResult:
        reflib = inp.get("reference_library") or {}
        models = reflib.get("models") or []
        if not models:
            return ToolResult.error("No reference library for this task (de-novo "
                                    "mode) — decide the structure from first principles.")
        name = (args.get("compound") or "").strip()
        if not name:
            return ToolResult.error(
                "Name an analogue to read. Available: "
                + ", ".join(x.get("compound", "") for x in models))
        m = next((x for x in models
                  if (x.get("compound") or "").lower() == name.lower()), None)
        if not m:
            return ToolResult.error(
                f"Unknown analogue '{name}'. Available: "
                + ", ".join(x.get("compound", "") for x in models))
        # the analogue's FULL modeling report (its own published write-up), if handed
        report_md = None
        rp = context_reports.get(m.get("compound"))
        if rp and os.path.isfile(rp):
            try:
                report_md = open(rp, encoding="utf-8").read()[:20000]
            except OSError:
                report_md = None
        return ToolResult.success(
            f"finished model of the analogue {m.get('compound')} (a leave-one-out "
            "reference, NOT the target) — read its full report and structure and "
            "reuse it by analogy where the chemistry matches; fit the numbers to "
            "THIS compound.",
            compound=m.get("compound"),
            calculation_methods=m.get("calculation_methods"),
            processes=m.get("processes"),
            estimated_parameters=m.get("estimated_parameters"),
            report_markdown=report_md,
            reference_markdown=_reference_markdown(m))

    registry.register(Tool(
        name="osp_read_reference",
        description=(
            "OPEN and read one analogue's finished model from the reference library "
            "(by compound name, from the index osp_inspect returns). Returns that "
            "analogue's FULL modeling report (report_markdown — its own published "
            "write-up: methods, rationale, fitted values) plus a structured summary "
            "(distribution/permeability method, processes, fitted parameter names). "
            "You MUST open and READ the closest analogues this way before choosing "
            "your structure; the inspect index alone is not enough."),
        input_schema={"type": "object", "properties": {
            "compound": {"type": "string",
                         "description": "analogue name from the reference_library index"}},
            "required": ["compound"]},
        handler=read_reference, phase="observe"))

    # -- act (numerical parameter identification) ----------------------- #
    def optimize(args: dict, session) -> ToolResult:
        from ..engines import osp_optimize as OO
        estimate = dict(args.get("estimate") or {})
        # link_scale: fit ONE shared multiplier over a group of collinear params
        # (per-enzyme clearances) so the identifiable TOTAL is fit while the split
        # (their ratio) is held - resolve each member's current value as its base.
        link_scale = []
        cur = None
        for g in (args.get("link_scale") or []):
            members = g.get("members") or []
            if cur is None:
                cur = {p["name"]: p["value"] for p in _current_model(snapshot_path)["parameters"]}
            based = {m: float(cur.get(m, 1.0)) for m in members}
            link_scale.append({"members": based,
                               "bounds": g.get("bounds", [0.1, 10.0]),
                               "name": g.get("name") or ("total:" + "+".join(members))})
        if not estimate and not link_scale:
            return ToolResult.error(
                "provide 'estimate': {parameter: [lo, hi]} - the parameters to fit "
                "numerically (choose 2-4 identifiable, uncertain ones); or "
                "'link_scale' to fit the shared magnitude of a collinear group.")

        # --- the measured-quantity principle (enforced, not just advised) --- #
        #  * a measured physical CONSTANT (MW, pKa, reference pH) can never be
        #    estimated - reject outright.
        #  * a MEASURED-SOFT parameter (fraction unbound, solubility) may be
        #    refined, but only WITHIN its measured uncertainty range - clamp the
        #    requested bounds to that range so the optimizer cannot drag a measured
        #    value outside what the measurement supports to rescue a bad structure.
        lit = (inp.get("given_data", {}) or {}).get("literature_physicochemical", [])
        # if the input GAVE a direct physchem (e.g. lipophilicity), fix it - never re-fit a measured
        # value (the Tizanidine failure: it fit a given logP, freeing the method's identifiability).
        estimate, given_fix, constraint_notes = _fix_given_physchem(estimate, args.get("fix"), lit)
        for name in list(estimate.keys()):
            base = name.split("@", 1)[0]   # qualified per-process name -> base
            tier = osp_catalog.param_tier(base)
            if tier == "constant":
                return ToolResult.error(
                    f"'{name}' is a measured physical constant and cannot be "
                    "estimated (fitting it would model a different molecule). Fix "
                    "it at its literature value and remove it from 'estimate'.")
            if tier == "measured_soft":
                mr = osp_catalog.measured_range(base, lit)
                if mr and isinstance(estimate.get(name), (list, tuple)) \
                        and len(estimate[name]) == 2:
                    lo, hi = float(estimate[name][0]), float(estimate[name][1])
                    clo, chi = max(lo, mr[0]), min(hi, mr[1])
                    if clo >= chi:            # requested band entirely outside
                        clo, chi = mr
                    if not (chi > clo > 0):
                        # the measured range has no positive WIDTH to fit within
                        # (a single measured value) -> you cannot estimate it; FIX it
                        # at the measured value rather than fail with a cryptic bounds
                        # error. (This is what made co-fitting fu with logP fail.)
                        estimate.pop(name)
                        fixval = mr[0] if mr[0] > 0 else (lo if lo > 0 else mr[1])
                        given_fix.setdefault(name, fixval)
                        constraint_notes.append(
                            f"{name}: its measured value has no uncertainty band to fit "
                            f"within (range [{mr[0]:.3g}, {mr[1]:.3g}]) — FIXED at "
                            f"{fixval:.3g} instead of estimated.")
                    elif (clo, chi) != (lo, hi):
                        estimate[name] = [clo, chi]
                        constraint_notes.append(
                            f"{name}: a measured quantity - bounds constrained to "
                            f"its measured range [{mr[0]:.3g}, {mr[1]:.3g}] "
                            f"(requested [{lo:.3g}, {hi:.3g}]).")

        def _progress(i, values, sse, error=None):
            vs = ", ".join(f"{k}={v:.3g}" for k, v in values.items())
            if sse is not None:
                msg = f"       eval {i}: log_sse={sse} [{vs}]"
            else:
                why = f" — {error}" if error else ""
                msg = f"       eval {i}: run FAILED{why} [{vs}]"
            print(msg, flush=True)

        if not estimate and not link_scale:
            return ToolResult.error(
                "after fixing given measurements, no parameter is left to estimate - the input "
                f"already supplies these values ({'; '.join(constraint_notes)}). Choose a different "
                "uncertain parameter (e.g. Intrinsic clearance) to fit.")
        # pre-flight: reject configs PK-Sim can't run, with an actionable reason,
        # instead of burning eval attempts on an opaque 'no result CSVs' crash.
        probs = _structure_problems(args.get("structure"), estimate,
                                    [g.get("members") or [] for g in (args.get("link_scale") or [])])
        if probs:
            return ToolResult.error(
                "this model configuration cannot run (fix the STRUCTURE, not the "
                "bounds): " + " | ".join(probs))

        r = OO.run_optimization(
            cli, snapshot_path, observed, estimate=estimate,
            fix=given_fix, structure=args.get("structure"),
            fit_simulations=args.get("fit_simulations"),
            max_evals=_cap_evals(int(args.get("max_evals") or 30)),
            on_eval=_progress if getattr(config, "stream_optimizer", True) else None,
            link_scale=link_scale or None)
        if not r.get("ok"):
            return ToolResult.error(f"optimization failed: {r.get('message')}")

        plist = [{"parameter": k, "value": v, "unit": ""}
                 for k, v in r["optimized"].items()]
        flags = osp_score.plausibility(plist)
        gmfe = r["fit"].get("gmfe")
        hist = session.get("osp_history") or []
        hist.append({"estimate": list(estimate), "structure": args.get("structure"),
                     "gmfe": gmfe})
        session.put("osp_history", hist)
        best = session.get("osp_best_gmfe")
        if gmfe is not None and (best is None or gmfe < best):
            session.put("osp_best_gmfe", gmfe)
            # store the FIXED parameters alongside the optimized ones - they are
            # part of the model (e.g. GFR fraction=0), so the report's re-run must
            # apply them too or it reproduces a different model than was fitted.
            session.put("osp_best_edits",
                        {"parameters": r["optimized"], "fix": args.get("fix") or {},
                         **(args.get("structure") or {})})
            session.put("osp_best_sensitivity", r.get("sensitivity") or {})
        recs = r.get("recommendations") or []
        # surface the actionable identifiability findings first in the message so
        # the agent acts on them (fix unidentifiable params, refit) rather than
        # re-floating the same parameters next round.
        rec_line = ""
        if recs:
            hi = [a for a in recs if a.get("severity") == "high"]
            lead = hi[0] if hi else recs[0]
            rec_line = (f" | NEXT: {lead['action']}"
                        + (f" (+{len(recs)-1} more identifiability note(s))"
                           if len(recs) > 1 else ""))
        images = _fit_images(r.get("series") or [],
                             f"observed vs simulated — GMFE {gmfe}")
        struct = args.get("structure") or {}
        defaulted = _defaulted_process_params(
            struct.get("add_processes"),
            list(estimate.keys()) + list(given_fix.keys()) + list((args.get("fix") or {}).keys()))
        dflt_note = ""
        if defaulted:
            dflt_note = (" | ⚠ added-process parameter(s) left at their DEFAULT (not fitted "
                         "or set): "
                         + ", ".join(f"{d['parameter']}"
                                     + (f"@{d['molecule']}" if d['molecule'] else "")
                                     + f"={d['default']}" for d in defaulted)
                         + " — a default silently defines the model; fit or set these.")
        out = ToolResult.success(
            f"optimized {list(r['optimized'])} on {len(r['fit_simulations'])} "
            f"study(ies) -> GMFE {gmfe} "
            f"(best so far {session.get('osp_best_gmfe')}){rec_line}{dflt_note}"
            + (" | a fit overlay is attached — inspect the curve shape, not just "
               "the GMFE." if images else ""),
            optimized=r["optimized"], fit=r["fit"], by_route=r["by_route"],
            worst_datasets=r["worst_datasets"],
            params_at_bound=r["params_at_bound"],
            sensitivity=r.get("sensitivity"),
            link_scales=r.get("link_scales"),
            recommendations=recs,
            parameter_flags=flags,
            defaulted_process_params=defaulted or None,   # params silently at their default
            measured_constraints=constraint_notes or None,
            n_evals=r["n_evals"], fit_simulations=r["fit_simulations"],
            iteration=len(hist),
            saw_fit_curve=bool(images))               # transparency marker for the UI
        out.images = images
        return out

    registry.register(Tool(
        name="osp_optimize",
        description=(
            "ACT (numerical fit): you decide WHICH parameters to estimate and "
            "their plausible bounds; a derivative-free optimizer fits them to the "
            "observed data (like PK-Sim Parameter Identification, headless). "
            "estimate={parameter:[lo,hi]} (choose 2-4 identifiable, uncertain "
            "parameters - clearance, permeabilities, effective lipophilicity - "
            "not everything). fix={parameter:value} pins parameters at literature "
            "values. structure={calculation_methods:.., processes:.., "
            "add_processes:[{type,molecule,parameters}]} sets the model structure "
            "(not optimized; see osp_options for legal choices). "
            "link_scale=[{members:[p1,p2,..], bounds:[lo,hi]}] fits ONE shared "
            "multiplier over a collinear group (e.g. per-enzyme clearances the "
            "plasma data cannot split): it recovers their identifiable TOTAL while "
            "holding their RATIO fixed at the members' current values - use it "
            "instead of fixing one member at a guess, which corrupts the total. "
            "Returns the "
            "optimized values, the "
            "full-set GMFE + per-route bias, and params_at_bound (a parameter "
            "pinned to a bound = unidentifiable or wrong structure - reason about "
            "it). Fits against a representative subset of studies for speed. The "
            "result also attaches a SEMILOG PLOT of your predicted curves vs the "
            "observed data (log concentration vs time): LOOK AT IT - a good GMFE "
            "can still hide a missing distribution phase, a wrong terminal slope, "
            "or an over/undershot Cmax or shifted tmax. Diagnose the shape, not "
            "just the number, before your next edit."),
        input_schema={
            "type": "object",
            "properties": {
                "estimate": {"type": "object",
                             "description": "{parameter: [lo, hi]} to fit"},
                "fix": {"type": "object",
                        "description": "{parameter: value} pinned at literature"},
                "structure": {"type": "object",
                              "description": "calculation_methods / processes"},
                "link_scale": {"type": "array", "items": {"type": "object"},
                               "description": "[{members:[p1,p2,..], bounds:[lo,hi]}] "
                               "- fit ONE shared multiplier over a collinear group, "
                               "recovering their total while holding their ratio"},
                "fit_simulations": {"type": "array", "items": {"type": "string"},
                                    "description": "optional exact simulation names"},
                "max_evals": {"type": "integer",
                              "description": "optimizer budget (default 30). Each "
                              "eval is a full PK-Sim build+run (tens of seconds), "
                              "so keep it modest: 20-30 is plenty for 2-4 "
                              "parameters. Prefer a good structure + tight bounds "
                              "over a large budget."},
            },
            "required": ["estimate"],
        },
        handler=optimize, phase="act"))

    def sweep_methods(args, session):
        """DETERMINISTIC structure sweep by COORDINATE DESCENT: sweep the 5 partition methods at an
        anchor permeability, adopt the best, then sweep the 3 permeability methods on that winner
        ONLY when permeability can matter (large molecule, poor partition fit, or caller asked) -
        instead of the exhaustive 5 x 3 = 15 grid, most of which is redundant for perfusion-limited
        drugs. Re-fits the given physchem fresh under each method so a method is never judged at
        frozen (wrong) physchem. full_grid=true forces the exhaustive 15. The DATA picks the method,
        no answer is used."""
        from ..engines import osp_optimize as OO
        from ..engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS

        def _sweep_progress(i, values, sse, error=None):    # live per-eval line during a sweep
            if not getattr(config, "stream_optimizer", True):
                return
            vs = ", ".join(f"{k}={v:.3g}" for k, v in values.items())
            tag = f"log_sse={sse}" if sse is not None else f"run FAILED{(' — ' + error) if error else ''}"
            print(f"       eval {i}: {tag} [{vs}]", flush=True)

        estimate = dict(args.get("estimate") or {})
        if not estimate:
            return ToolResult.error(
                "provide 'estimate': {param:[lo,hi]} - the physchem to RE-FIT under each method "
                "(e.g. Lipophilicity, Fraction unbound). The sweep fits these fresh for every "
                "partition x permeability combo, so a method is never judged at frozen physchem.")
        # A method sweep MUST re-fit the CLEARANCE alongside distribution, or the methods
        # are judged at a wrong default clearance and the correct one (e.g. Rodgers and
        # Rowland) is mis-ranked - the live Triazolam failure, where sweeping only
        # Lipophilicity ranked R&R worst (3.03) and the agent adopted the wrong method.
        _CLEARANCE_KW = ("kcat", "clint", "intrinsic clearance", "specific clearance",
                         "clspec", "plasma clearance", "in vitro vmax", "vmax")
        _is_clear = lambda nm: any(k in (nm or "").lower() for k in _CLEARANCE_KW)
        if not any(_is_clear(n) for n in estimate):
            try:
                model_params = _current_model(snapshot_path)["parameters"]
                clr = [p["name"] for p in model_params if _is_clear(p["name"])]
            except Exception:                       # model unreadable (e.g. unit tests) -> skip guard
                clr = []
            if clr:
                return ToolResult.error(
                    "method ranking needs a CLEARANCE parameter free in `estimate` too. "
                    "With clearance held at a naive default, each distribution method is "
                    "judged at the wrong systemic exposure, so the correct method is "
                    "mis-ranked (e.g. Rodgers and Rowland looked worst on Triazolam only "
                    "because kcat was at its default). Add your clearance parameter to "
                    f"estimate and re-sweep - candidates in this model: {clr[:4]}. Keep a "
                    "GIVEN measured physchem (e.g. measured logP) in `fix`, not in the sweep.")
        fix = args.get("fix") or {}
        parts = args.get("partition_methods") or PARTITION_METHODS
        perms = args.get("permeability_methods") or PERMEABILITY_METHODS
        # Coordinate descent instead of the full 5x3 grid: the permeability method
        # only moves the fit for PERMEABILITY-limited drugs; the OSP library is
        # mostly lipophilic, perfusion-limited small molecules where it is inert.
        # So sweep the partition axis at an anchor permeability first, then sweep
        # the permeability axis ONLY on the winning partition and ONLY when it
        # can matter (large molecule, or the partition fit is still poor, or the
        # caller asked). full_grid=True restores the exhaustive 5x3 for auditing.
        full_grid = bool(args.get("full_grid"))
        perm_requested = bool(args.get("permeability_methods"))
        ANCHOR_PERM = perms[0] if perms else "PK-Sim Standard"
        # Only escalate to the permeability axis when the partition fit is GENUINELY
        # broken (permeability-limited). A merely-mediocre fit (e.g. 1.6-1.9, which
        # just needs more parameter tuning) does NOT indicate a permeability problem
        # - live Triazolam showed every permeability method giving an identical GMFE,
        # so 1.5 was too aggressive and wasted 10 combos. 2.0 still catches a truly
        # permeability-limited drug (e.g. vancomycin fits ~2.35 without the right
        # permeability method). Tunable via perm_threshold.
        PERM_TRIGGER_GMFE = float(args.get("perm_threshold") or 2.0)
        try:
            with open(snapshot_path, encoding="utf-8") as _fh:
                _c0 = (json.load(_fh).get("Compounds") or [{}])[0]
            is_small = _c0.get("IsSmallMolecule", True)
        except (OSError, ValueError):
            is_small = True
        # Carry the WHOLE mechanism the agent passed (add_processes AND processes), not just
        # 'processes' - dropping add_processes ran every method on a no-clearance model, so the
        # clearance estimate had no process to attach to and froze, making the grid meaningless.
        # Only calculation_methods are swept, so strip them from the held-fixed base structure.
        base_structure = dict(args.get("structure") or {})
        base_structure.pop("calculation_methods", None)
        max_evals = _cap_evals(int(args.get("max_evals") or 12))
        # a given physchem (e.g. lipophilicity) must be FIXED, not swept-with-a-free-value: fitting it
        # frees a knob that compensates the method, making the sweep unable to tell methods apart.
        sweep_lit = (inp.get("given_data", {}) or {}).get("literature_physicochemical", [])
        estimate, fix, given_notes = _fix_given_physchem(estimate, fix, sweep_lit)
        for _n in given_notes:
            print(f"  sweep: {_n}", flush=True)
        if not estimate:
            return ToolResult.error(
                "after fixing given measurements there is no free physchem to re-fit under each "
                "method - the method is now identifiable, so just set it with osp_optimize (the "
                "given values are fixed) rather than sweeping.")

        # MEMOIZE: the sweep is expensive (dozens of PK-Sim fits). Choosing the distribution method
        # is a ONE-TIME decision - re-sweeping the SAME grid wastes hours. If an identical sweep
        # (same physchem set + processes) already ran this session, return its result and tell the
        # agent to refine with osp_optimize instead of re-sweeping.
        sig = json.dumps({"e": sorted(estimate), "p": base_structure,
                          "pm": list(parts), "pe": list(perms), "fg": full_grid},
                         default=str, sort_keys=True)
        cache = session.get("osp_sweep_cache") or {}
        if sig in cache:
            c = cache[sig]
            return ToolResult.success(
                f"already swept this grid ({c['partition']} / {c['permeability']} -> GMFE "
                f"{c['gmfe']} was best). The distribution method is chosen; do NOT re-sweep - refine "
                f"with osp_optimize (fit the free physchem / clearance further) or finish.",
                ranked=c.get("ranked"), best=c, cached=True)
        # Rank all 15 method combos in FAST mode (fewer studies, no sensitivity/full-data run) at the
        # FULL eval budget - the speedup comes from fast mode, NOT from cutting evals. Cutting evals
        # too far leaves the free parameter (e.g. CLint) under-fit, so every method looks equally bad
        # and ties - then the ranking is meaningless. Enough evals per combo so the fit reaches a
        # point where the methods actually differ; the winner is then re-fit at full fidelity.
        coarse = max(max_evals, 10)
        lit = (inp.get("given_data", {}) or {}).get("literature_physicochemical", [])

        def _run_pairs(pairs, est, budget, early_stop=False):
            """Coarse-fit a list of explicit (partition, permeability) pairs in FAST
            mode. early_stop applies the COMPENSATION signature (used on the
            single-axis partition phase): >=2 combos tie in GMFE AND a free
            physchem took different values across them, meaning the free param
            absorbed the method difference - a structural degeneracy, safe to stop."""
            jobs = max(1, int(getattr(config, "parallel_jobs", 1) or 1))
            if jobs > 1 and len(pairs) > 1:
                # PARALLEL: each combo's PK-Sim runs in its own temp dir, so the
                # fits are independent. Run them concurrently; prints stay on THIS
                # thread (workers don't print) so stdout capture is not interleaved.
                # early_stop / per-eval streaming are dropped in parallel mode.
                from concurrent.futures import ThreadPoolExecutor
                print(f"  sweep: fitting {len(pairs)} method combos in parallel "
                      f"(jobs={jobs})…", flush=True)

                def _one(pmpe):
                    pm, pe = pmpe
                    st = dict(base_structure)
                    st["calculation_methods"] = {"partition": pm, "permeability": pe}
                    return pm, pe, OO.run_optimization(cli, snapshot_path, observed,
                                                       estimate=est, fix=fix, structure=st,
                                                       max_evals=budget, fast=True)
                out = []
                with ThreadPoolExecutor(max_workers=jobs) as ex:
                    for pm, pe, r in ex.map(_one, list(pairs)):
                        if r.get("ok") and r["fit"].get("gmfe") is not None:
                            out.append({"partition": pm, "permeability": pe,
                                        "gmfe": r["fit"]["gmfe"], "optimized": r["optimized"],
                                        "params_at_bound": r.get("params_at_bound")})
                            fit_str = ", ".join(f"{k}={v:.3g}" for k, v in r["optimized"].items())
                            print(f"  sweep [{pm} / {pe}] -> GMFE {r['fit']['gmfe']}  "
                                  f"(fit: {fit_str})", flush=True)
                return out, False
            out = []
            done = 0
            for pm, pe in pairs:
                structure = dict(base_structure)          # add_processes/processes held fixed
                structure["calculation_methods"] = {"partition": pm, "permeability": pe}
                # announce the combo BEFORE the (multi-minute) fit so the live view
                # shows what is running, not silence; stream per-eval progress too.
                print(f"  sweep [{pm} / {pe}] fitting… ({done + 1}/{len(pairs)})", flush=True)
                r = OO.run_optimization(cli, snapshot_path, observed, estimate=est,
                                        fix=fix, structure=structure, max_evals=budget,
                                        fast=True,   # ranking only: fewer studies, no sens/full
                                        on_eval=_sweep_progress)
                if r.get("ok") and r["fit"].get("gmfe") is not None:
                    out.append({"partition": pm, "permeability": pe,
                                "gmfe": r["fit"]["gmfe"], "optimized": r["optimized"],
                                "params_at_bound": r.get("params_at_bound")})
                    fit_str = ", ".join(f"{k}={v:.3g}" for k, v in r["optimized"].items())
                    print(f"  sweep [{pm} / {pe}] -> GMFE {r['fit']['gmfe']}  (fit: {fit_str})",
                          flush=True)
                done += 1
                if early_stop:
                    gm = [x["gmfe"] for x in out]
                    tie = len(gm) >= 2 and (max(gm) - min(gm)) <= 0.02 * min(gm)
                    compensated = False
                    if done >= 2 and tie:
                        for pname, bnd in est.items():
                            vals = [x["optimized"].get(pname) for x in out
                                    if isinstance(x["optimized"].get(pname), (int, float))]
                            width = abs(bnd[1] - bnd[0]) if isinstance(bnd, (list, tuple)) else 0
                            if len(vals) >= 2 and width and (max(vals) - min(vals)) > 0.05 * width:
                                compensated = True
                                break
                    if compensated:
                        print(f"  sweep: partition method is UNIDENTIFIABLE here - a free physchem "
                              f"compensates it (GMFE ties at ~{min(gm):.3g} while its fitted value "
                              f"shifts across methods); stopping early, the lever is elsewhere",
                              flush=True)
                        return out, True
            return out, False

        def _widen(est, best_row):
            """If the best fit railed on a FREE (non-given, non-measured) parameter, widen that
            parameter to its PHYSICAL range and return the new estimate - so a too-tight self-imposed
            bound (e.g. an effective lipophilicity capped at a measured logP that was never given)
            cannot box the sweep out of the answer. General: only params with NO given measurement
            are widened, and only to the physical plausibility range."""
            widened = dict(est)
            changed = {}
            for ab in best_row.get("params_at_bound") or []:
                pname, side = ab.get("parameter"), ab.get("bound")
                if pname not in est:
                    continue
                base = pname.split("@", 1)[0]
                # respect anything the input GAVE as a measurement, or that has a measured range
                if osp_catalog.measured_range(base, lit) or _given_measurement(base, lit):
                    continue
                pb = osp_score.physical_bounds(base)
                if not pb:
                    continue
                lo, hi = est[pname]
                nlo, nhi = lo, hi
                if side == "upper" and pb[1] > hi:
                    nhi = pb[1]
                if side == "lower" and pb[0] < lo:
                    nlo = pb[0]
                if (nlo, nhi) != (lo, hi):
                    widened[pname] = [nlo, nhi]
                    changed[pname] = [nlo, nhi]
            return (widened, changed) if changed else (None, None)

        if full_grid:
            # exhaustive 5x3 (auditing / permeability-limited compounds you want fully covered)
            pairs_a = [(pm, pe) for pm in parts for pe in perms]
            results, compensated = _run_pairs(pairs_a, estimate, coarse, early_stop=False)
        else:
            # PHASE A: partition axis at the anchor permeability (the axis that actually moves Vd)
            pairs_a = [(pm, ANCHOR_PERM) for pm in parts]
            results, compensated = _run_pairs(pairs_a, estimate, coarse, early_stop=True)
        if not results:
            return ToolResult.error("sweep produced no successful fit across the method grid")
        results.sort(key=lambda x: x["gmfe"])
        est2, changed = _widen(estimate, results[0])
        if est2:
            print(f"  sweep: free non-given param(s) railed -> widening to physical range and "
                  f"re-sweeping: {changed}", flush=True)
            more, comp2 = _run_pairs(pairs_a, est2, coarse, early_stop=not full_grid)
            compensated = compensated or comp2
            results = sorted(results + more, key=lambda x: x["gmfe"])
        # PHASE B: the permeability axis. Coordinate descent (winner partition only) is safe ONLY
        # when the partition fit at the anchor is already GOOD - that means permeability is inert
        # (perfusion-limited) and does not interact with the partition ranking. When the fit is
        # still poor (or it is a large molecule, or the caller asked), permeability CAN interact
        # with partition, so the anchor-based partition ranking is untrustworthy and we must cover
        # the FULL remaining grid (all partitions x non-anchor permeabilities) to not miss the
        # optimum - exactly the Vancomycin (permeability-limited) case.
        if not full_grid and len(perms) > 1 and not compensated:
            best_gmfe = results[0]["gmfe"]
            need_perm = (perm_requested or not is_small
                         or best_gmfe > PERM_TRIGGER_GMFE)
            if need_perm:
                reason = ("caller requested permeability sweep" if perm_requested else
                          "large molecule (permeability-limited)" if not is_small else
                          f"partition fit still poor (GMFE {best_gmfe:.3g} > {PERM_TRIGGER_GMFE}) "
                          f"- axes may interact")
                print(f"  sweep: permeability MATTERS here ({reason}) - covering the full remaining "
                      f"grid so a partition x permeability interaction is not missed", flush=True)
                done = {(x["partition"], x["permeability"]) for x in results}
                pairs_b = [(pm, pe) for pm in parts for pe in perms
                           if pe != ANCHOR_PERM and (pm, pe) not in done]
                more_b, _ = _run_pairs(pairs_b, est2 or estimate, coarse)
                results = sorted(results + more_b, key=lambda x: x["gmfe"])
            else:
                print(f"  sweep: permeability axis SKIPPED (perfusion-limited small molecule and "
                      f"partition fit already good, GMFE {best_gmfe:.3g}); anchored at "
                      f"'{ANCHOR_PERM}'. Pass full_grid=true to force the exhaustive 5x3.",
                      flush=True)
        # FINE refine the winning method at FULL fidelity (all studies + sensitivity + full-data run):
        # the grid ran in fast mode (subset studies, no sensitivity), so re-fit the winner properly
        # for the reported GMFE and identifiability.
        if results:
            top = results[0]
            est_final = est2 or estimate           # est2 (if any) only widens bounds - safe to reuse
            structure = dict(base_structure)          # add_processes/processes held fixed
            structure["calculation_methods"] = {"partition": top["partition"],
                                                 "permeability": top["permeability"]}
            print(f"  sweep: refining winner {top['partition']} / {top['permeability']} "
                  "at full budget…", flush=True)
            rr = OO.run_optimization(cli, snapshot_path, observed, estimate=est_final, fix=fix,
                                     structure=structure, max_evals=max_evals,
                                     on_eval=_sweep_progress)
            if rr.get("ok") and rr["fit"].get("gmfe") is not None:
                top["gmfe"] = rr["fit"]["gmfe"]
                top["optimized"] = rr["optimized"]
                top["params_at_bound"] = rr.get("params_at_bound")
                results.sort(key=lambda x: x["gmfe"])
                print(f"  sweep: refined winner {top['partition']} / {top['permeability']} at full "
                      f"budget -> GMFE {top['gmfe']}", flush=True)
        best = results[0]
        prev = session.get("osp_best_gmfe")
        if prev is None or best["gmfe"] < prev:
            session.put("osp_best_gmfe", best["gmfe"])
            session.put("osp_best_edits",
                        {"parameters": best["optimized"], "fix": fix,
                         "calculation_methods": {"partition": best["partition"],
                                                 "permeability": best["permeability"]}})
        cache[sig] = {**best, "ranked": results}          # remember this grid so it is not re-swept
        session.put("osp_sweep_cache", cache)
        # return a COMPACT ranking (top 5, method + gmfe only) - the full per-combo fit dicts
        # are heavy and not needed by the agent, which only picks the winning method.
        ranked_compact = [{"partition": x["partition"], "permeability": x["permeability"],
                           "gmfe": x["gmfe"]} for x in results[:5]]
        # verdict for the "quick screen -> if none good, adjust and re-sweep" flow.
        # ~GMFE 2.0 is a SOFT target (overridable via perm_threshold), reported as
        # EVIDENCE, not a pass/fail command: below it the method is a good pick; above
        # it the method is likely not the lever - but on a genuinely hard compound this
        # may be its achievable best, so the call is the agent's, not the harness's.
        railed = bool(best.get("params_at_bound"))
        adequate = best["gmfe"] is not None and best["gmfe"] <= PERM_TRIGGER_GMFE and not railed
        advice = None
        if not adequate:
            advice = (f"Best screened GMFE {best['gmfe']} is above the ~{PERM_TRIGGER_GMFE:g} "
                      "soft target"
                      + (", and the winner railed on a bound (so the bound, not the data, is "
                         "setting a value)" if railed else "")
                      + ". That usually means the distribution method is not the remaining "
                      "lever: consider freeing a measured physchem WITHIN its measured range, "
                      "or revisiting the mechanism (add/disable a process), then re-sweeping. "
                      "But if you judge this a genuinely hard compound, this GMFE may be its "
                      "achievable best - your call, not a hard fail.")
            print(f"  sweep: best GMFE {best['gmfe']} above the ~{PERM_TRIGGER_GMFE:g} "
                  f"soft target - {advice}", flush=True)
        msg = (f"screened {len(results)} method combo(s), re-fitting {list(estimate)} under each; "
               f"BEST = {best['partition']} / {best['permeability']} -> GMFE {best['gmfe']} "
               f"(best so far {session.get('osp_best_gmfe')}). ")
        msg += (advice if not adequate else
                "The distribution method is now chosen - refine with osp_optimize, "
                "do NOT re-sweep the same grid.")
        return ToolResult.success(
            msg, ranked_top=ranked_compact, screen_adequate=adequate, advice=advice,
            best={"partition": best["partition"], "permeability": best["permeability"],
                  "gmfe": best["gmfe"], "optimized": best["optimized"]})

    registry.register(Tool(
        name="osp_sweep_methods",
        description=(
            "ACT (deterministic structure sweep by COORDINATE DESCENT): sweep the 5 partition "
            "methods at an anchor permeability, adopt the best, then sweep the 3 permeability "
            "methods on that winner ONLY when it can matter (large molecule, poor partition fit, or "
            "you pass permeability_methods) - typically ~5-8 fits, not the full 5x3=15, since "
            "permeability is inert for perfusion-limited drugs. RE-FITs your physchem "
            "(estimate={param:[lo,hi]}) under each method, so a method is never judged at frozen "
            "physchem. Use whenever distribution / Vd is off. Your `estimate` MUST include a "
            "CLEARANCE parameter (e.g. kcat@CYP3A4 / Intrinsic clearance) alongside the "
            "distribution physchem - a method judged at a wrong default clearance is mis-ranked. "
            "Do NOT distort a measured physchem parameter to fix Vd before you have swept the "
            "methods. Pass the SAME "
            "structure={add_processes:[...], processes:{...}} you use with osp_optimize - the "
            "sweep holds your mechanism fixed and only varies the calculation methods. Optional: "
            "fix={param:value}, partition_methods/permeability_methods to restrict the axes, "
            "max_evals (per combo, default 12), full_grid=true to force the exhaustive 5x3."),
        input_schema={
            "type": "object",
            "properties": {
                "estimate": {"type": "object",
                             "description": "{param:[lo,hi]} physchem to re-fit under each method"},
                "fix": {"type": "object", "description": "{param:value} pinned at literature"},
                "structure": {"type": "object",
                              "description": "{processes:..} to hold the mechanism fixed"},
                "partition_methods": {"type": "array", "items": {"type": "string"}},
                "permeability_methods": {"type": "array", "items": {"type": "string"}},
                "max_evals": {"type": "integer", "description": "optimizer budget per combo (12)"},
                "full_grid": {"type": "boolean",
                              "description": "force the exhaustive 5x3 grid instead of coordinate "
                                             "descent (auditing / permeability-limited drugs)"},
            },
            "required": ["estimate"],
        },
        handler=sweep_methods, phase="act"))

    registry.register(Tool(
        name="osp_try_model",
        description=(
            "ACT / DRY-RUN a model: apply an edit spec, run PK-Sim headless once "
            "(no fitting), and score the fit. edits = {parameters:{name:value}, "
            "calculation_methods:{partition:..,permeability:..}, "
            "processes:{Molecule:true/false}, add_processes:[{type,molecule,"
            "parameters}]}. It APPLIES add_processes too, so use it to VALIDATE a "
            "structure you are about to fit — build it here first with default "
            "parameters; if it runs you get a GMFE, if it can't you get PK-Sim's "
            "real error — before committing to a full osp_optimize. Returns GMFE "
            "and %-within-2-fold overall and PER ROUTE, plus the per-route "
            "geometric BIAS (>1 = over-predicts -> raise clearance; <1 = under). "
            "Use the bias and worst datasets to decide the next edit. It also "
            "attaches a SEMILOG overlay of predicted vs observed curves — read the "
            "shape (distribution phase, terminal slope, Cmax/tmax), not just the "
            "GMFE."),
        input_schema={
            "type": "object",
            "properties": {
                "edits": {
                    "type": "object",
                    "description": "parameters / calculation_methods / processes / add_processes",
                    "properties": {
                        "parameters": {"type": "object"},
                        "calculation_methods": {"type": "object"},
                        "processes": {"type": "object"},
                        "add_processes": {"type": "array", "items": {"type": "object"}},
                    },
                },
            },
            "required": ["edits"],
        },
        handler=try_model, phase="act"))

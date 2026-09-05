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
        comp = (json.load(fh).get("Compounds") or [{}])[0]
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
        return ToolResult.success(
            "task description: objective, known biology, literature priors, "
            "current model, and the observed clinical data",
            objective=inp.get("objective"),
            background=bg.get("description"),
            known_biology=bg.get("literature_facts"),
            compound_identity=gd.get("compound_identity"),
            literature_physicochemical=gd.get("literature_physicochemical")
            or inp.get("literature_physicochemical"),
            unknowns_guidance=inp.get("unknowns_guidance"),
            parameters_to_determine=to_determine,
            given_input_parameters=given_in_model,
            value_status_note=(
                "current_model values carry a 'value_status': 'given' = a "
                "measured/published input (trust it); 'placeholder' = a naive "
                "benchmark default that is NOT a measurement - determine it from "
                "the data or from established physchem knowledge, do not trust the "
                "number shown. The starting value of a placeholder carries no "
                "information about the answer."),
            evaluation_rubric=inp.get("evaluation_rubric"),
            current_model=model,
            observed_overview=_observed_overview(observed),
            valid_partition_methods=PARTITION_METHODS,
            valid_permeability_methods=PERMEABILITY_METHODS,
        )

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
        return ToolResult.success(
            f"GMFE {overall.get('gmfe')} overall "
            f"(within2fold {overall.get('within_2fold_pct')}%); "
            f"best so far {session.get('osp_best_gmfe')}",
            gmfe_overall=overall.get("gmfe"),
            within_2fold_pct=overall.get("within_2fold_pct"),
            bias_overall=overall.get("bias"),
            by_route=score["by_route"],
            worst_datasets=worst,
            parameter_flags=flags,
            edits_applied=applied,
            not_found=applied.get("not_found"),
            n_matched=len(predicted), n_total=len(observed),
            best_gmfe_so_far=session.get("osp_best_gmfe"),
            iteration=len(hist),
        )

    # -- observe (authoritative action space) --------------------------- #
    def options(args: dict, session) -> ToolResult:
        model = _current_model(snapshot_path)
        expressed = _expressed_molecules(snapshot_path)
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
            if "@" in p["name"]:
                entry["note"] = ("per-process parameter - estimate/set it with this "
                                 "exact qualified name to target only this process")
            # the CURRENT value's trust status: a placeholder is a naive default,
            # NOT a measurement - the shown number carries no information.
            if p.get("value_status") == "placeholder":
                entry["value_note"] = ("the shown value is a PLACEHOLDER default "
                                       "(not measured); DETERMINE this parameter - "
                                       "do not trust the starting number. If you "
                                       "have a literature/textbook value, fix it "
                                       "there; otherwise estimate it from the data.")
            elif p.get("value_status") == "given":
                entry["value_note"] = "the shown value is a measured/given input - trust it"
            if tier == "constant":
                entry["rule"] = "measured constant - cannot be estimated; fix it"
            elif tier == "measured_soft":
                mr = osp_catalog.measured_range(p["name"], lit)
                entry["measured_range"] = list(mr) if mr else None
                entry["rule"] = ("measured - fix by default; may be estimated only "
                                 "when justified and only within measured_range")
            editable.append(entry)
        return ToolResult.success(
            "authoritative action space: what you may edit and legal choices",
            molecule_type=molecule_type,
            identification_principle=(
                "Measured physical constants (tier=constant: MW, pKa, reference "
                "pH) are never estimated. Measured-but-refinable quantities "
                "(tier=measured_soft: fraction unbound, solubility) are fixed by "
                "default and may be estimated only when a residual misfit justifies "
                "it, constrained to their measured_range. Estimate the minimal "
                "identifiable set first (tier=estimate); free a measured quantity "
                "only if the data demand it."),
            editable_parameters=editable,
            calculation_methods={
                "partition": {"current": next(
                    (m for m in model["calculation_methods"]
                     if "partition" in m.lower()), None),
                    "options": osp_catalog.PARTITION_METHOD_INFO},
                "permeability": {"current": next(
                    (m for m in model["calculation_methods"]
                     if "permeability" in m.lower()), None),
                    "options": osp_catalog.PERMEABILITY_METHOD_INFO}},
            processes_present=model["processes"],
            expressed_molecules=expressed,
            addable_process_types=osp_catalog.addable_process_types(expressed),
            interaction_process_types=osp_catalog.interaction_process_types(),
            process_catalog_notes=(
                "addable_process_types are single-compound mechanisms you can add "
                "via add_processes. interaction_process_types (DDI: inhibition / "
                "induction) are part of PK-Sim's library but are NOT addable here - "
                "they link a perpetrator to a victim's enzyme and require a multi-"
                "compound DDI setup. Each entry has 'validated' (confirmed to run "
                "through PKSim.CLI) and 'provenance'; prefer validated mechanisms, "
                "and treat validated=false as needing a confirmation run."),
            edit_spec_help={
                "parameters": "{name: value} - names above",
                "calculation_methods": "{partition: <opt>, permeability: <opt>}",
                "processes": "{molecule: false} to disable an existing process",
                "add_processes": "[{type, molecule, parameters}] - attach a NEW "
                "mechanism (enzyme process needs an expressed enzyme)"},
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
        constraint_notes = []
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
                    if (clo, chi) != (lo, hi):
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

        r = OO.run_optimization(
            cli, snapshot_path, observed, estimate=estimate,
            fix=args.get("fix"), structure=args.get("structure"),
            fit_simulations=args.get("fit_simulations"),
            max_evals=int(args.get("max_evals") or 30),
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
        return ToolResult.success(
            f"optimized {list(r['optimized'])} on {len(r['fit_simulations'])} "
            f"study(ies) -> GMFE {gmfe} "
            f"(best so far {session.get('osp_best_gmfe')}){rec_line}",
            optimized=r["optimized"], fit=r["fit"], by_route=r["by_route"],
            worst_datasets=r["worst_datasets"],
            params_at_bound=r["params_at_bound"],
            sensitivity=r.get("sensitivity"),
            link_scales=r.get("link_scales"),
            recommendations=recs,
            sensitivity_hint=("'relative' is each parameter's local influence on "
                              "the fit, normalised to the most influential (1.0). "
                              "A parameter with a small 'relative' is weakly "
                              "constrained by the data - its fitted value is "
                              "uncertain. ACT on 'recommendations': fix a weakly-"
                              "identified or collinear parameter to a known/"
                              "literature value and refit the identifiable set - "
                              "do not keep floating a parameter the data cannot "
                              "pin (it lands on an artifact)."),
            parameter_flags=flags,
            measured_constraints=constraint_notes or None,
            n_evals=r["n_evals"], fit_simulations=r["fit_simulations"],
            iteration=len(hist))

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
            "it). Fits against a representative subset of studies for speed."),
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
        """DETERMINISTIC structure sweep: try EVERY partition x permeability method and RE-FIT the
        given physchem under each, then adopt the best by GMFE. This is the enumerable grid (5 x 3 =
        15) that should be brute-forced, not sampled by the LLM: it removes the failure where the
        agent tries a method at frozen (wrong) physchem, sees no improvement, and wrongly concludes
        'distribution insensitive'. General: the DATA picks the method, no answer is used."""
        from ..engines import osp_optimize as OO
        from ..engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
        estimate = dict(args.get("estimate") or {})
        if not estimate:
            return ToolResult.error(
                "provide 'estimate': {param:[lo,hi]} - the physchem to RE-FIT under each method "
                "(e.g. Lipophilicity, Fraction unbound). The sweep fits these fresh for every "
                "partition x permeability combo, so a method is never judged at frozen physchem.")
        fix = args.get("fix") or {}
        parts = args.get("partition_methods") or PARTITION_METHODS
        perms = args.get("permeability_methods") or PERMEABILITY_METHODS
        processes = (args.get("structure") or {}).get("processes")
        max_evals = int(args.get("max_evals") or 12)
        lit = (inp.get("given_data", {}) or {}).get("literature_physicochemical", [])

        def _grid(est):
            out = []
            for pm in parts:
                for pe in perms:
                    structure = {"calculation_methods": {"partition": pm, "permeability": pe}}
                    if processes:
                        structure["processes"] = processes
                    r = OO.run_optimization(cli, snapshot_path, observed, estimate=est,
                                            fix=fix, structure=structure, max_evals=max_evals)
                    if r.get("ok") and r["fit"].get("gmfe") is not None:
                        out.append({"partition": pm, "permeability": pe,
                                    "gmfe": r["fit"]["gmfe"], "optimized": r["optimized"],
                                    "params_at_bound": r.get("params_at_bound")})
                        print(f"  sweep [{pm} / {pe}] -> GMFE {r['fit']['gmfe']}", flush=True)
            return out

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

        results = _grid(estimate)
        if not results:
            return ToolResult.error("sweep produced no successful fit across the method grid")
        results.sort(key=lambda x: x["gmfe"])
        est2, changed = _widen(estimate, results[0])
        if est2:
            print(f"  sweep: free non-given param(s) railed -> widening to physical range and "
                  f"re-sweeping: {changed}", flush=True)
            more = _grid(est2)
            results = sorted(results + more, key=lambda x: x["gmfe"])
        best = results[0]
        prev = session.get("osp_best_gmfe")
        if prev is None or best["gmfe"] < prev:
            session.put("osp_best_gmfe", best["gmfe"])
            session.put("osp_best_edits",
                        {"parameters": best["optimized"], "fix": fix,
                         "calculation_methods": {"partition": best["partition"],
                                                 "permeability": best["permeability"]}})
        return ToolResult.success(
            f"swept {len(results)} method combos, re-fitting {list(estimate)} under each; BEST = "
            f"{best['partition']} / {best['permeability']} -> GMFE {best['gmfe']} "
            f"(best so far {session.get('osp_best_gmfe')}). Adopt the best method, then refine.",
            ranked=results, best=best)

    registry.register(Tool(
        name="osp_sweep_methods",
        description=(
            "ACT (deterministic structure sweep): try EVERY partition x permeability calculation "
            "method (the full 5x3 grid) and RE-FIT your physchem (estimate={param:[lo,hi]}) under "
            "each, then adopt the best-GMFE method. Use this whenever distribution / Vd is off - it "
            "is cheaper and more reliable than guessing methods one at a time, and it avoids the "
            "trap of judging a method at frozen physchem. Do NOT distort a measured physchem "
            "parameter to fix Vd before you have swept the methods. Optional: fix={param:value}, "
            "structure={processes:..} to hold your mechanism fixed, partition_methods/"
            "permeability_methods to restrict the grid, max_evals (per combo, default 12)."),
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
            },
            "required": ["estimate"],
        },
        handler=sweep_methods, phase="act"))

    registry.register(Tool(
        name="osp_try_model",
        description=(
            "ACT: apply an edit spec to the model, run PK-Sim headless, and score "
            "the fit against the observed data. edits = {parameters:{name:value}, "
            "calculation_methods:{partition:..,permeability:..}, "
            "processes:{Molecule:true/false}}. Returns GMFE and %-within-2-fold "
            "overall and PER ROUTE, plus the per-route geometric BIAS (>1 = the "
            "model over-predicts -> lower exposure/raise clearance; <1 = "
            "under-predicts). Use the bias and worst datasets to decide the next "
            "edit. Iterate to minimize GMFE."),
        input_schema={
            "type": "object",
            "properties": {
                "edits": {
                    "type": "object",
                    "description": "parameters / calculation_methods / processes",
                    "properties": {
                        "parameters": {"type": "object"},
                        "calculation_methods": {"type": "object"},
                        "processes": {"type": "object"},
                    },
                },
            },
            "required": ["edits"],
        },
        handler=try_model, phase="act"))

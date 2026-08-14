"""Generate a scientific evaluation report from an agent modeling run.

Mirrors an OSP evaluation report, auto-filled from the run:

  Objective / Data / Model / Model choice / Parameters / Pharmacological
  rationale / Concentration-time analysis / Comparison with the ground-truth
  model / Conclusion / the full LLM trajectory.

Two renderers share one assembled ``ReportData``:
  * ``write_html`` - self-contained HTML with inline SVG concentration-time
    plots; pure standard library, works anywhere (open it and Print -> PDF).
  * ``write_pdf``  - a real multi-page PDF via matplotlib (if installed).

The narrative sections (model-choice rationale, why the parameters are
pharmacologically sound, conclusion) are written by an optional Claude synthesis
pass; without an API key they fall back to a deterministic write-up built from
the parameters, their plausible ranges/priors, and the fit.
"""

from __future__ import annotations

import html
import math
from dataclasses import dataclass, field
from typing import Any

from .engines import osp_catalog
from .engines import osp_score


def _nca(o: dict) -> dict[str, Any]:
    pts = sorted((t, c) for t, c in zip(o["time_h"], o["conc_mg_L"])
                 if isinstance(t, (int, float)) and isinstance(c, (int, float)))
    if not pts:
        return {"study": o.get("study") or o["dataset"], "n": 0}
    t = [p[0] for p in pts]; c = [p[1] for p in pts]
    cmax = max(c); tmax = t[c.index(cmax)]
    auc = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in range(1, len(t)))
    tail = [(t[i], c[i]) for i in range(len(t)) if c[i] > 0][-3:]
    thalf = None
    if len(tail) >= 2 and tail[0][1] > tail[-1][1] and tail[-1][0] > tail[0][0]:
        k = (math.log(tail[0][1]) - math.log(tail[-1][1])) / (tail[-1][0] - tail[0][0])
        thalf = round(math.log(2) / k, 2) if k > 0 else None
    return {"study": (o.get("study") or o["dataset"])[:40], "c_max_mg_L": round(cmax, 4),
            "t_max_h": round(tmax, 2), "auc_mg_h_L": round(auc, 3), "t_half_h": thalf}


def _trajectory(session) -> list[dict[str, Any]]:
    from .state import Decision, Observation
    steps, pending = [], None
    for ev in session.transcript:
        if isinstance(ev, Decision):
            acts = ", ".join(c.name for c in ev.calls) or "finish"
            pending = {"reason": ev.text or "", "action": acts, "result": ""}
            steps.append(pending)
        elif isinstance(ev, Observation) and pending is not None:
            c = ev.content
            bits = [c.get("message", "")]
            if c.get("optimized"):
                bits.append(f"optimized={c['optimized']}")
            if c.get("params_at_bound"):
                bits.append(f"at_bound={c['params_at_bound']}")
            pending["result"] = (pending["result"] + " | " if pending["result"] else "") \
                + " ".join(str(b) for b in bits if b)
    return steps


def _diagnostics(session, best_edits, fit) -> dict[str, Any]:
    """Honest self-assessment of the run, so the report cannot silently present
    an un-fitted or failed model as a finished result. Reads the transcript for
    what actually happened (did an optimize succeed? did the engine error?)."""
    from .state import Observation
    opt_ok = False
    engine_errors = 0
    for ev in getattr(session, "transcript", []):
        if not isinstance(ev, Observation):
            continue
        if not ev.ok:
            engine_errors += 1
            continue
        if getattr(ev, "tool", None) == "osp_optimize" and \
                isinstance(ev.content, dict) and ev.content.get("optimized"):
            opt_ok = True
    params_fitted = bool((best_edits or {}).get("parameters"))
    g = fit.get("gmfe")
    if g is None:
        verdict = "no runnable model was produced"
    elif g <= 1.5:
        verdict = "good"
    elif g <= 2.0:
        verdict = "acceptable"
    else:
        verdict = "poor"
    return {"optimization_succeeded": opt_ok, "engine_errors": engine_errors,
            "params_fitted": params_fitted, "gmfe": g, "fit_verdict": verdict}


def _status_banner(diag: dict) -> str:
    """A one-line, prominent status shown at the top of the report."""
    g = diag.get("gmfe")
    if g is None:
        return ("⚠ NO RUNNABLE MODEL — the loop did not produce a scored model; "
                "results below are incomplete.")
    if not diag.get("optimization_succeeded"):
        why = (f" ({diag['engine_errors']} engine error(s) during the run)"
               if diag.get("engine_errors") else "")
        return ("⚠ NOT OPTIMIZED — the numerical optimizer did not complete, so "
                "parameters below were NOT fitted to the data; the reported "
                f"GMFE {g} reflects an un-fitted / hand-set model{why}.")
    if diag.get("fit_verdict") == "poor":
        return (f"⚠ POOR FIT — overall GMFE {g} (>2-fold typical error). The "
                "model is not yet an adequate description of the data.")
    return (f"Fit is {diag['fit_verdict']} (overall GMFE {g}); parameters were "
            "fitted by the numerical optimizer.")


def _ode_section(structure: dict, parameters: list, fixed: dict | None) -> dict:
    """The governing PBPK equations, annotated for THIS model.

    PK-Sim instantiates the full system (~15 perfused organs, each split into
    plasma / blood-cell / interstitial / intracellular sub-compartments) at build
    time - hundreds of formula nodes. What belongs in an evaluation report is the
    governing mass balance plus the elimination/absorption terms that are
    actually active here, and where each fitted parameter enters."""
    procs = [str(p) for p in (structure.get("processes") or [])]
    fixed = fixed or {}
    has_metab = any("CYP" in p or "Metaboli" in p or "AADAC" in p for p in procs)
    # renal active only if a GFR/renal process is present AND GFR fraction != 0
    gfr_val = fixed.get("GFR fraction")
    renal_present = any("Glomerul" in p or "GFR" in p or "Kidney" in p for p in procs)
    renal_active = renal_present and (gfr_val is None or float(gfr_val or 0) > 0)

    fitted = {p["name"] for p in parameters}
    perm_limited = any("permeab" in n.lower() and "intestinal" not in n.lower()
                       for n in fitted) or "Permeability" in fitted

    eqs = [
        ("Perfused tissue (distribution)",
         "V_t · dC_t/dt = Q_t · ( C_art − C_t / (K_p,t / R_bp) )",
         "Each non-eliminating organ t (brain, heart, muscle, skin, adipose, "
         "bone, ...): rate in from arterial blood at flow Q_t, rate out with the "
         "venous effluent in equilibrium with tissue via the tissue:plasma "
         "partition coefficient K_p,t (R_bp = blood:plasma ratio)."),
        ("Venous blood pool",
         "V_ven · dC_ven/dt = Σ_t Q_t · C_t / (K_p,t / R_bp) − Q_co · C_ven + R_iv(t)",
         "Collects the effluent of every organ; Q_co = cardiac output. IV doses "
         "enter here as an input rate R_iv(t) (bolus or infusion)."),
        ("Arterial blood pool",
         "V_art · dC_art/dt = Q_co · ( C_lung − C_art )",
         "Blood is oxygenated/mixed through the lung and distributed to the organs."),
    ]
    if perm_limited:
        eqs.append((
            "Permeability-limited exchange (cellular)",
            "V_cell · dC_cell/dt = P · SA · ( f_u·C_int − f_u,cell·C_cell / K_p,t )",
            "Where cellular permeability limits uptake, the intracellular space is "
            "a separate state exchanging with interstitium across the membrane "
            "(permeability P × surface area SA), driven by the unbound gradient."))
    if has_metab:
        eqs.append((
            "Liver (hepatic metabolism)",
            "V_liv · dC_liv/dt = Q_liv·C_art + Q_po·C_po − Q_hv·C_liv/(K_p,liv/R_bp) "
            "− CL_int · f_u · C_liv",
            "The liver receives arterial + portal (Q_po, first-pass) inflow; "
            "CYP3A4 metabolism removes drug as a first-order intrinsic-clearance "
            "sink CL_int·f_u·C_liv on the unbound concentration."))
    eqs.append((
        "Kidney (renal filtration)",
        "V_kid · dC_kid/dt = Q_kid·(C_art − C_kid/(K_p,kid/R_bp)) − GFR · f_u · C_kid",
        "Glomerular filtration removes unbound drug at rate GFR·f_u·C_kid."
        + ("" if renal_active else "  In THIS model this term is ZERO (GFR "
           "fraction = 0): renal excretion of unchanged drug is negligible, so "
           "the renal sink is switched off.")))
    eqs.append((
        "Oral absorption (PO doses)",
        "dA_lum/dt = −P_int · SA_int · C_lum,u ;   J_abs = P_int · SA_int · C_lum,u  →  gut wall → portal vein → liver",
        "Dissolved drug in the intestinal lumen permeates the mucosa "
        "(specific intestinal permeability P_int × surface area), enters the gut "
        "wall, and drains via the portal vein to the liver (first-pass)."))

    # tissue:plasma partition coefficient definition (the chosen method)
    methods = structure.get("calculation_methods") or []
    part = next((m.split(" - ")[-1] for m in methods if "partition" in m.lower()),
                "the selected method")
    kp = ("K_p,t is not a fitted number but is computed per organ from tissue "
          f"composition by the {part} method — a function of effective "
          "lipophilicity, the unbound fraction f_u, drug ionization (pKa) and each "
          "tissue's lipid / water / phospholipid / protein content.")

    # where each fitted parameter enters the ODEs
    where = {
        "Lipophilicity": "sets every tissue partition coefficient K_p,t (the "
                         "distribution terms in all organ ODEs).",
        "Intrinsic clearance": "the hepatic elimination sink CL_int in the liver ODE.",
        "Permeability": "the cellular permeability P in the permeability-limited "
                        "exchange term.",
        "Specific intestinal permeability (transcellular)":
            "P_int in the oral-absorption flux J_abs.",
        "Fraction unbound (plasma, reference value)":
            "f_u — scales every unbound-driven term (K_p,t, CL_int, GFR).",
        "GFR fraction": "the renal filtration sink GFR in the kidney ODE.",
    }
    param_map = []
    for p in parameters:
        w = where.get(p["name"])
        if w:
            param_map.append((p["name"], w))
    for name in (fixed or {}):
        if name in where and name not in {p["name"] for p in parameters}:
            param_map.append((f"{name} (fixed)", where[name]))

    compartments = (
        "PK-Sim's standard human whole-body structure: ~15 perfused organs "
        "(lung, brain, heart, kidney, gut, spleen, pancreas, liver, stomach, "
        "muscle, skin, bone, adipose, gonads) plus arterial and venous blood "
        "pools. Each organ is further divided into plasma, blood-cell, "
        "interstitial and intracellular sub-compartments, so the solved system is "
        "one mass-balance ODE per sub-compartment (order ~10^2 states). The "
        "physiological volumes V and flows Q are fixed by the individual; only the "
        "drug terms below are model choices.")
    caveat = (
        "These are the governing equations; PK-Sim generates the fully "
        "instantiated system (every organ and sub-compartment, and the exact "
        f"{part} K_p formula) at build time and integrates it numerically. The "
        "concentrations plotted are the peripheral-venous-plasma state.")
    return {"compartments": compartments, "equations": eqs, "kp": kp,
            "param_map": param_map, "caveat": caveat}


def assemble(session, config, cli, input_dict, snapshot_path, best_edits,
             ref_snapshot_path=None, answer_edits=None, run_models=True):
    """Build ReportData from a finished run (re-runs the best + reference models
    via PK-Sim to get the concentration-time profiles)."""
    from .engines.snapshot_edit import apply_edits
    observed = input_dict["given_data"]["clinical_observed_data"]
    best_edits = best_edits or {}
    # the model actually fitted = optimized params + the FIXED params (e.g. GFR=0)
    # + structure. Merge fix into parameters so the re-run reproduces the fit.
    fixed = dict(best_edits.get("fix") or {})
    run_edits = {k: v for k, v in best_edits.items() if k != "fix"}
    run_edits["parameters"] = {**fixed, **(best_edits.get("parameters") or {})}

    fit, ref, profiles = {}, {}, []
    if run_models and cli is not None:
        res = cli.build_and_run(snapshot_path, edits=run_edits)
        pred, _ = osp_score.map_predictions(res.get("profiles", []), observed)
        score = osp_score.score_fit(observed, pred)
        fit = dict(score["overall"]); fit["by_route"] = score["by_route"]
        predmap = {p["dataset"]: p for p in pred}
        refmap = {}
        if ref_snapshot_path:
            rres = cli.build_and_run(ref_snapshot_path)
            rpred, _ = osp_score.map_predictions(rres.get("profiles", []), observed)
            ref["gmfe"] = osp_score.score_fit(observed, rpred)["overall"]["gmfe"]
            refmap = {p["dataset"]: p for p in rpred}
        for o in observed:
            pm, rm = predmap.get(o["dataset"]), refmap.get(o["dataset"])
            profiles.append({
                "study": o["dataset"],
                "obs": list(zip(o["time_h"], o["conc_mg_L"])),
                "pred": list(zip(pm["time_h"], pm["pred_conc_mg_L"])) if pm else [],
                "ref": list(zip(rm["time_h"], rm["pred_conc_mg_L"])) if rm else None})

    # final structure = snapshot + best edits
    import json as _json
    with open(snapshot_path, encoding="utf-8") as _fh:
        _snap0 = _json.load(_fh)
    final, _ = apply_edits(_snap0, run_edits)
    comp = (final.get("Compounds") or [{}])[0]
    structure = {
        "calculation_methods": comp.get("CalculationMethods") or [],
        "processes": [p.get("Molecule") or p.get("InternalName")
                      for p in comp.get("Processes") or []]}

    ref_params = (answer_edits or {}).get("parameters", {})
    estimated = best_edits.get("parameters") or {}
    params = []
    for name, val in estimated.items():
        cat = osp_catalog.describe_parameter(name)
        # role reflects what ACTUALLY happened this run: everything in the
        # optimized set was estimated, even if the catalog's default role is
        # "measured" (e.g. fraction unbound moved into the fit). Do not let the
        # static catalog label contradict the trajectory.
        params.append({"name": name, "value": val, "unit": cat.get("unit", ""),
                       "role": "estimated",
                       "plausible_range": cat.get("range"),
                       "reference": ref_params.get(name)})
    # fixed parameters (held at literature values, e.g. GFR fraction = 0) shown
    # as their own rows so the table is complete and unambiguous.
    fixed_rows = []
    for name, val in fixed.items():
        if name in estimated:
            continue
        cat = osp_catalog.describe_parameter(name)
        fixed_rows.append({"name": name, "value": val, "unit": cat.get("unit", ""),
                           "role": "fixed",
                           "plausible_range": cat.get("range"),
                           "reference": ref_params.get(name)})

    diag = _diagnostics(session, best_edits, fit)
    bg = input_dict.get("background") or {}
    d = ReportData(
        title=f"PBPK evaluation report — {input_dict.get('compound','compound')}",
        objective=input_dict.get("objective", ""),
        background=bg.get("description", ""),
        known_biology=bg.get("literature_facts", []),
        data_overview={"n_datasets": len(observed),
                       "routes": sorted({o.get("route") for o in observed if o.get("route")})},
        nca_rows=[_nca(o) for o in observed],
        structure=structure, parameters=params + fixed_rows, fit=fit, reference=ref,
        profiles=profiles, narrative={}, trajectory=_trajectory(session),
        diagnostics=diag, status=_status_banner(diag),
        odes=_ode_section(structure, params, fixed),
        literature=(input_dict.get("given_data", {}) or {}).get(
            "literature_physicochemical", []) or [])
    d.narrative = (llm_narrative(d, config) if config and config.anthropic_key_present()
                   else deterministic_narrative(d))
    return d


# --------------------------------------------------------------------------- #
# assembled report data
# --------------------------------------------------------------------------- #

@dataclass
class ReportData:
    title: str
    objective: str
    background: str
    known_biology: list[str]
    data_overview: dict[str, Any]
    nca_rows: list[dict[str, Any]]
    structure: dict[str, Any]              # methods, processes
    parameters: list[dict[str, Any]]       # name, value, unit, role, range, prior
    fit: dict[str, Any]                    # gmfe, within2fold, by_route
    reference: dict[str, Any]              # reference params + gmfe (ground truth)
    profiles: list[dict[str, Any]]         # per study: observed + predicted (+ref)
    narrative: dict[str, str]              # model_choice, parameter_rationale, conclusion
    trajectory: list[dict[str, Any]]       # steps: reason, action, result
    diagnostics: dict[str, Any] = field(default_factory=dict)  # honest run self-assessment
    status: str = ""                       # one-line banner shown at the top
    odes: dict[str, Any] = field(default_factory=dict)  # governing equations + param map
    literature: list = field(default_factory=list)  # literature physchem anchors


# --------------------------------------------------------------------------- #
# deterministic narrative (fallback when no LLM synthesis)
# --------------------------------------------------------------------------- #

def _lit_anchor(param_name: str, literature: list) -> tuple[float, float] | None:
    """Map a model parameter to the literature physchem range [lo, hi] if there
    is a comparable measured value, so the fitted value can be checked against it."""
    key = param_name.lower()
    def vals(entry):
        out = []
        if entry.get("value") is not None:
            out.append(float(entry["value"]))
        out += [float(x) for x in entry.get("reported_values", [])]
        rp = entry.get("reported_range_percent")
        if rp:
            out += [float(x) / 100.0 for x in rp]
        return out
    for e in literature or []:
        pn = str(e.get("parameter", "")).lower()
        if "lipophil" in key and ("logd" in pn or "logp" in pn):
            v = vals(e)
            if v:
                return (min(v), max(v))
        if "unbound" in key and "unbound" in pn:
            v = vals(e)
            if v:
                return (min(v), max(v))
    return None


def deterministic_narrative(d: "ReportData") -> dict[str, str]:
    methods = d.structure.get("calculation_methods", [])
    procs = d.structure.get("processes", [])
    mc = (f"Distribution was described with {', '.join(methods) or 'the default methods'}. "
          f"The active processes were {', '.join(procs) or 'none'}. Parameters that "
          "cannot be reliably transferred from in-vitro to in-vivo were estimated "
          "against the clinical data; measured physicochemical inputs were fixed.")
    lines = []
    for p in d.parameters:
        rng = p.get("plausible_range")
        v = p.get("value")
        if p.get("role") == "fixed":
            continue   # rationale is about the estimated parameters
        msg = f"- {p['name']} = {v:.4g} {p.get('unit','')}".rstrip()
        # compare against the independent literature anchor when one exists
        anchor = _lit_anchor(p["name"], d.literature) if v is not None else None
        if anchor:
            lo, hi = anchor
            if lo <= v <= hi:
                msg += f" — consistent with the literature value ({lo:g}–{hi:g})."
            else:
                fold = v / ((lo + hi) / 2) if (lo + hi) else None
                near = lo if v < lo else hi
                msg += (f" — DEPARTS from the literature value ({lo:g}–{hi:g}): "
                        f"fitted {'below' if v < near else 'above'} the measured "
                        f"value{f' (~{fold:.2g}× the midpoint)' if fold else ''}; "
                        "likely effective-vs-measured difference or optimizer "
                        "compensation — interpret with caution.")
        elif "permeab" in p["name"].lower():
            msg += (" — no independent measurement; identified only from plasma "
                    "data and typically WEAKLY IDENTIFIABLE, so its absolute value "
                    "is uncertain.")
        elif rng and v is not None and rng[0] <= v <= rng[1]:
            msg += f" — within the physiological range {rng}."
        elif rng:
            msg += f" — OUTSIDE the expected range {rng}; review."
        lines.append(msg)
    pr = ("Each estimated parameter, checked against independent literature where "
          "available:\n" + "\n".join(lines))
    # blind conclusion: judged on the model's own fit + plausibility, not the
    # reference (the ground-truth comparison is a separate, factual section).
    # Tone MUST match the fit - do not call a poor fit adequate.
    g = d.fit.get("gmfe")
    diag = d.diagnostics or {}
    verdict = diag.get("fit_verdict")
    if not diag.get("optimization_succeeded", True):
        concl = (f"The numerical optimizer did not complete for this run, so the "
                 f"parameters were not fitted to the data. The resulting model "
                 f"gives an overall GMFE of {g}, which should be read as a "
                 "starting point, not a fitted result. Re-running the parameter "
                 "identification is required before drawing PK conclusions.")
    elif verdict == "poor":
        concl = (f"The model does NOT yet adequately reproduce the observed plasma "
                 f"concentrations: overall GMFE is {g} (beyond the ~2-fold typical "
                 "error). This indicates a structural or parameter problem to "
                 "resolve (check per-route bias and any parameters pinned to a "
                 "bound) before the model can be considered fit for purpose.")
    else:
        concl = (f"The model reproduces the observed plasma concentrations with an "
                 f"overall GMFE of {g} ({verdict}). The fit and the parameter "
                 "values are consistent with the drug's known disposition, "
                 "supporting the model as an adequate description of the "
                 "pharmacokinetics.")
    return {"model_choice": mc, "parameter_rationale": pr, "conclusion": concl}


def llm_narrative(d: "ReportData", config) -> dict[str, str]:
    """Optional Claude synthesis of the scientific narrative sections."""
    try:
        import anthropic
    except Exception:
        return deterministic_narrative(d)
    # BLIND: the narrative is written without the ground-truth model, so the
    # scientific rationale is justified on physics + data, not on "it matches the
    # reference". The reference stays only in the factual comparison section.
    blind_params = [{k: v for k, v in p.items() if k != "reference"}
                    for p in d.parameters]
    diag = d.diagnostics or {}
    payload = {
        "objective": d.objective, "known_biology": d.known_biology,
        "structure": d.structure, "parameters": blind_params,
        # the independent LITERATURE anchors (measured physchem) so the rationale
        # can COMPARE each fitted value to what is known, not just rubber-stamp it
        # as "within range". This is what turns the section from narration into a
        # critique. (These are public literature values, NOT the fitted reference.)
        "literature_physicochemical": d.literature,
        "fit": {k: v for k, v in d.fit.items() if k != "reference"},
        # honest facts the narrative MUST respect - do not describe fitting that
        # did not happen, and match the tone to the fit quality.
        "run_facts": {
            "optimization_completed": diag.get("optimization_succeeded", True),
            "parameters_were_fitted": diag.get("params_fitted", True),
            "engine_errors": diag.get("engine_errors", 0),
            "fit_verdict": diag.get("fit_verdict")},
    }
    tool = {"name": "write_sections", "description": "Write the report narrative.",
            "input_schema": {"type": "object", "properties": {
                "model_choice": {"type": "string"},
                "parameter_rationale": {"type": "string"},
                "conclusion": {"type": "string"}},
                "required": ["model_choice", "parameter_rationale", "conclusion"]}}
    try:
        client = anthropic.Anthropic()
        import json as _json
        msg = client.messages.create(
            model=config.model, max_tokens=3500,
            system=(
                "You are a pharmacometrician CRITICALLY evaluating a PBPK model - "
                "not advertising it. Be precise and mechanistic.\n"
                "model_choice: justify the structure (distribution/permeability "
                "method, processes) from the biology.\n"
                "parameter_rationale: for EACH estimated parameter, COMPARE the "
                "fitted value against the literature_physicochemical anchor when "
                "one exists (e.g. fitted lipophilicity vs the reported logD/logP; "
                "fitted fraction unbound vs the reported range). If the fitted "
                "value AGREES with literature, say so - that is real support. If it "
                "DEPARTS from literature (e.g. lipophilicity fitted well below the "
                "measured logD), you MUST flag it explicitly and explain it "
                "(effective vs measured lipophilicity, or the optimizer "
                "compensating across correlated parameters) - do NOT rationalize a "
                "departure as fine just because it is inside a wide plausible "
                "range. Parameters identified only from plasma data and lacking a "
                "literature anchor (tissue/intestinal permeabilities) are typically "
                "WEAKLY IDENTIFIABLE - state that their absolute values are "
                "uncertain and should not be over-interpreted. A good fit with a "
                "parameter far from its expected value is a caution, not a success.\n"
                "Describe ONLY what happened: honor run_facts - if "
                "optimization_completed or parameters_were_fitted is false, do NOT "
                "claim parameters were fitted; call them un-fitted/preliminary. "
                "Match the tone to fit_verdict - never call a 'poor' fit adequate. "
                "Each parameter is labelled role=estimated or role=fixed - describe "
                "it as such (never call an estimated parameter 'measured'). Keep "
                "each section focused: 1-2 tight paragraphs; ALWAYS fill the "
                "conclusion. Call write_sections."),
            tools=[tool],
            messages=[{"role": "user", "content": _json.dumps(payload, default=str)}])
        for b in msg.content:
            if getattr(b, "type", None) == "tool_use":
                out = dict(b.input)
                # never ship an empty section (e.g. a long rationale starving the
                # conclusion): fill any blank from the deterministic write-up.
                det = None
                for key in ("model_choice", "parameter_rationale", "conclusion"):
                    if not str(out.get(key) or "").strip():
                        det = det or deterministic_narrative(d)
                        out[key] = det[key]
                return out
    except Exception:
        pass
    return deterministic_narrative(d)


# --------------------------------------------------------------------------- #
# inline SVG concentration-time plot (log y), pure stdlib
# --------------------------------------------------------------------------- #

def _svg_profile(study: str, obs, pred, ref=None, w=360, h=240) -> str:
    pts = [(t, c) for t, c in obs if c and c > 0] + \
          [(t, c) for t, c in pred if c and c > 0]
    if ref:
        pts += [(t, c) for t, c in ref if c and c > 0]
    if not pts:
        return f'<svg width="{w}" height="{h}"></svg>'
    xs = [t for t, _ in pts]; ys = [c for _, c in pts]
    x0, x1 = min(xs), max(xs) or 1
    ly = [math.log10(c) for c in ys]
    y0, y1 = min(ly), max(ly)
    if y1 == y0:
        y1 += 1
    pad = 34

    def px(t):
        return pad + (t - x0) / (x1 - x0 or 1) * (w - 2 * pad)

    def py(c):
        return h - pad - (math.log10(c) - y0) / (y1 - y0) * (h - 2 * pad)

    def path(series, color, dash=""):
        s = [(t, c) for t, c in series if c and c > 0]
        if len(s) < 2:
            return ""
        dd = " ".join(("M" if i == 0 else "L") + f"{px(t):.1f},{py(c):.1f}"
                      for i, (t, c) in enumerate(sorted(s)))
        da = f' stroke-dasharray="{dash}"' if dash else ""
        return f'<path d="{dd}" fill="none" stroke="{color}" stroke-width="1.8"{da}/>'

    dots = "".join(f'<circle cx="{px(t):.1f}" cy="{py(c):.1f}" r="2.6" '
                   f'fill="#111"/>' for t, c in obs if c and c > 0)
    axes = (f'<line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#999"/>'
            f'<line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#999"/>')
    labels = (f'<text x="{w/2}" y="{h-6}" font-size="10" text-anchor="middle">'
              f'time (h)</text>'
              f'<text x="10" y="{h/2}" font-size="10" text-anchor="middle" '
              f'transform="rotate(-90 10 {h/2})">conc (mg/L, log)</text>')
    return (f'<svg width="{w}" height="{h}" style="background:#fff">'
            f'<text x="{w/2}" y="14" font-size="11" text-anchor="middle" '
            f'font-weight="600">{html.escape(study[:44])}</text>'
            f'{axes}{path(pred,"#c0392b")}{path(ref,"#2980b9","4 3")}{dots}{labels}'
            f'</svg>')


# --------------------------------------------------------------------------- #
# HTML renderer
# --------------------------------------------------------------------------- #

def write_html(d: ReportData, path: str) -> None:
    def esc(x):
        return html.escape(str(x))

    def table(headers, rows):
        th = "".join(f"<th>{esc(h)}</th>" for h in headers)
        trs = "".join("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in r) + "</tr>"
                      for r in rows)
        return f"<table><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>"

    param_rows = [[p["name"], f"{p.get('value'):.4g}" if isinstance(p.get("value"), (int, float)) else p.get("value"),
                   p.get("unit", ""), p.get("role", ""),
                   p.get("plausible_range", ""),
                   f"{p.get('reference'):.4g}" if isinstance(p.get("reference"), (int, float)) else "-"]
                  for p in d.parameters]
    plots = "".join(f'<div class="plot">{_svg_profile(pf["study"], pf["obs"], pf["pred"], pf.get("ref"))}</div>'
                    for pf in d.profiles)
    traj = ""
    for i, s in enumerate(d.trajectory, 1):
        traj += (f'<div class="step"><div class="sh">Step {i} · {esc(s.get("action",""))}</div>'
                 + (f'<div class="reason">{esc(s.get("reason",""))}</div>' if s.get("reason") else "")
                 + (f'<pre class="res">{esc(s.get("result",""))}</pre>' if s.get("result") else "")
                 + "</div>")
    ode = d.odes or {}
    if ode:
        eqs_html = "".join(
            f'<div class="eq"><div class="eqn">{esc(name)}</div>'
            f'<pre class="eqf">{esc(formula)}</pre>'
            f'<div class="eqnote">{esc(note)}</div></div>'
            for name, formula, note in ode.get("equations", []))
        pmap = "".join(f"<li><b>{esc(n)}</b> — {esc(w)}</li>"
                       for n, w in ode.get("param_map", []))
        ode_html = (
            f'<p class="ode-comp">{esc(ode.get("compartments",""))}</p>'
            f'{eqs_html}'
            f'<p><b>Tissue partitioning.</b> {esc(ode.get("kp",""))}</p>'
            f'<p><b>Where the fitted / fixed parameters enter:</b></p><ul>{pmap}</ul>'
            f'<p class="eqnote"><i>{esc(ode.get("caveat",""))}</i></p>')
    else:
        ode_html = "<p>(not available)</p>"
    fit = d.fit
    ref = d.reference
    by_route = "".join(f"<li>{esc(r)}: GMFE {esc(m.get('gmfe'))}, bias {esc(m.get('bias'))}, "
                       f"within-2fold {esc(m.get('within_2fold_pct'))}%</li>"
                       for r, m in (fit.get("by_route") or {}).items())
    doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>{esc(d.title)}</title>
<style>
 body{{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;color:#1a1a1a;max-width:860px;margin:24px auto;padding:0 18px}}
 h1{{font-size:24px;border-bottom:2px solid #333;padding-bottom:6px}}
 h2{{font-size:18px;margin-top:28px;border-bottom:1px solid #ddd;padding-bottom:4px}}
 table{{border-collapse:collapse;width:100%;margin:8px 0;font-size:13px}}
 th,td{{border:1px solid #ccc;padding:4px 8px;text-align:left}} th{{background:#f4f4f4}}
 .plots{{display:flex;flex-wrap:wrap;gap:8px}} .plot{{border:1px solid #eee}}
 .legend{{font-size:12px;color:#555}} .legend b.o{{color:#111}} .legend b.p{{color:#c0392b}} .legend b.r{{color:#2980b9}}
 .step{{border-left:3px solid #c0392b;padding:4px 10px;margin:8px 0;background:#fafafa}}
 .sh{{font-weight:600}} .reason{{color:#333;white-space:pre-wrap;margin:4px 0}}
 pre.res{{background:#f0f0f0;padding:6px;overflow-x:auto;font-size:12px;white-space:pre-wrap}}
 ul{{margin:6px 0}}
 .banner{{padding:10px 14px;margin:12px 0;border-radius:5px;font-weight:600}}
 .banner.warn{{background:#fdecea;border:1px solid #e74c3c;color:#922}}
 .banner.ok{{background:#eafaf1;border:1px solid #27ae60;color:#1a6b3c}}
 .ode-comp{{color:#333}}
 .eq{{margin:10px 0;border-left:3px solid #2980b9;padding:2px 12px}}
 .eqn{{font-weight:600;font-size:13px}}
 pre.eqf{{background:#f4f7fb;border:1px solid #e1e8f0;padding:8px 10px;margin:4px 0;
   overflow-x:auto;font-size:13px;font-family:"Cambria Math",Consolas,monospace}}
 .eqnote{{color:#555;font-size:12.5px}}
 @media print{{.step{{break-inside:avoid}} .plot{{break-inside:avoid}}}}
</style></head><body>
<h1>{esc(d.title)}</h1>
{f'<div class="banner {"warn" if d.status.startswith(chr(0x26A0)) else "ok"}">{esc(d.status)}</div>' if d.status else ''}

<h2>1. Objective</h2><p>{esc(d.objective)}</p>
<p><b>Background.</b> {esc(d.background)}</p>
<b>Known biology:</b><ul>{''.join(f'<li>{esc(x)}</li>' for x in d.known_biology)}</ul>

<h2>2. Data</h2>
<p>{esc(d.data_overview.get('n_datasets'))} observed plasma datasets;
routes {esc(', '.join(str(r) for r in (d.data_overview.get('routes') or [])))}.</p>
{table(['Study','C_max (mg/L)','T_max (h)','AUC (mg·h/L)','t½ (h)'],
       [[r.get('study'), r.get('c_max_mg_L'), r.get('t_max_h'), r.get('auc_mg_h_L'), r.get('t_half_h')] for r in d.nca_rows])}

<h2>3. Model &amp; model choice</h2>
<p><b>Distribution / permeability:</b> {esc(', '.join(d.structure.get('calculation_methods') or []))}<br>
<b>Processes:</b> {esc(', '.join(d.structure.get('processes') or []))}</p>
<p>{esc(d.narrative.get('model_choice',''))}</p>

<h2>4. Model equations (ODE system)</h2>
{ode_html}

<h2>5. Parameters</h2>
{table(['Parameter','Value','Unit','Role','Plausible range','Reference (truth)'], param_rows)}

<h2>6. Pharmacological rationale</h2>
<p style="white-space:pre-wrap">{esc(d.narrative.get('parameter_rationale',''))}</p>

<h2>7. Concentration-time analysis</h2>
<p>Overall <b>GMFE {esc(fit.get('gmfe'))}</b>, within-2-fold {esc(fit.get('within_2fold_pct'))}%.</p>
<ul>{by_route}</ul>
<p class="legend"><b class="o">●</b> observed &nbsp; <b class="p">— agent model</b> &nbsp; <b class="r">--- reference model</b></p>
<div class="plots">{plots}</div>

<h2>8. Comparison with the ground-truth model</h2>
<p>Reference (published) overall GMFE: <b>{esc(ref.get('gmfe'))}</b> vs this model {esc(fit.get('gmfe'))}.</p>
{table(['Parameter','This model','Reference','Fold vs reference'],
       [[p['name'], f"{p.get('value'):.4g}" if isinstance(p.get('value'),(int,float)) else p.get('value'),
         f"{p.get('reference'):.4g}" if isinstance(p.get('reference'),(int,float)) else '-',
         f"{(p['value']/p['reference']):.2g}" if isinstance(p.get('reference'),(int,float)) and p.get('reference') and isinstance(p.get('value'),(int,float)) else '-']
        for p in d.parameters])}

<h2>9. Conclusion</h2><p>{esc(d.narrative.get('conclusion',''))}</p>

<h2>10. Full modeling trajectory (LLM)</h2>{traj}
</body></html>"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


# --------------------------------------------------------------------------- #
# PDF renderer (matplotlib; on the user's machine)
# --------------------------------------------------------------------------- #

def write_pdf(d: ReportData, path: str) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except Exception:
        return False
    import textwrap

    def text_page(pdf, title, blocks):
        # blocks: list of strings; paginate
        lines = []
        for b in blocks:
            for ln in b.split("\n"):
                lines += textwrap.wrap(ln, 100) or [""]
        per = 46
        for start in range(0, max(1, len(lines)), per):
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.08, 0.95, title, fontsize=14, fontweight="bold")
            fig.text(0.08, 0.92, "\n".join(lines[start:start + per]), fontsize=9,
                     va="top", family="monospace")
            pdf.savefig(fig); plt.close(fig)

    with PdfPages(path) as pdf:
        # title + summary
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.08, 0.92, d.title, fontsize=18, fontweight="bold")
        y = 0.87
        if d.status:
            warn = d.status.startswith("⚠")
            fig.text(0.08, y, "\n".join(textwrap.wrap(d.status, 90)), fontsize=10,
                     va="top", fontweight="bold",
                     color=("#b03030" if warn else "#1a6b3c"),
                     bbox=dict(boxstyle="round", facecolor=("#fdecea" if warn
                               else "#eafaf1"), edgecolor=("#e74c3c" if warn
                               else "#27ae60")))
            y -= 0.07
        summ = (f"Objective:\n{d.objective}\n\n"
                f"Overall GMFE {d.fit.get('gmfe')} "
                f"(reference {d.reference.get('gmfe')}); "
                f"within-2-fold {d.fit.get('within_2fold_pct')}%.\n\n"
                f"Structure: {', '.join(d.structure.get('calculation_methods') or [])}; "
                f"processes {', '.join(d.structure.get('processes') or [])}.")
        fig.text(0.08, y, "\n".join(textwrap.wrap(summ, 95)), fontsize=10, va="top")
        pdf.savefig(fig); plt.close(fig)

        # narrative + parameters
        pr = ["MODEL CHOICE", d.narrative.get("model_choice", ""), "",
              "PARAMETERS"]
        for p in d.parameters:
            pr.append(f"  {p['name']} = {p.get('value')} {p.get('unit','')}  "
                      f"[{p.get('role','')}] ref={p.get('reference','-')}")
        pr += ["", "PHARMACOLOGICAL RATIONALE", d.narrative.get("parameter_rationale", "")]
        text_page(pdf, "Model, parameters & rationale", pr)

        # model equations (ODE system)
        ode = d.odes or {}
        if ode:
            ob = ["COMPARTMENTS / STATE VARIABLES", ode.get("compartments", ""), "",
                  "GOVERNING EQUATIONS"]
            for name, formula, note in ode.get("equations", []):
                ob += [f"  {name}:", f"      {formula}", f"      ({note})", ""]
            ob += ["TISSUE PARTITIONING", ode.get("kp", ""), "",
                   "WHERE THE FITTED / FIXED PARAMETERS ENTER"]
            ob += [f"  - {n}: {w}" for n, w in ode.get("param_map", [])]
            ob += ["", ode.get("caveat", "")]
            text_page(pdf, "Model equations (ODE system)", ob)

        # concentration-time plots (grid)
        profs = d.profiles
        for start in range(0, len(profs), 6):
            fig, axes = plt.subplots(2, 3, figsize=(8.27, 11.69))
            for ax, pf in zip(axes.flat, profs[start:start + 6]):
                o = [(t, c) for t, c in pf["obs"] if c and c > 0]
                pr_ = [(t, c) for t, c in pf["pred"] if c and c > 0]
                if o:
                    ax.scatter([t for t, _ in o], [c for _, c in o], s=10, c="k")
                if pr_:
                    pr_.sort()
                    ax.plot([t for t, _ in pr_], [c for _, c in pr_], "r-")
                if pf.get("ref"):
                    rf = sorted((t, c) for t, c in pf["ref"] if c and c > 0)
                    if rf:
                        ax.plot([t for t, _ in rf], [c for _, c in rf], "b--")
                ax.set_yscale("log"); ax.set_title(pf["study"][:30], fontsize=7)
                ax.tick_params(labelsize=6)
            for ax in axes.flat[len(profs[start:start + 6]):]:
                ax.axis("off")
            fig.suptitle("Observed (●) vs agent (—) vs reference (--)", fontsize=10)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # comparison + conclusion
        text_page(pdf, "Ground-truth comparison & conclusion",
                  [f"Reference GMFE {d.reference.get('gmfe')} vs this model "
                   f"{d.fit.get('gmfe')}", "",
                   "CONCLUSION", d.narrative.get("conclusion", "")])

        # trajectory
        tb = []
        for i, s in enumerate(d.trajectory, 1):
            tb.append(f"[Step {i}] {s.get('action','')}")
            if s.get("reason"):
                tb.append("  reason: " + s["reason"])
            if s.get("result"):
                tb.append("  -> " + s["result"])
            tb.append("")
        text_page(pdf, "Full modeling trajectory (LLM)", tb)
    return True

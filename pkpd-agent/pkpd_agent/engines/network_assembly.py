r"""Assemble the FULL coupled immune network - every cell (template A) and every cytokine
(template B) dynamic at once - into ONE SBML, calibrated to the joint steady state. This is the
scale-integration step the isolated hub/cell/pair builders never exercised: does a ~20-species
network of the model's OWN edges, each free rate pinned by its own target, hold together as one
self-consistent ODE system?

Coupling (both directions live, nothing clamped):
    d(Cyt)/dt  = SUM_cells [ ksec_Cyt * (Cell/target_Cell) * PROD(secretion modifiers) ]
                 - kcl_Cyt * Cyt
    d(Cell)/dt = kIn * PROD(influx eff) + kprolif * Cell * PROD(prolif eff)
                 - kd * Cell * PROD(death eff)

Secretion scales with the secreting cell's abundance (Cell/target_Cell), so a change in any cell
feeds back into the cytokines it makes, which in turn drive the cells - the loop the single-node
probes could only bound.

CALIBRATION is well-posed because at the joint target state every species sits at its target, so
each balance equation contains exactly ONE free rate and is solved independently:
    ksec_Cyt  (one per cytokine)  <- its clearance balance
    kprolif_Cell (one per cell)   <- its life-cycle balance
The PER-CELL secretion split (which of a cytokine's secreting cells makes how much) is NOT pinned
by a single cytokine target - that is the same identifiability bottleneck seen for the K
half-effects, so secretion is lumped into one cell-count-weighted scale per cytokine and flagged.

Species without a steady-state target (GMCSF, AutoAb) cannot be calibrated; their edges are dropped
and reported. Drug/PK/score/latch species (the rest of the paper's 59) are out of scope here - they
need the MATLAB SimBiology engine, not this dependency-free integrator.
"""

from __future__ import annotations

import re

from pkpd_agent.engines import cell_lifecycle as CL

_SEC = re.compile(r"(?i)^([A-Za-z0-9]+)Sec([A-Za-z0-9]+?)_Maxby([A-Za-z0-9]+)$")


def _kcl_of(prov, cyt):
    """The clearance param+value for a cytokine, tolerating name variants (kcl_IL6 vs kcl_IL-6)."""
    key = cyt.replace("-", "").lower()
    for n, p in prov.items():
        m = re.match(r"(?i)^kcl_(.+)$", n)
        if m and m.group(1).replace("-", "").lower() == key and p.get("value_from_reference") \
                is not None:
            return n, float(p["value_from_reference"])
    return None, None


def discover_secretion(prov, cell_token_to_canonical):
    """{cytokine: {cell: {modifiers: {mod:(max,lit)}}}} from <Cyt>Sec<Cell>_Maxby<Mod> params.
    Cell tokens are mapped to canonical cell species via the provided map; unknown cells are kept
    under their raw token so the caller can drop them."""
    out = {}
    for n, p in prov.items():
        m = _SEC.match(n)
        if not m:
            continue
        cyt, cell_tok, mod = m.group(1), m.group(2), m.group(3)
        cell = cell_token_to_canonical.get(cell_tok, cell_tok)
        out.setdefault(cyt, {}).setdefault(cell, {})[mod] = (
            p.get("value_from_reference"), bool(p.get("from_literature")))
    return out


def cell_token_map(aliases):
    """Invert canonical->[prefixes] into secreting-cell-token -> canonical."""
    inv = {}
    for canon, prefs in aliases.items():
        for pfx in prefs:
            inv[pfx] = canon
    return inv


def _sec_effect(mods, levels, prior=1.5):
    """PROD_mod (1 + (Max-1)*0.5) at baseline, skipping modifiers with no level."""
    eff = 1.0
    for mod, (mx, _lit) in mods.items():
        if mod not in levels:
            continue
        m = prior if mx is None else float(mx)
        eff *= 1.0 + (m - 1.0) * 0.5
    return eff


def assemble_network(prov, levels, cells, aliases, prior=1.5, sec_override=None):
    """Build the full coupled spec + calibrate every free rate. Returns (spec, meta) where meta
    records the dropped (uncalibratable) edges, the free rates, and per-cytokine secreting cells.

    By default the STRUCTURE (which cell secretes which cytokine, which cytokine modulates each
    secretion and each cell flux) is the model's own, read from the params. Pass ``sec_override``
    (same shape as discover_secretion) and a ``cells`` whose prolif/influx/apop sets have been
    replaced to assemble from an AGENT's proposed structure instead - the constant VALUES are still
    looked up from the model where the edge exists (prior otherwise) and every free rate is
    re-calibrated to the targets, so any structure yields a self-consistent network and the quality
    difference shows up in recall/precision and the held-out response, not in the calibration."""
    sec = sec_override if sec_override is not None else \
        discover_secretion(prov, cell_token_map(aliases))
    cyt_levels = {c: levels[c] for c in levels}                    # cytokines-with-level are dynamic
    dyn_cells = set(cells)
    cell_target = {c: cells[c]["target"] for c in cells}           # cells carry their own target
    # cytokines we can make dynamic: have a level AND a clearance
    dyn_cyts, kcl = {}, {}
    for c in list(sec) + list(cyt_levels):
        if c in dyn_cyts or c not in levels:
            continue
        kp, kv = _kcl_of(prov, c)
        if kp is not None:
            dyn_cyts[c] = levels[c]; kcl[c] = (kp, kv)

    species, params, reactions = [], {}, []
    dropped = {"no_level_or_clearance": [], "secreting_cell_not_dynamic": []}

    # ---- cytokine secretion + clearance ----
    sec_cells = {}
    for cyt in sorted(dyn_cyts):
        contribs = []          # (cell, modifiers)
        for cell, mods in sec.get(cyt, {}).items():
            if cell not in dyn_cells:
                dropped["secreting_cell_not_dynamic"].append(f"{cyt}<-{cell}")
                continue
            contribs.append((cell, mods))
        sec_cells[cyt] = sorted(c for c, _ in contribs)
        ksec = f"ksec_{cyt}"
        # calibrate: ksec * SUM_cells eff_cell = kcl * level   (at target, Cell/target=1)
        denom = sum(_sec_effect(mods, levels, prior) for _, mods in contribs)
        kcl_param, kcl_val = kcl[cyt]
        params[kcl_param] = kcl_val
        if denom > 0:
            params[ksec] = kcl_val * levels[cyt] / denom
            for cell, mods in contribs:
                factors = []
                for mod, (mx, _l) in mods.items():
                    if mod not in levels:
                        continue
                    m = prior if mx is None else float(mx)
                    mp, kp = f"Msec_{cyt}_{cell}_{mod}", f"Ksec_{cyt}_{mod}"
                    params[mp] = m; params.setdefault(kp, levels[mod])
                    factors.append(f"(1 + ({mp} - 1) * {mod} / ({kp} + {mod}))")
                prodf = " * ".join(factors) or "1"
                reactions.append({"id": f"{cyt}_sec_{cell}", "reactants": [], "products": [cyt],
                                  "rate": f"{ksec} * ({cell} / {cell_target[cell]:.10g}) * "
                                          f"({prodf})"})
        else:                                                       # no dynamic secreting cell
            params[ksec] = kcl_val * levels[cyt]                    # constant zeroth-order fallback
            reactions.append({"id": f"{cyt}_sec_const", "reactants": [], "products": [cyt],
                              "rate": ksec})
        reactions.append({"id": f"{cyt}_clr", "reactants": [cyt], "products": [],
                          "rate": f"{kcl_param} * {cyt}"})
        species.append({"name": cyt, "initial": levels[cyt]})

    # ---- cell life-cycle (regulators are the dynamic cytokines) ----
    free_kprolif = {}
    for cell in sorted(dyn_cells):
        info = cells[cell]
        kpp = f"kprolif_{cell}"
        kp, marg = CL.fit_base_prolif(info, levels, prior=prior)
        rxns, vals = CL.cell_reactions(cell, info, levels, kpp, prior=prior)
        vals[kpp] = kp
        free_kprolif[cell] = (kp, marg)
        params.update(vals)
        reactions += rxns
        species.append({"name": cell, "initial": info["target"]})

    spec = {"name": "ra_immune_network", "species": species,
            "parameters": [{"name": k, "value": v} for k, v in params.items()],
            "reactions": reactions, "rules": []}
    meta = {"dynamic_cytokines": sorted(dyn_cyts), "dynamic_cells": sorted(dyn_cells),
            "secreting_cells": sec_cells, "dropped": dropped,
            "free_ksec": [f"ksec_{c}" for c in sorted(dyn_cyts)],
            "free_kprolif": {c: free_kprolif[c] for c in sorted(free_kprolif)}}
    return spec, meta


def apply_structure(model_sec, cells, sec_struct, cell_struct, prior_flag=(None, False)):
    """Substitute an AGENT's chosen structure into the model-shaped dicts, keeping the constant
    VALUES from the model where the edge exists and marking the rest as uncited (prior at assembly).

      sec_struct:  {cytokine: {cell: [modulator cytokines]}}   - agent's secretion edges
      cell_struct: {cell: {'prolif'|'influx'|'apop': [cytokines]}} - agent's cell-flux edges

    Returns (sec2, cells2) in the same shapes discover_secretion / discover_cells produce, so
    assemble_network(..., sec_override=sec2) with cells2 builds the agent's network."""
    sec2 = {}
    for cyt, per_cell in sec_struct.items():
        for cell, mods in per_cell.items():
            for mod in mods:
                val = model_sec.get(cyt, {}).get(cell, {}).get(mod, prior_flag)
                sec2.setdefault(cyt, {}).setdefault(cell, {})[mod] = val
    cells2 = {}
    for cell, info in cells.items():
        new = dict(info)
        chosen = cell_struct.get(cell, {})
        for flux in ("prolif", "influx", "apop"):
            new[flux] = {cyt: info[flux].get(cyt, prior_flag) for cyt in chosen.get(flux, [])}
        cells2[cell] = new
    return sec2, cells2


def _score(chosen, truth):
    chosen, truth = set(chosen), set(truth)
    hit = chosen & truth
    return {"recall": round(len(hit) / len(truth), 3) if truth else 0.0,
            "precision": round(len(hit) / len(chosen), 3) if chosen else 0.0,
            "chosen": sorted(chosen), "truth": sorted(truth),
            "missed": sorted(truth - chosen), "extra": sorted(chosen - truth)}


_SECRETE_SYS = (
    "You are deciding a QSP model's STRUCTURE from immunology, before looking up any values. Given "
    "a cytokine and a list of CELL TYPES, name which of those cell types are meaningful cellular "
    "SOURCES of that cytokine in this disease. Reason from biology; do NOT assume a parameter "
    "table. Use ONLY the exact cell-type names given. JSON only.")


def _propose_secreting_cells(cyt, cell_candidates, call):
    """Which cell types secrete this cytokine - a source (not regulator) question, so it gets its
    own cell-appropriate prompt rather than the cytokine-regulator one."""
    from pkpd_agent.engines.llm_structure import _parse_json
    user = (f"Cytokine: {cyt}. Candidate cell types: " + ", ".join(cell_candidates) +
            f'\n\nWhich of these cell types secrete {cyt}? Return JSON '
            '{"cells": ["name", ...]}. Only names from the list.')
    d = _parse_json(call(_SECRETE_SYS, user))
    cset = set(cell_candidates)
    return [c for c in (d.get("cells") or []) if c in cset]


def propose_structure(prov, levels, cells, aliases, call, log=lambda *a: None):
    """The FULL agent structural pass: for every given node the agent proposes its edges, scored
    against the model's own. Returns (sec_struct, cell_struct, scores). Nodes are GIVEN (not the
    agent's job); the agent decides topology only. Constants are looked up afterwards in assembly."""
    from pkpd_agent.engines import model_assembly as MA
    model_sec = discover_secretion(prov, cell_token_map(aliases))
    dyn_cells = sorted(cells)
    dyn_cyts = sorted(c for c in levels if _kcl_of(prov, c)[0] and c not in dyn_cells)
    scores = {"secreting_cells": {}, "secretion_mods": {}, "cell_flux": {}}

    # --- per cytokine: which cells secrete it, and which cytokines modulate that secretion ---
    sec_struct = {}
    for cyt in dyn_cyts:
        cell_cands = [c for c in dyn_cells]
        secreting = _propose_secreting_cells(cyt, cell_cands, call)
        truth_cells = [c for c in model_sec.get(cyt, {}) if c in dyn_cells]
        scores["secreting_cells"][cyt] = _score(secreting, truth_cells)
        mod_cands = [c for c in dyn_cyts if c != cyt]
        mregs = MA.propose_regulators(cyt, mod_cands, "secretion (which cytokines up/down-regulate "
                                      "how much of it is secreted)", call)
        mods = [r["cytokine"] for r in mregs if r["cytokine"] in mod_cands]
        truth_mods = sorted({m for c in model_sec.get(cyt, {}).values() for m in c if m in levels})
        scores["secretion_mods"][cyt] = _score(mods, truth_mods)
        if secreting:
            sec_struct[cyt] = {cell: list(mods) for cell in secreting}
        log(f"  cyt {cyt:7} cells r={scores['secreting_cells'][cyt]['recall']:.2f} "
            f"p={scores['secreting_cells'][cyt]['precision']:.2f}  mods "
            f"r={scores['secretion_mods'][cyt]['recall']:.2f} "
            f"p={scores['secretion_mods'][cyt]['precision']:.2f}")

    # --- per cell: proliferation / influx / apoptosis regulators ---
    cell_struct = {}
    proc = {"prolif": "proliferation", "influx": "recruitment (influx from blood)",
            "apop": "apoptosis (death)"}
    cyt_cands = list(dyn_cyts)
    for cell in dyn_cells:
        cell_struct[cell] = {}
        for flux, pname in proc.items():
            regs = MA.propose_regulators(cell, cyt_cands, pname, call)
            picks = [r["cytokine"] for r in regs if r["cytokine"] in cyt_cands]
            truth = [c for c in cells[cell][flux] if c in levels]
            scores["cell_flux"][f"{cell}.{flux}"] = _score(picks, truth)
            cell_struct[cell][flux] = picks
        r = [scores["cell_flux"][f"{cell}.{f}"]["recall"] for f in proc]
        log(f"  cell {cell:11} flux recall {['%.2f'%x for x in r]}")
    return sec_struct, cell_struct, scores


def integrate_network(spec, clamp=None, t_end=60.0, dt=5e-3, diverge_fold=1e9):
    """RK4 on a spec directly, with each reaction rate COMPILED once (not re-parsed per step) - fast
    enough for the full stiff network. Clamped/boundary species are held fixed; the rest advance by
    sum(production) - sum(consumption). Returns the final {species: value} state.

    dt must satisfy the explicit-stability bound of the fastest clearance (here dt<~0.008 for the
    332/time TGFb clearance); the default 5e-3 is safe for this model. If any species runs away past
    ``diverge_fold`` x its starting value or goes non-finite, integration stops early and the result
    carries ``state['__diverged__'] = True`` - an UNSTABLE network (e.g. one an over-inclusive
    topology made positive-feedback-unstable), not a numeric artifact to read as a real level."""
    import math
    clamp = dict(clamp or {})
    params = {p["name"]: float(p["value"]) for p in spec.get("parameters", [])}
    held = {s["name"] for s in spec["species"] if s.get("boundary")} | set(clamp)
    dyn = [s["name"] for s in spec["species"] if s["name"] not in held]
    state = {s["name"]: float(s.get("initial", 0.0)) for s in spec["species"]}
    state.update(clamp)
    rxns = []
    for r in spec["reactions"]:
        code = compile((r.get("rate") or "0").replace("^", "**"), "<rate>", "eval")
        rxns.append(([p for p in r.get("products", []) if p in set(dyn)],
                     [x for x in r.get("reactants", []) if x in set(dyn)], code))
    env = {"min": min, "max": max, **params}

    def deriv(st):
        env.update(st)
        d = {k: 0.0 for k in dyn}
        for prod, reac, code in rxns:
            v = eval(code, {"__builtins__": {}}, env)      # our own compiled expression
            for p in prod:
                d[p] += v
            for x in reac:
                d[x] -= v
        return d

    start = {k: state[k] for k in dyn}
    limit = {k: diverge_fold * (abs(start[k]) or 1.0) for k in dyn}
    diverged = False
    for _ in range(int(t_end / dt)):
        try:
            k1 = deriv(state)
            s2 = {**state, **{k: state[k] + 0.5 * dt * k1[k] for k in dyn}}
            k2 = deriv(s2)
            s3 = {**state, **{k: state[k] + 0.5 * dt * k2[k] for k in dyn}}
            k3 = deriv(s3)
            s4 = {**state, **{k: state[k] + dt * k3[k] for k in dyn}}
            k4 = deriv(s4)
            for k in dyn:
                state[k] += dt / 6.0 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k])
        except (OverflowError, ValueError):
            diverged = True; break
        env.update(clamp)
        if any(not math.isfinite(state[k]) or abs(state[k]) > limit[k] for k in dyn):
            diverged = True; break
    out = {k: state[k] for k in state}
    if diverged:
        out["__diverged__"] = True
    return out

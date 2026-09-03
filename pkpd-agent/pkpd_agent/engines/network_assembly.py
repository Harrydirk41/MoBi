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


def assemble_network(prov, levels, cells, aliases, prior=1.5):
    """Build the full coupled spec + calibrate every free rate. Returns (spec, meta) where meta
    records the dropped (uncalibratable) edges, the free rates, and per-cytokine secreting cells."""
    sec = discover_secretion(prov, cell_token_map(aliases))
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


def integrate_network(spec, clamp=None, t_end=60.0, dt=5e-3):
    """RK4 on a spec directly, with each reaction rate COMPILED once (not re-parsed per step) - fast
    enough for the full stiff network. Clamped/boundary species are held fixed; the rest advance by
    sum(production) - sum(consumption). Returns the final {species: value} state.

    dt must satisfy the explicit-stability bound of the fastest clearance (here dt<~0.008 for the
    332/time TGFb clearance); the default 5e-3 is safe for this model."""
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

    for _ in range(int(t_end / dt)):
        k1 = deriv(state)
        s2 = {**state, **{k: state[k] + 0.5 * dt * k1[k] for k in dyn}}
        k2 = deriv(s2)
        s3 = {**state, **{k: state[k] + 0.5 * dt * k2[k] for k in dyn}}
        k3 = deriv(s3)
        s4 = {**state, **{k: state[k] + dt * k3[k] for k in dyn}}
        k4 = deriv(s4)
        for k in dyn:
            state[k] += dt / 6.0 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k])
        env.update(clamp)
    return {k: state[k] for k in state}

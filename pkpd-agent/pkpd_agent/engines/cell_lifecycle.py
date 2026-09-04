r"""Template A - the cell life-cycle ODE, general and node-agnostic.

The Vantage RA model has two ODE templates. Template B (cytokine secretion) is handled by the
build_general / model_assembly path: d(Cyt)/dt = production - clearance. Template A is the cell
life-cycle, and this module builds it from the model's OWN naming conventions, with no per-cell
hardcoding:

    d(Cell)/dt =  kIn_base * PROD(migration effects)          # influx / recruitment  (zeroth order)
                + kprolif   * Cell * PROD(prolif effects)      # proliferation        (first order)
                - kd_base   * Cell * PROD(apoptosis effects)   # death / apoptosis    (first order)

    each effect PROD = PROD_r ( 1 + (Max_r - 1) * Cyt_r/(K_r + Cyt_r) )   (K_r = baseline level)

So a cell has three regulator SETS (proliferation, influx, apoptosis) discovered from the params
    <Cell>Prolif_Maxby<Cyt>,  <Cell>Influx_Maxby<Cyt>,  <Cell>Apop_Maxby<Cyt>
and two literature baseline rates
    kd_<Cell>_Baseline  (death, always present),  kIn_<Cell>_Baseline  (influx, some cells)

WHICH rate is fitted vs looked up (the honest provenance split, same story as template B's kg):
    kd_base  - literature (cell turnover / half-life)                 -> LOOKED UP
    kIn_base - literature where cited, else fitted                    -> MIXED
    kprolif  - never in the table (no citation) -> the FREE rate      -> FITTED to the target
The single steady-state target pins exactly one free rate (kprolif), as one balance equation pins
one unknown. Cells with NO influx (FLS, PlasmaCells) are first-order birth-death: their steady
state is marginal (birth=death fixes the RATE, not the level), so the target is imposed as the
initial condition and kprolif is set to balance death there - flagged in the fit provenance.

Cell prefixes are irregular in this model (abbreviations, a 'Macrophase' typo); the canonical->
prefixes map lives in the project's cell_prefix_aliases.json so this engine stays general.
"""

from __future__ import annotations

import json
import os
import re

_FLUX_PATTERNS = {"prolif": "prolif", "influx": "influx", "apop": "apop"}


def load_cell_aliases(project_dir):
    """{canonical_cell: [name-prefixes]} from the project's config, or {} if none (identity)."""
    path = os.path.join(project_dir, "data", "cell_prefix_aliases.json")
    if not os.path.isfile(path):
        return {}
    return json.load(open(path)).get("canonical_to_prefixes", {})


def _prefix_regex(prefixes):
    return "|".join(re.escape(p) for p in sorted(prefixes, key=len, reverse=True))


def _maxby_set(prov, prefixes, flux):
    """{cyt: (max_value, from_literature)} for <prefix><flux>_Maxby<Cyt> params (case-insensitive).
    Nested modifiers <prefix><flux>_by<x>_Maxby<Cyt> are skipped (ambiguous secondary effects)."""
    pat = re.compile(rf"(?i)^(?:{_prefix_regex(prefixes)}){flux}_Maxby([A-Za-z0-9]+)$")
    out = {}
    for n, p in prov.items():
        m = pat.match(n)
        if m:
            out[m.group(1)] = (p.get("value_from_reference"), bool(p.get("from_literature")))
    return out


def _baseline(prov, prefixes, kind):
    """The kd_<Cell>_Baseline or kIn_<Cell>_Baseline param name + (value, from_literature)."""
    pat = re.compile(rf"(?i)^{kind}_(?:{_prefix_regex(prefixes)})_Baseline$")
    for n, p in prov.items():
        if pat.match(n):
            return n, p.get("value_from_reference"), bool(p.get("from_literature"))
    return None, None, None


def discover_cells(prov, targets, aliases):
    """Every cell that has a steady-state target and at least a death baseline is buildable.
    Returns {canonical: {prolif, influx, apop (each {cyt:(max,lit)}), kd_param, kd_val,
    kin_param, kin_val, kin_lit, target}} - all from the model's own conventions."""
    cell_targets = {t["model_species"]: float(t["target_model_unit"]) for t in targets
                    if t.get("kind") == "cell" and t.get("model_species")
                    and t.get("target_model_unit") is not None}
    out = {}
    for cell, target in cell_targets.items():
        prefixes = aliases.get(cell, [cell])
        kd_param, kd_val, _ = _baseline(prov, prefixes, "kd")
        if not kd_param:
            continue                                  # no death rate -> not a life-cycle node
        kin_param, kin_val, kin_lit = _baseline(prov, prefixes, "kIn")
        out[cell] = {
            "prolif": _maxby_set(prov, prefixes, "prolif"),
            "influx": _maxby_set(prov, prefixes, "influx"),
            "apop": _maxby_set(prov, prefixes, "apop"),
            "kd_param": kd_param, "kd_val": float(kd_val) if kd_val is not None else None,
            "kin_param": kin_param, "kin_val": float(kin_val) if kin_val is not None else None,
            "kin_lit": kin_lit, "target": target}
    return out


def _effect_at(regs, levels, chosen=None, prior=1.5):
    """PROD_r (1 + (Max_r-1)*hill) evaluated at baseline levels (hill=0.5 when K=level). ``regs`` is
    {cyt:(max,lit)}; restrict to ``chosen`` if given; an uncited (None) max uses the generic prior."""
    eff = 1.0
    for cyt, (mx, _lit) in regs.items():
        if chosen is not None and cyt not in chosen:
            continue
        if cyt not in levels:
            continue
        m = prior if mx is None else float(mx)
        eff *= 1.0 + (m - 1.0) * 0.5               # hill(level, K=level) = 0.5
    return eff


def fit_base_prolif(info, levels, chosen=None, prior=1.5):
    """Solve d(Cell)/dt = 0 at the target for the one free rate kprolif:
        0 = kIn*mig + kprolif*target*prolif - kd*target*apop
    -> kprolif = (kd*target*apop - kIn*mig) / (target*prolif).
    Returns (kprolif, marginal) where marginal=True means no influx (birth-death only: kprolif set
    so births balance deaths at the target IC; the level itself is not pinned by steady state)."""
    target, kd = info["target"], info["kd_val"]
    mig = _effect_at(info["influx"], levels, chosen, prior)
    pro = _effect_at(info["prolif"], levels, chosen, prior)
    apo = _effect_at(info["apop"], levels, chosen, prior)
    kin = info["kin_val"] if (info["kin_param"] and info["kin_val"] is not None) else 0.0
    marginal = not (info["kin_param"] and info["kin_val"] is not None and info["kin_val"] > 0)
    if marginal:
        kprolif = kd * apo / pro                    # births balance deaths at the target IC
    else:
        kprolif = (kd * target * apo - kin * mig) / (target * pro)
    return kprolif, marginal


def cell_reactions(cell, info, levels, kprolif_param, chosen=None, prior=1.5,
                   max_of=None, k_of=None):
    """Build the three life-cycle reactions for a cell as {id,reactants,products,rate} dicts, plus
    the parameter values they reference. ``chosen`` restricts each flux's regulators to an agent's
    picks (default: the model's full set). ``max_of``/``k_of`` name the per-edge Max/K params so
    the caller can keep them unique across a multi-cell network; defaults are M_<cell>_<flux>_<cyt>
    and K_<cell>_<cyt>. Returns (reactions, values)."""
    max_of = max_of or (lambda flux, cyt: f"M_{cell}_{flux}_{cyt}")
    k_of = k_of or (lambda cyt: f"K_{cell}_{cyt}")
    values = {}

    def factors(regs, flux):
        specs = []
        for cyt, (mx, _lit) in regs.items():
            if (chosen is not None and cyt not in chosen) or cyt not in levels:
                continue
            m = prior if mx is None else float(mx)
            mp, kp = max_of(flux, cyt), k_of(cyt)
            values[mp] = m
            values.setdefault(kp, levels[cyt])
            specs.append(f"(1 + ({mp} - 1) * {cyt} / ({kp} + {cyt}))")
        return " * ".join(specs) or "1"

    rxns = []
    # influx (zeroth order) - only when the model gives this cell an influx baseline
    if info["kin_param"] and info["kin_val"] is not None and info["kin_val"] > 0:
        values[info["kin_param"]] = info["kin_val"]
        rxns.append({"id": f"{cell}_influx", "reactants": [], "products": [cell],
                     "rate": f"{info['kin_param']} * ({factors(info['influx'], 'influx')})"})
    # proliferation (first order in Cell) - the fitted free rate
    rxns.append({"id": f"{cell}_prolif", "reactants": [], "products": [cell],
                 "rate": f"{kprolif_param} * {cell} * ({factors(info['prolif'], 'prolif')})"})
    # death / apoptosis (first order in Cell) - literature baseline
    values[info["kd_param"]] = info["kd_val"]
    rxns.append({"id": f"{cell}_death", "reactants": [cell], "products": [],
                 "rate": f"{info['kd_param']} * {cell} * ({factors(info['apop'], 'apop')})"})
    return rxns, values


def synthesize_influx(info, levels, frac=0.5, prior=1.5):
    """WHAT-IF (not the honest build): give a marginal cell - one the model gives no literature
    influx rate, so it was built as pure birth-death - a synthesized influx baseline equal to a
    fraction of its death flux, making it influx-pinned (non-marginal). Used to test whether a
    drug whose mechanism suppresses cell influx (e.g. MTX) can act on that cell once it has an
    influx arm. Returns True if applied (the cell had an influx param but no value)."""
    if not info["kin_param"] or info["kin_val"] is not None:
        return False
    apo = _effect_at(info["apop"], levels, prior=prior)
    mig = _effect_at(info["influx"], levels, prior=prior)      # influx flux is kIn * mig
    # influx flux = frac * death flux -> kIn * mig = frac * kd * target * apo
    info["kin_val"] = frac * info["kd_val"] * info["target"] * apo / mig
    return True


def all_regulators(info, chosen=None):
    """The union of a cell's proliferation, influx and apoptosis regulators (for candidate/truth
    scoring), restricted to ``chosen`` if given."""
    regs = set(info["prolif"]) | set(info["influx"]) | set(info["apop"])
    return {c for c in regs if chosen is None or c in chosen}

r"""Honest build order (no cheating): the agent decides the STRUCTURE from biology FIRST - which
cytokines regulate FLS proliferation - then looks up values ONLY for what it proposed, in the
real data table (MOESM2). We never hand it the answer model's parameter list. Then we grade the
structure it chose and check how close the resulting model is - does building from biology +
data lookup give a 'roughly right' model?

    python -m examples.run_qsp_build_honest         (needs ANTHROPIC_API_KEY)

Pure: reads the two JSONs; the LLM chooses regulators from the cytokine node list (with
distractors); values are looked up by (process, cytokine), not by grabbing the model's params.
"""

from __future__ import annotations

import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_tasks as LT, model_assembly as MA

_DATA = os.path.join(os.path.dirname(__file__), "..", "projects", "vantage_ra", "data")


def effect(regs, levels, maxes, ks, cap=10.0):
    s = sum((maxes[c] - 1) * levels[c] / (ks[c] + levels[c]) for c in regs)
    return 1.0 + min(cap, s)


def main() -> None:
    cfg = AgentConfig(mock=False)
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set."); return

    prov = {p["name"]: p for p in json.load(open(os.path.join(_DATA, "param_provenance.json")))}
    tg = {t["model_species"]: t for t in json.load(open(os.path.join(_DATA,
          "steady_state_targets.json"))) if t.get("model_species")}
    # cytokine node list given to the agent (includes distractors it might wrongly pick)
    cytokines = sorted(t["model_species"] for t in tg.values()
                       if t.get("kind") == "cytokine" and t.get("model_species"))
    # answer key (for GRADING only): the model's actual FLS-proliferation regulators
    truth = {n.split("Maxby")[-1] for n in prov if n.startswith("FLSProlif_Maxby")}

    print(f"== agent decides FLS-proliferation regulators from biology "
          f"(choosing among {len(cytokines)} cytokines, with distractors) ==")
    regs = MA.propose_regulators("FLS", cytokines, "proliferation", LT.default_call(cfg))
    proposed = [r["cytokine"] for r in regs]
    for r in regs:
        print(f"    {r['cytokine']:6} {r.get('direction'):4}  {r.get('basis') or ''}")

    hit = [c for c in proposed if c in truth]
    prec = len(hit) / len(proposed) if proposed else 0
    rec = len(hit) / len(truth) if truth else 0
    print(f"\n  structure vs truth ({sorted(truth)}): "
          f"precision {prec:.2f}, recall {rec:.2f}; "
          f"missed {sorted(truth - set(proposed))}, extra {sorted(set(proposed) - truth)}")

    # FAIR: keep EVERY proposed regulator. Use the real value where MOESM2 has it; for the extras
    # (proposed but not in the table) a real modeller would find literature or use a prior - so
    # assign a direction-based prior (up -> 1.5, down -> 0.6), NOT drop them. This keeps the
    # agent's over-inclusion in the model, so held-out shows its true cost.
    dir_of = {r["cytokine"]: (r.get("direction") or "up") for r in regs}
    maxes, levels, from_data, from_prior = {}, {}, [], []
    for c in proposed:
        if c not in tg:                                # cytokine has no level -> can't include
            continue
        levels[c] = float(tg[c]["target_model_unit"])
        p = prov.get(f"FLSProlif_Maxby{c}")
        if p and p.get("from_literature"):
            maxes[c] = float(p["value_from_reference"]); from_data.append(c)
        else:
            maxes[c] = 1.5 if dir_of.get(c) == "up" else 0.6   # prior, would refine from lit
            from_prior.append(c)
    print(f"\n  values: {from_data} from MOESM2; {from_prior} from a direction prior "
          "(kept, not dropped - the agent would find literature for these)")

    kept = list(maxes)
    ks = {c: levels[c] for c in kept}                  # K = level (Hill 0.5), one plausible fit
    kd = float(prov["kd_FLS_Baseline"]["value_from_reference"])
    fls_t = float(tg["FLS"]["target_model_unit"])
    eff = effect(kept, levels, maxes, ks)
    kg = fls_t * kd / eff

    # full-structure reference model (all TRUE regulators) fit the same way
    tmax = {c: float(prov[f"FLSProlif_Maxby{c}"]["value_from_reference"]) for c in truth if c in tg}
    tlev = {c: float(tg[c]["target_model_unit"]) for c in truth if c in tg}
    tks = {c: tlev[c] for c in tmax}
    teff = effect(tmax.keys(), tlev, tmax, tks); tkg = fls_t * kd / teff

    # held-out operating point (anti-IL6 therapy) - where structural errors show
    held = dict(levels); held["IL6"] = levels.get("IL6", 0) * 0.1
    theld = dict(tlev); theld["IL6"] = tlev.get("IL6", 0) * 0.1
    a_fls = kg / kd * effect(kept, held, maxes, ks)
    t_fls = tkg / kd * effect(tmax.keys(), theld, tmax, tks)
    print(f"\n== honestly-built (with the agent's extra edges kept) vs the real structure ==")
    print(f"  training (steady state): both match the target {fls_t:g} by construction")
    print(f"  held-out (anti-IL6):  agent-built = {a_fls:.3g}   real-structure = {t_fls:.3g}   "
          f"diff {abs(a_fls-t_fls)/max(a_fls,t_fls):.0%}")
    print("\n  -> the agent chose structure from biology (recall of the true regulators), kept "
          "its\n     over-inclusions with priors, and fit the rest. The held-out gap is the cost "
          "of its\n     structural choices (mostly OVER-inclusion) - not of cheating. Sloppiness "
          "keeps it small.")


if __name__ == "__main__":
    main()

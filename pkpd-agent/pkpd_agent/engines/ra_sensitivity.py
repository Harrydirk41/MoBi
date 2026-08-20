"""Stage-1 (a different skill): does the LLM know WHICH parameters matter?

Recalling biology (topology, scope) is one skill; knowing which knobs actually move the
clinical readout is another - it is dynamical-systems reasoning about leverage, not
literature recall. The paper's global sensitivity analysis (Fig 9, Sobol indices on
DAS28-CRP) gives the ground-truth top-20 most influential parameters. We hand the agent a
pool of real model parameters (the top-20 hidden among plausible distractors) and ask it
to pick/rank the most sensitive; scored by overlap with the GSA top-20, against the random
baseline (picking K of the pool blind).
"""

from __future__ import annotations

# The model's top-20 global-sensitivity parameters for DAS28-CRP, in rank order (Fig 9,
# npj Syst Biol Appl 2024). kg_/kIn_ = baseline growth/influx rates; F_ = fractional
# disease drivers; the rest are apoptosis/influx max-effect terms.
GSA_TOP20 = [
    "kg_FLS_Baseline", "F_CAM", "kg_BCells_Baseline", "kg_Macrophage_Baseline",
    "F_GMCSF", "kg_Th1_Baseline", "F_TNFa", "kIn_Th1_Baseline", "LeukoInflux_MaxbyCAM",
    "F_IL1b", "kIn_BCells_Baseline", "F_VEGF", "kIn_Macrophage_Baseline", "F_IL10",
    "F_IL17", "F_IL6", "kg_CTL_Baseline", "F_BAFF", "BCellApop_MaxbyBAFF",
    "MacroApop_MaxbyGMCSF",
]

# Real model parameters that are NOT in the top-20 - plausible-looking distractors
# (effect strengths, secretion rates, clearances, drug binding).
DISTRACTORS = [
    "kd_FLS_Baseline", "kd_Endo_Baseline", "kcl_VEGF", "kcl_RANTES", "kcl_TGFb",
    "FLSProlif_MaxbyTNFa", "FLSProlif_MaxbyIL6", "FLSProlif_MaxbyIL1b",
    "FLSProlif_MaxbyIL17", "FLSProlif_MaxbyTGFb", "VEGFSecFLS_MaxbyIL1b",
    "VEGFSecFLS_MaxbyIL6", "VEGFSecFLS_MaxbyTNFa", "IL6SecFLS_MaxbyIL1b",
    "IL6SecFLS_MaxbyIL17", "GMCSFSecFLS_MaxbyTNFa", "MIP3SecFLS_MaxbyTNFa",
    "EndoProlif_MaxbyVEGF", "Endoinflux_MaxbyTNFa", "EndoApop_MaxbyTNFa",
    "MacroProlif_MaxbyGMCSF", "MacroProlif_MaxbyTNFa", "IL1bSecFLS_MaxbyIL10",
    "TNFaSecFLS_MaxbyIL10", "RANTESSecFLS_MaxbyTNFa", "KD_ADA", "Koff_TCZ", "KD_TCZ",
    "IFNgSecTh1_MaxbyIL10", "IL17SecTh17_MaxbyIL6",
]

# Deterministic scrambled presentation order (interleave, no rank leak, no RNG).
def pool() -> list[str]:
    merged = []
    a, b = GSA_TOP20, DISTRACTORS
    for i in range(max(len(a), len(b))):
        if i < len(b):
            merged.append(b[i])
        if i < len(a):
            merged.append(a[i])
    return sorted(merged)                       # alphabetical hides the rank entirely


def _spearman(order: list[str]) -> float | None:
    """Rank correlation between the agent's order and the GSA rank, on their overlap."""
    gsa_rank = {p: i for i, p in enumerate(GSA_TOP20)}
    hits = [p for p in order if p in gsa_rank]
    if len(hits) < 3:
        return None
    agent_rank = {p: i for i, p in enumerate(hits)}
    truth_rank = {p: i for i, p in enumerate(sorted(hits, key=lambda x: gsa_rank[x]))}
    n = len(hits)
    d2 = sum((agent_rank[p] - truth_rank[p]) ** 2 for p in hits)
    return round(1 - 6 * d2 / (n * (n * n - 1)), 3)


def score_sensitivity(ranked: list[str]) -> dict:
    """Overlap of the agent's picks with the GSA top-20, vs a blind-pick random baseline."""
    top = set(GSA_TOP20)
    picks = [p for p in dict.fromkeys(ranked or [])]      # dedupe, keep order
    hit = [p for p in picks if p in top]
    n_pool = len(pool())
    prec = len(hit) / len(picks) if picks else 0.0
    rec = len(hit) / len(top)
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    # random baseline: pick len(picks) of the pool blind -> expected hits
    rand_hits = len(picks) * len(top) / n_pool if n_pool else 0.0
    rand_recall = rand_hits / len(top) if top else 0.0
    return {
        "n_picked": len(picks), "pool_size": n_pool, "hit": len(hit),
        "precision": round(prec, 3), "recall": round(rec, 3), "f1": round(f1, 3),
        "spearman_on_hits": _spearman(picks),
        "random_baseline_recall": round(rand_recall, 3),
        "beats_random": rec > rand_recall,
        "missed_top20": sorted(top - set(hit)),
        "picked_distractors": [p for p in picks if p not in top],
    }

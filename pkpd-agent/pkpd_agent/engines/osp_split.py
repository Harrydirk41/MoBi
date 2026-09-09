"""Building vs verification study split for held-out predictive grading (phase 2).

A fair test of a built model is how well it PREDICTS data it was not fitted on - exactly
how OSP qualifies its models (model-building vs model-verification). This splits a
model's observed datasets, BY WHOLE STUDY (never splitting a study across the two), into:

  * BUILDING     - given to the agent, with the concentration data, to fit on.
  * VERIFICATION - held out; the agent gets only the study DESIGN (dose/route/...), never
                   the concentrations, and is graded on how well it predicts them.

Only datasets the reference model actually maps to a simulation (the OutputMappings
linkage) are split; unmapped extras (e.g. DDI-network control arms) are dropped.

Preference order for the split:
  1. OSP's OWN tag - a SimulationClassifications group whose name marks verification
     ('model verification', 'validation', 'test'). That is the modeller's real split.
  2. A deterministic study-level holdout - IV studies stay in building (they anchor
     disposition), a spread ~1/3 of the remaining studies are held out. Reproducible.
  3. If too few studies to split meaningfully, no holdout (graded on all, flagged).
"""

from __future__ import annotations

from collections import defaultdict

from . import osp_score

_HOLDOUT_FRAC = 0.33
_MIN_STUDIES = 5
_VERIF_WORDS = ("verif", "validation", "valid.", " test", "test ", "evaluation")


def _study_key(o: dict) -> str:
    """Stable per-study id for grouping - author+year when present, else a
    route+dose signature, else the dataset name."""
    ds = o.get("dataset", "")
    return (osp_score._study_token(ds) or osp_score._study_token(o.get("study"))
            or f"{o.get('route')}|{o.get('dose')}" or ds)


def _route(o: dict) -> str:
    r = (o.get("route") or "").upper()
    return "IV" if "IV" in r else ("PO" if r in ("PO", "ORAL") else r)


def _spread_indices(n: int, k: int) -> list[int]:
    """k indices spread evenly across 0..n-1 (so a holdout samples the dose range,
    not a clustered tail)."""
    if k <= 0 or n <= 0:
        return []
    return sorted({min(n - 1, int((i + 0.5) * n / k)) for i in range(k)})


def _osp_verification_datasets(snap: dict, link: dict) -> "set[str] | None":
    """Datasets the reference's OWN SimulationClassifications mark as verification, via
    a group whose name signals it. Returns the verification dataset set, or None."""
    sim_to_ds: dict[str, list[str]] = defaultdict(list)
    for ds, sim in link.items():
        sim_to_ds[sim].append(ds)
    verif: set[str] = set()
    saw_tag = False
    for c in snap.get("SimulationClassifications") or []:
        name = (c.get("Name") or "").lower()
        if any(w in name for w in _VERIF_WORDS):
            saw_tag = True
            for sim in c.get("Classifiables") or []:
                verif.update(sim_to_ds.get(sim, []))
    return verif if (saw_tag and verif) else None


def split_studies(snap: dict, observed: list[dict],
                  holdout_frac: float = _HOLDOUT_FRAC,
                  min_studies: int = _MIN_STUDIES) -> dict:
    """Return {'building': [...], 'verification': [...], 'method': str,
    'n_studies': int, 'held_out_studies': [...]}. Dataset-name lists."""
    link = osp_score.linkage_from_snapshot(snap)
    linked = [o for o in observed if o.get("dataset") in link]
    groups: dict[str, list[str]] = defaultdict(list)
    routes: dict[str, set] = defaultdict(set)
    for o in linked:
        st = _study_key(o)
        groups[st].append(o["dataset"])
        routes[st].add(_route(o))

    all_ds = [ds for g in groups.values() for ds in g]

    # 1. OSP's own tag
    tag = _osp_verification_datasets(snap, link)
    if tag:
        verif = [ds for ds in all_ds if ds in tag]
        build = [ds for ds in all_ds if ds not in tag]
        if verif and build:
            held = sorted({s for s, g in groups.items() if any(d in tag for d in g)})
            return {"building": build, "verification": verif, "method": "osp-tag",
                    "n_studies": len(groups), "held_out_studies": held}

    # 2. deterministic study-level holdout (IV studies stay in building)
    if len(groups) >= min_studies:
        po = sorted(s for s in groups if "IV" not in routes[s])
        if po:
            k = max(1, round(len(po) * holdout_frac))
            held = {po[i] for i in _spread_indices(len(po), k)}
            # building must stay the MAJORITY by dataset count - if the held studies are
            # data-rich, return the largest to building until building >= verification.
            n_all = len(all_ds)
            while held and sum(len(groups[s]) for s in held) > n_all / 2 and len(held) > 1:
                biggest = max(held, key=lambda s: len(groups[s]))
                held.discard(biggest)
            if held and len(held) < len(groups):
                verif = [ds for s in held for ds in groups[s]]
                build = [ds for s in groups if s not in held for ds in groups[s]]
                return {"building": build, "verification": verif,
                        "method": "study-holdout", "n_studies": len(groups),
                        "held_out_studies": sorted(held)}

    # 3. too few studies to split - grade on all, flagged
    return {"building": all_ds, "verification": [], "method": "none-too-few",
            "n_studies": len(groups), "held_out_studies": []}

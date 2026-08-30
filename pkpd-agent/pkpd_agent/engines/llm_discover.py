"""Structure-discovery benchmark (B): can the agent REDISCOVER a load-bearing edge that
reading cannot supply? Ablate one regulatory edge from the calibrated model; the model then
misses an observable (a cytokine/cell level shifts). Give the LLM the SYMPTOM (what moved and
which way) plus a candidate edge set, and ask which missing edge explains it. Grade against
the edge actually removed. This tests whether data + reasoning can recover the wiring that the
literature omits - the complement to benchmark A (which tests reading).

Pure here: candidate-set construction and the LLM's symptom-to-edge reasoning (pluggable call).
The ablation + observable read live in the runner (MATLAB).
"""

from __future__ import annotations

import random

from .llm_structure import _parse_json


def candidate_set(true_edge: tuple, all_edges: list, n_distractors: int = 8,
                  seed: int = 0) -> list:
    """Build the hypothesis space: the true (removed) edge plus ``n_distractors`` other
    plausible edges, shuffled. ``all_edges`` is the pool of (src, dst) pairs to draw distractors
    from. Returns a shuffled list of (src, dst)."""
    rng = random.Random(seed)
    pool = [e for e in all_edges if e != true_edge]
    rng.shuffle(pool)
    cand = [true_edge] + pool[:n_distractors]
    rng.shuffle(cand)
    return cand


_SYS = ("You are debugging a QSP model that is missing one regulatory edge. You are given the "
        "SYMPTOM - which species levels are off, and in which direction, versus the calibrated "
        "target - and a list of candidate missing edges (source -> destination). Reason from "
        "the biology and the symptom to rank the candidates by how likely each is the missing "
        "edge. A species that is too LOW suggests a missing activator INTO it; too HIGH suggests "
        "a missing inhibitor. Output JSON only.")


def rank_candidates(symptom: list, candidates: list, call) -> list:
    """The LLM ranks candidate missing edges from the symptom. ``symptom`` is
    [{species, direction, rel_change}]; ``candidates`` is [(src, dst)]. Returns the candidates
    reordered most-likely-first (filtered to the given set). Pluggable ``call`` for tests."""
    sym = "\n".join(f"  {s['species']} is {s['direction']} "
                    f"({s.get('rel_change', '?')})" for s in symptom)
    cand = "\n".join(f"  {i+1}. {s} -> {d}" for i, (s, d) in enumerate(candidates))
    user = ("SYMPTOM (level vs calibrated target):\n" + sym +
            "\n\nCANDIDATE missing edges:\n" + cand +
            '\n\nReturn JSON {"ranking": [{"src": ..., "dst": ...}], "reason": "one phrase"} '
            "with the candidates ordered most-likely-missing first.")
    d = _parse_json(call(_SYS, user))
    cset = set(candidates)
    out = []
    for e in (d.get("ranking") or []):
        if isinstance(e, dict):
            pair = (e.get("src"), e.get("dst"))
            if pair in cset and pair not in out:
                out.append(pair)
    for pair in candidates:                            # append any the LLM dropped
        if pair not in out:
            out.append(pair)
    return out


def rank_of_true(ranking: list, true_edge: tuple) -> int:
    """1-indexed position of the true edge in the ranking (len+1 if absent)."""
    return ranking.index(true_edge) + 1 if true_edge in ranking else len(ranking) + 1

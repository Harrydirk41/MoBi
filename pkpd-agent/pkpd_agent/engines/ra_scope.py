"""Stage-2 benchmark: SCOPE selection - which cells and mediators belong in the model?

Before wiring any edges, a modeler decides the model's CAST: which cell types and soluble
mediators to include, and - just as important - which to leave out to keep the model
parsimonious. That is a judgment task (what matters for THIS disease and endpoint vs the
whole immunology textbook), and it is exactly the step the earlier tasks skipped by handing
the agent the cast. Here the agent is given only the disease and the modeling goal and must
propose the cast; we score it against the Vantage RA model's actual 26 biological nodes.

The interesting failure mode is OVER-INCLUSION: RA involves dozens of mediators, but the
model deliberately excludes many (the paper names IL-2, IL-8, IL-15, IL-18, IL-32, ...).
So recall is easy and precision is the real test - does the agent know the model's scope
discipline, or does it dump the textbook? We flag which over-inclusions are real RA
mediators the model chose to omit, so the "extra" list is interpretable.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .ra_network import CELLS, CYTOKINES, NODES, canon_node

MODEL_CELLS = list(CELLS)
MODEL_CYTOKINES = list(CYTOKINES)
MODEL_NODES = list(NODES)                       # the 26-node biological cast (answer key)

# Real RA-relevant mediators the model deliberately EXCLUDED - used only to make the
# "extra" list interpretable (these are defensible over-inclusions, not hallucinations).
_KNOWN_EXCLUDED = {
    "IL2", "IL8", "IL15", "IL18", "IL32", "IL4", "IL13", "IL21", "IL33", "IL22",
    "IL5", "IL9", "IL25", "IL34", "IL35", "IL37", "RANKL", "OPG", "CXCL9", "CXCL10",
    "CXCL13", "CCL5", "CX3CL1", "S100", "COMPLEMENT", "PGE2", "MMP", "RF", "IL23R",
    "NEUTROPHIL", "DENDRITIC", "MASTCELL", "NK", "OSTEOCLAST", "CHONDROCYTE",
}


def _norm(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(s or "")).upper()


@dataclass
class ScopeScore:
    hit: int
    missed: int
    extra: int
    precision: float
    recall: float
    f1: float
    missed_nodes: list
    extra_nodes: list
    extra_known_mediators: list


def score_scope(proposed: list[str]) -> ScopeScore:
    """Score a proposed cast (list of cell/mediator names) against the model's 26 nodes."""
    model_hits: set[str] = set()
    extras: set[str] = set()
    for raw in proposed or []:
        c = canon_node(raw)
        if c is not None:
            model_hits.add(c)
        else:
            n = _norm(raw)
            if n:
                extras.add(n)
    truth = set(MODEL_NODES)
    hit = len(model_hits & truth)
    missed = truth - model_hits
    n_prop = len(model_hits) + len(extras)
    prec = hit / n_prop if n_prop else 0.0
    rec = hit / len(truth) if truth else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    known = sorted(e for e in extras if e in _KNOWN_EXCLUDED)
    return ScopeScore(
        hit=hit, missed=len(missed), extra=len(extras),
        precision=round(prec, 3), recall=round(rec, 3), f1=round(f1, 3),
        missed_nodes=sorted(missed), extra_nodes=sorted(extras),
        extra_known_mediators=known)

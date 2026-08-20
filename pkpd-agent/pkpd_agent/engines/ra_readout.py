"""Stage-3 (the 'design' part): how does biology map to the CLINICAL READOUT?

The most creative modeling step is the bridge from mechanism to endpoint - defining
DAS28-CRP / ACR as a function of the model's physiological variables. Recovering it tests
whether the LLM can INVENT that connection, not recall a pathway. The answer key is which
model species the DAS28-CRP (and ACR) rule actually depends on, extracted by walking the
algebraic rule graph from the readout down to species (species are terminal state
variables, so the walk stops at them and never descends into the dynamics).

We then ask the agent which nodes drive DAS28-CRP and score node-recovery. Needs the full
`network.json` (sb_network_json.m dump) for the real rule graph.
"""

from __future__ import annotations

import json
import re

from .ra_network import canon_node
from .ra_scope import resolve_node

_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
DEFAULT_TARGETS = ("DAS28_CRP", "ACR_Perc", "ACR20", "DAS28")


def _rule_map(rules: list) -> dict[str, str]:
    """name -> RHS expression, for every 'LHS = RHS' rule. The LHS may carry a
    compartment prefix (e.g. 'Synovium.DAS28_CRP') - strip it so bare target names match."""
    out = {}
    for r in rules or []:
        expr = r.get("rule", "") if isinstance(r, dict) else str(r)
        if "=" in expr:
            lhs, rhs = expr.split("=", 1)
            lhs = lhs.strip().split(".")[-1]          # drop 'Synovium.' etc.
            out[lhs] = rhs
    return out


def readout_drivers(network_json: str, targets=DEFAULT_TARGETS) -> dict:
    """Which model species does the clinical readout depend on? Walk the algebraic rule
    graph from each target; collect canonical nodes (species are terminal). Returns the
    node set plus the raw target rules (for inspection) and which targets were found."""
    with open(network_json, encoding="utf-8") as fh:
        data = json.load(fh)
    rules = _rule_map(data.get("rules", []))

    found = [t for t in targets if t in rules]
    drivers: set[str] = set()
    raw = {t: rules[t].strip() for t in found}
    visited: set[str] = set()
    frontier = list(found)
    while frontier:
        name = frontier.pop()
        if name in visited:
            continue
        visited.add(name)
        rhs = rules.get(name)
        if rhs is None:
            continue
        for tok in _TOKEN.findall(rhs):
            c = canon_node(tok)
            if c is not None:
                drivers.add(c)
            elif tok in rules and tok not in visited:
                frontier.append(tok)
    return {"targets_found": found, "drivers": sorted(drivers), "target_rules": raw,
            "all_rule_names": sorted(rules)}


def score_readout(proposed: list[str], drivers: list[str]) -> dict:
    """Node-recovery of the readout drivers: precision/recall/F1 of the agent's proposed
    DAS28 drivers vs the species the readout actually depends on."""
    truth = {canon_node(d) or d for d in drivers}
    # tolerant matching for the agent's free-text (Th1 cell, endothelial cell, plasma
    # cells, ...) - the strict canon_node under-credits synonyms (as it did for scope).
    picks = {resolve_node(p) for p in (proposed or []) if resolve_node(p)}
    junk = {re.sub(r"[^A-Za-z0-9]", "", str(p)).upper()
            for p in (proposed or []) if resolve_node(p) is None}
    hit = len(picks & truth)
    false_pos = picks - truth                      # valid nodes that aren't drivers
    n_prop = len(picks) + len(junk)
    prec = hit / n_prop if n_prop else 0.0
    rec = hit / len(truth) if truth else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"n_truth": len(truth), "hit": hit, "precision": round(prec, 3),
            "recall": round(rec, 3), "f1": round(f1, 3),
            "missed": sorted(truth - picks),
            "extra": sorted(false_pos) + sorted(junk)}

"""Stage-1 topology reconstruction: draft the regulatory network from the model's NODES
plus literature, and score it against the model's OWN structure (the answer key).

The hard, creative part of model building is the wiring - which node influences which.
This benchmark isolates it: give an LLM only the node list (species) and the paper's
references, have it propose the signed influence edges (A activates / inhibits B), and
compare that draft to the edges extracted from network.json. Same discipline as
llm_structure / llm_tasks (generate a draft, regress against known-good), but for the
model STRUCTURE itself.

  * ``ground_truth_edges``  - the answer key: parse network.json into influence edges.
  * ``draft_topology``      - the LLM proposes edges from nodes (+ references).
  * ``compare_topology``    - precision / recall + which edges hit / missed / extra.

The drafter is pluggable via ``call(system, user) -> str`` (tests use a stub). Nothing
here is disease-specific: species are identified by matching the model's own node list.
"""

from __future__ import annotations

import re

from .llm_structure import _parse_json

_TOK = re.compile(r"[A-Za-z_]\w*")


def _tokens(expr: str) -> set:
    return set(_TOK.findall(expr or ""))


def ground_truth_edges(network: dict) -> set:
    """Extract the influence edges (src, dst) from a network.json dump. A species that
    appears in a reaction's rate law - directly, or via an intermediate parameter a rule
    defines from it - INFLUENCES that reaction's products; reactants flow to products.
    Returns a set of (src, dst) name pairs (unsigned). General: species are recognized by
    the model's own species list, no disease vocabulary."""
    species = {s.get("name") for s in network.get("species", []) if s.get("name")}

    def species_in(expr: str) -> set:
        return _tokens(expr) & species

    # pass 1: parameter -> the species that appear in the rule that defines it
    #         (captures the Pro_<x>_by<Y>_effect = f(Y) regulatory intermediates)
    param_srcs: dict[str, set] = {}
    for ru in network.get("rules", []):
        expr = ru.get("rule", "") if isinstance(ru, dict) else str(ru)
        if "=" not in expr:
            continue
        lhs, rhs = expr.split("=", 1)
        pname = lhs.strip().split(".")[-1]
        srcs = species_in(rhs)
        if srcs and pname and pname not in species:
            param_srcs.setdefault(pname, set()).update(srcs)

    # pass 2: reactions -> edges (rate-law modifiers + expanded params + mass flow)
    edges: set = set()
    for rx in network.get("reactions", []):
        prods = [p for p in (rx.get("products") or []) if p in species]
        reacts = [r for r in (rx.get("reactants") or []) if r in species]
        targets = prods or reacts
        rate = rx.get("rate", "")
        srcs = species_in(rate) - set(prods) - set(reacts)   # direct modifier species
        for tok in _tokens(rate):                             # params -> their source species
            if tok in param_srcs:
                srcs |= param_srcs[tok]
        for s in srcs:
            for d in targets:
                if s != d:
                    edges.add((s, d))
        for r in reacts:                                      # mass flow reactant -> product
            for p in prods:
                if r != p:
                    edges.add((r, p))
    return edges


_SYS = ("You are a systems-biology modeller. Given a list of model NODES (species) and a "
        "paper's references, propose the REGULATORY NETWORK: which node influences which, "
        "and whether it activates or inhibits. Use ONLY the exact node names given; never "
        "invent a node. Ground every edge in the biology (name the reference when you can). "
        "Output JSON only, no prose.")


def draft_topology(nodes: list, references: str, call) -> list:
    """The LLM proposes the influence edges among the given nodes (optionally drawing on the
    paper's references). Returns [{src, dst, sign, basis}], filtered to real node pairs."""
    user = ("NODES (use these EXACT names):\n" + "\n".join(nodes) +
            (("\n\nPAPER REFERENCES to draw on:\n" + references) if references else "") +
            '\n\nReturn JSON {"edges": [{"src": node, "dst": node, "sign": '
            '"activate" | "inhibit", "basis": "one phrase (+ reference)"}]}. Only edges '
            "between the nodes above.")
    d = _parse_json(call(_SYS, user))
    nset = set(nodes)
    out = []
    for e in (d.get("edges") or []):
        if not isinstance(e, dict):
            continue
        s, t = e.get("src"), e.get("dst")
        if s in nset and t in nset and s != t:
            out.append({"src": s, "dst": t, "sign": e.get("sign"),
                        "basis": e.get("basis")})
    return out


def compare_topology(draft_edges: list, truth_edges: set) -> dict:
    """Score drafted edges (unsigned) against the ground-truth influence set. Returns
    precision / recall / f1 and the hit / missed / extra edge lists - the map of what the
    literature reconstructs vs what is a modelling choice."""
    draft = {(e["src"], e["dst"]) for e in draft_edges}
    truth = set(truth_edges)
    hit = draft & truth
    p = len(hit) / len(draft) if draft else 0.0
    r = len(hit) / len(truth) if truth else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"precision": round(p, 3), "recall": round(r, 3), "f1": round(f1, 3),
            "n_draft": len(draft), "n_truth": len(truth), "hit": len(hit),
            "missed": sorted(truth - draft), "extra": sorted(draft - truth)}

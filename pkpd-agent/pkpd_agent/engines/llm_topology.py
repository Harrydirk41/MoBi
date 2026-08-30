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


def edge_provenance(network: dict) -> dict:
    """Extract the influence edges AND, for each, the model knobs that create it - the map
    from a (src, dst) edge to the rule-parameters and reactions responsible for it. This is
    what makes a functional ablation possible: to sever a regulatory edge you freeze the very
    parameter its name encodes (the ``Pro_<dst>_by<src>_effect`` intermediate).

    Returns ``{(src, dst): {"rule_params": set[str], "reactions": set[str], "mass_flow": bool}}``.
    General: species are recognized by the model's own species list; the ``_by<SRC>_`` naming is
    exploited only opportunistically (via the rule graph), never required."""
    species = {s.get("name") for s in network.get("species", []) if s.get("name")}

    def species_in(expr: str) -> set:
        return _tokens(expr) & species

    # pass 1: rule parameter -> the species in the rule that defines it (the regulatory
    #         intermediates, e.g. Pro_<x>_by<Y>_effect = f(Y))
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

    # pass 2: reactions -> edges, recording the knob(s) behind each edge
    prov: dict = {}

    def add(src, dst, rname, rule_param=None, mass=False):
        if src == dst:
            return
        e = prov.setdefault((src, dst),
                            {"rule_params": set(), "reactions": set(), "mass_flow": False})
        e["reactions"].add(rname)
        if rule_param:
            e["rule_params"].add(rule_param)
        if mass:
            e["mass_flow"] = True

    for rx in network.get("reactions", []):
        rname = rx.get("reaction") or rx.get("name") or rx.get("rate", "")
        prods = [p for p in (rx.get("products") or []) if p in species]
        reacts = [r for r in (rx.get("reactants") or []) if r in species]
        targets = prods or reacts
        rate = rx.get("rate", "")
        for s in species_in(rate) - set(prods) - set(reacts):    # direct modifier species
            for d in targets:
                add(s, d, rname)
        for tok in _tokens(rate):                                # rule params -> source species
            if tok in param_srcs:
                for s in param_srcs[tok]:
                    for d in targets:
                        add(s, d, rname, rule_param=tok)
        for r in reacts:                                         # mass flow reactant -> product
            for p in prods:
                add(r, p, rname, mass=True)
    return prov


def ground_truth_edges(network: dict) -> set:
    """The influence edge set (src, dst) - a species that appears in a reaction's rate law,
    directly or via a rule-defined intermediate, influences that reaction's products; reactants
    flow to products. Unsigned name pairs; the answer key for structural scoring."""
    return set(edge_provenance(network))


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


def functional_weights(provenance: dict, knockout, baseline: float,
                       progress=None) -> dict:
    """Weight each regulatory edge by its CLINICAL importance: knock the edge out of the REAL,
    calibrated model and measure how far the readout moves. This bridges structure to function -
    a topology error only matters if the edge it drops carries clinical signal.

    ``knockout(rule_params) -> float`` freezes the given regulatory intermediate parameters at
    their baseline (severing the edge but preserving the operating point) and returns the model's
    readout under the flagship protocol; ``baseline`` is that readout with nothing knocked out.
    Only edges carrying a rule-parameter knob are weighable (the ``_by<SRC>_`` regulatory
    intermediates); mass-flow / direct-modifier edges have no isolable knob and are skipped.
    Returns ``{(src, dst): abs(readout - baseline)}``. ``knockout`` is injected so this is
    testable without MATLAB."""
    weights: dict = {}
    weighable = [(e, p["rule_params"]) for e, p in provenance.items() if p["rule_params"]]
    for i, (edge, params) in enumerate(weighable):
        if progress:
            progress(i + 1, len(weighable), edge)
        weights[edge] = abs(knockout(sorted(params)) - baseline)
    return weights


def score_topology_functional(draft_edges: list, weights: dict) -> dict:
    """Score a drafted topology by FUNCTIONAL-WEIGHT recall, not edge count: of the total
    clinical signal carried by the (weighable) edges, how much did the draft capture? Answers
    'the LLM found 60% of edges but they account for 90% of the readout impact' (or the reverse).
    Only edges present in ``weights`` are scored (those with a measurable knockout effect)."""
    draft = {(e["src"], e["dst"]) for e in draft_edges}
    total = sum(weights.values())
    hit_edges = {e: w for e, w in weights.items() if e in draft}
    hit = sum(hit_edges.values())
    missed = sorted(((w, e) for e, w in weights.items() if e not in draft), reverse=True)
    return {"weight_recall": round(hit / total, 3) if total else 0.0,
            "edge_recall": round(len(hit_edges) / len(weights), 3) if weights else 0.0,
            "n_weighable": len(weights), "total_weight": round(total, 4),
            "hit_weight": round(hit, 4),
            "missed_ranked": [(round(w, 4), e) for w, e in missed[:20]]}

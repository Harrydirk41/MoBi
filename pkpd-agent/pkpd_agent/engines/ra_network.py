"""Stage-1 benchmark: reconstruct the RA disease NETWORK, scored against the real model.

This is the one task that is genuinely Stage-1 (model *structure*, not downstream use of
a finished model) and genuinely reasoning-heavy: the agent is given only the cast of
cells and cytokines and must propose the regulatory network - who up/down-regulates whose
secretion / proliferation / influx. There is no simulation to "run and compare"; the edge
space is combinatorial, so it cannot be brute-forced. We score the proposal against the
real Vantage RA model's own wiring, which is the answer key.

The real model encodes every regulatory edge in a naming convention:

    (Pro|Anti|Hill)_<Target><Process>[<Cell>]_by<Source>[_effect]

e.g. ``Pro_IL6Sec_byMacro_effect``  = Macrophage PROMOTES IL-6 secretion,
     ``Anti_TNFaSec_byFLS_effect``  = FLS DOWN-regulates TNF-a secretion,
     ``Hill_GMCSFSecMacro_byTNFa``  = TNF-a drives (Hill) macrophage GM-CSF secretion.

So the ground-truth edge list is recoverable from the model's parameter + rule NAMES.
``parse_truth`` works on any name list: the complete answer key comes from the MATLAB
dump ``sb_network_json.m`` (all parameters + rules); the SimBiology diagram gives a
partial bootstrap used for tests.
"""

from __future__ import annotations

import json
import re
import zipfile
from dataclasses import dataclass


# -- the cast (nodes) -------------------------------------------------------- #
# Cells and mediators the model tracks. Aliases map the model's spellings to a
# canonical node so scoring is not defeated by naming (IFNg vs IFNgamma, etc.).
CELLS = ["Macro", "Th1", "Th17", "CTL", "Treg", "BCell", "PlasmaCell", "FLS", "Endo"]
CYTOKINES = ["TNFa", "IL6", "IL17", "IL1b", "IFNg", "IL12", "IL23", "GMCSF",
             "VEGF", "BAFF", "MCP1", "MIP3", "CAM", "RANTES", "AutoAb", "TGFb", "IL10"]
NODES = CELLS + CYTOKINES

_ALIAS = {
    "macrophage": "Macro", "macrophages": "Macro", "macro": "Macro", "mac": "Macro",
    "tnf": "TNFa", "tnfa": "TNFa", "tnfalpha": "TNFa",
    "ifng": "IFNg", "ifngamma": "IFNg", "ifn": "IFNg",
    "il1b": "IL1b", "il1beta": "IL1b", "il1": "IL1b",
    "gmcsf": "GMCSF", "gm-csf": "GMCSF",
    "cd8": "CTL", "ctl": "CTL",
    "bcell": "BCell", "bcells": "BCell", "b": "BCell",
    "plasma": "PlasmaCell", "plasmacell": "PlasmaCell", "plasmacells": "PlasmaCell",
    "fls": "FLS", "endo": "Endo", "endothelial": "Endo",
    "autoab": "AutoAb", "autoantibody": "AutoAb", "mip3": "MIP3", "mcp1": "MCP1",
    "cam": "CAM", "rantes": "RANTES", "vegf": "VEGF", "baff": "BAFF",
    "tgfb": "TGFb", "tgfbeta": "TGFb", "tgf": "TGFb", "il10": "IL10",
}


def canon_node(raw: str) -> str | None:
    """Map a raw model token to a canonical node name, or None if unknown."""
    if not raw:
        return None
    k = re.sub(r"[^A-Za-z0-9]", "", raw).lower()
    if k in _ALIAS:
        return _ALIAS[k]
    for n in NODES:
        if n.lower() == k:
            return n
    return None


# process suffixes on a regulated target, longest first so 'SecMacro' strips before 'Sec'
_PROC = ["SecMacro", "SecFLS", "Sec", "Prolif", "Influx", "Apop", "Death", "Prod"]


def _split_target(core: str) -> tuple[str | None, str]:
    """'IL6Sec' -> (IL6, 'Sec');  'TNFaSecMacro' -> (TNFa, 'SecMacro')."""
    for p in _PROC:
        if core.endswith(p):
            return canon_node(core[: -len(p)]), p
    return canon_node(core), ""


@dataclass(frozen=True)
class Edge:
    source: str
    sign: int          # +1 promote, -1 inhibit
    target: str
    process: str = ""  # Sec / Prolif / Influx / ...

    def pair(self) -> tuple[str, str]:
        return (self.source, self.target)

    def signed(self) -> tuple[str, int, str]:
        return (self.source, self.sign, self.target)


_NAME_RE = re.compile(r"^(Pro|Anti|Hill)_([A-Za-z0-9]+?)_by([A-Za-z0-9]+)")
_SIGN = {"Pro": 1, "Anti": -1, "Hill": 1}


def edges_from_names(names: list[str]) -> list[Edge]:
    """Parse regulatory edges out of a list of parameter/rule names following the
    model's (Pro|Anti|Hill)_<TargetProcess>_by<Source> convention. Unknown nodes
    are dropped (keeps the edge set on the shared vocabulary)."""
    out: dict[tuple, Edge] = {}
    for nm in names:
        m = _NAME_RE.match(nm or "")
        if not m:
            continue
        sign = _SIGN[m.group(1)]
        tgt, proc = _split_target(m.group(2))
        src = canon_node(m.group(3))
        if src is None or tgt is None or src == tgt:
            continue
        e = Edge(src, sign, tgt, proc)
        out[e.signed()] = e            # dedupe on (src,sign,tgt), keep a process label
    return list(out.values())


# A regulatory rule reads  (Pro|Anti)_<TargetProcess>_effect = min(..., MM(SRC,...)+MM(SRC,...))
# where each MM()'s FIRST argument is a driving species (the real edge source). This is
# where the bulk of the network lives - the _by parameter NAMES only expose a fraction.
_RULE_RE = re.compile(r"^\s*(Pro|Anti|Hill)_([A-Za-z0-9]+?)_effect\s*=\s*(.*)$")
_MM_SRC = re.compile(r"MM\(\s*([A-Za-z_][A-Za-z0-9_]*)")


def edges_from_rules(rules: list) -> list[Edge]:
    """Parse regulatory edges from the model's repeatedAssignment rule expressions:
    each MM(source, ...) term in a Pro_/Anti_<TargetProcess>_effect rule is one edge
    source -> target with the rule's sign. This recovers the cell-process regulation
    (proliferation / influx / apoptosis / secretion drivers) that the name scan misses."""
    out: dict[tuple, Edge] = {}
    for r in rules or []:
        expr = r.get("rule", "") if isinstance(r, dict) else str(r)
        m = _RULE_RE.match(expr or "")
        if not m:
            continue
        sign = _SIGN[m.group(1)]
        tgt, proc = _split_target(m.group(2))
        if tgt is None:
            continue
        for sm in _MM_SRC.finditer(m.group(3)):
            src = canon_node(sm.group(1))
            if src is None or src == tgt:
                continue
            e = Edge(src, sign, tgt, proc)
            out.setdefault(e.signed(), e)
    return list(out.values())


def parse_truth(network_json: str) -> list[Edge]:
    """The complete answer key from the MATLAB dump (sb_network_json.m output):
    regulatory edges from the rule MM() terms UNION the _by parameter names."""
    with open(network_json, encoding="utf-8") as fh:
        s = json.load(fh)
    edges: dict[tuple, Edge] = {}
    for e in edges_from_rules(s.get("rules", [])):
        edges.setdefault(e.signed(), e)
    names = [p.get("name", "") for p in s.get("parameters", [])]
    for e in edges_from_names(names):
        edges.setdefault(e.signed(), e)
    return list(edges.values())


def parse_truth_from_diagram(sbproj: str) -> list[Edge]:
    """Partial bootstrap answer key from the SimBiology diagram (no MATLAB needed).
    The diagram draws a curated subset of the network - use for tests, not scoring."""
    z = zipfile.ZipFile(sbproj)
    diag = next(n for n in z.namelist() if n.startswith("diagram") and n.endswith(".json"))
    d = json.loads(z.read(diag))
    names = []
    for e in d["entries"][0]["content"]["entities"]:
        t = e["content"].get("title") or ""
        names.append(t)
    return edges_from_names(names)


# -- scoring ----------------------------------------------------------------- #
def score_network(proposed: list[Edge], truth: list[Edge],
                  sign_aware: bool = True) -> dict:
    """Precision / recall / F1 of a proposed edge set against the truth, plus the
    hit / missed / extra breakdown so 'extra' edges can be eyeballed for defensible
    biology the model simply omits (the known noise in full-network reconstruction).

    Reported twice by the caller: sign_aware=True (must get up-vs-down right) and
    sign_aware=False (topology only - did you find the interaction at all)."""
    def key(e: Edge):
        return e.signed() if sign_aware else e.pair()

    P = {key(e) for e in proposed}
    T = {key(e) for e in truth}
    hit = P & T
    missed = T - P
    extra = P - T
    prec = len(hit) / len(P) if P else 0.0
    rec = len(hit) / len(T) if T else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {
        "sign_aware": sign_aware,
        "n_proposed": len(P), "n_truth": len(T),
        "hit": len(hit), "missed": len(missed), "extra": len(extra),
        "precision": round(prec, 3), "recall": round(rec, 3), "f1": round(f1, 3),
        "missed_edges": sorted(str(x) for x in missed),
        "extra_edges": sorted(str(x) for x in extra),
    }


def score_signs(pred: dict, truth: list[Edge]) -> dict:
    """Isolated sign accuracy: given the TRUE (unsigned) edges, how many signs did the
    agent get right? pred maps (source, target) -> sign. Compared to the majority-class
    baseline (guess every edge the more common sign) - the real bar, since if most edges
    are activating, 'all +1' already scores high."""
    total = len(truth)
    if total == 0:
        return {"n": 0}
    correct = 0
    for e in truth:
        s = pred.get((e.source, e.target))
        if s is None:
            continue
        if (1 if s >= 0 else -1) == e.sign:
            correct += 1
    n_pos = sum(1 for e in truth if e.sign > 0)
    maj = max(n_pos, total - n_pos) / total
    acc = correct / total
    return {"n": total, "correct": correct, "accuracy": round(acc, 3),
            "majority_baseline": round(maj, 3), "beats_majority": acc > maj,
            "frac_positive": round(n_pos / total, 3)}


def edges_from_proposal(items: list[dict]) -> list[Edge]:
    """Turn an agent's proposed edges ({source, target, sign} dicts, sign as
    +1/-1 or 'promote'/'inhibit'/'activate'/'suppress') into canonical Edges."""
    out: dict[tuple, Edge] = {}
    for it in items or []:
        src = canon_node(str(it.get("source", "")))
        tgt = canon_node(str(it.get("target", "")))
        if src is None or tgt is None or src == tgt:
            continue
        s = it.get("sign", 1)
        if isinstance(s, str):
            s = -1 if s.lower() in ("-", "inhibit", "suppress", "anti", "down", "-1") else 1
        s = 1 if s >= 0 else -1
        e = Edge(src, int(s), tgt)
        out[e.signed()] = e
    return list(out.values())

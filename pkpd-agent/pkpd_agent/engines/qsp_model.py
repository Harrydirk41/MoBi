"""Model-agnostic core for the Stage-1 benchmarks: derive everything from the model dump.

The RA-specific engines (ra_network / ra_scope / ra_params / ra_readout) hardcode this
model's cast, naming, and answer keys. This module makes those DERIVED from a standardized
SimBiology dump (``sb_network_json.m`` -> network.json) plus a small per-model ``QSPModelSpec``
(patterns that say which species are drugs/PK/readouts, and the readout rule names). Point
it at a new model's network.json with a new spec and the same benchmarks run - no code
changes. Only the sensitivity task still needs an external GSA list (a figure/analysis, not
derivable from structure).

What is derived here:
  * nodes      - biological species (all species minus drug/PK/readout by the spec patterns)
  * edges      - signed regulatory edges from the rule expressions (Pro/Anti + MM() drivers)
  * params     - name/units/value, split dimensionless / physiological / model-scaling
  * readout    - the species the clinical-readout rule depends on (rule-graph walk)
  * a model-scoped tolerant matcher for the agent's free-text proposals

Cross-model generality is BUILT here but validated only on this RA model (backward-compat
regression); a second QSP model is needed to prove it truly transfers.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field

from .ra_network import Edge, _SIGN, _split_target  # reuse the edge dataclass + helpers
from .ra_params import _MODEL_SCALING_UNITS

_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_RULE_RE = re.compile(r"^\s*(Pro|Anti|Hill)_([A-Za-z0-9]+?)_effect\s*=\s*(.*)$")
_MM_SRC = re.compile(r"MM\(\s*([A-Za-z_][A-Za-z0-9_]*)")
_NAME_RE = re.compile(r"^(Pro|Anti|Hill)_([A-Za-z0-9]+?)_by([A-Za-z0-9]+)")
# unprefixed effect-strength params also name edges: 'IL6SecMacro_MaxbyTNFa' -> TNFa->IL6
_NAME_RE2 = re.compile(r"^([A-Za-z0-9]+?)_(?:Max)?by([A-Za-z0-9]+)$")


@dataclass
class QSPModelSpec:
    """The only per-model configuration needed to run the benchmarks on a new model."""
    name: str
    readout_targets: list[str]                         # rule LHS names of the clinical readout
    drug_patterns: list[str] = field(default_factory=list)     # species regexes = drugs/PK
    readout_patterns: list[str] = field(default_factory=list)  # species regexes = readouts/flags
    aliases: dict = field(default_factory=dict)        # extra free-text synonyms -> node
    readout_name: str = "the disease-severity score"   # what the readout is, for prompts
    gsa_top: list[str] = field(default_factory=list)   # global-sensitivity top params
                                                       #   (from an analysis/figure, external)


def _norm(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(s or "")).upper()


@dataclass
class QSPParam:
    name: str
    units: str
    value: float

    def dimensionless(self) -> bool:
        return self.units.lower() in ("dimensionless", "", "none", "fraction")

    def model_scaling(self) -> bool:
        return self.units.lower() in _MODEL_SCALING_UNITS

    def physiological(self) -> bool:
        return (not self.dimensionless()) and (not self.model_scaling())


class QSPModel:
    """A model derived from network.json + a QSPModelSpec. Provides the answer keys and a
    model-scoped tolerant matcher, so the benchmarks are model-agnostic."""

    def __init__(self, data: dict, spec: QSPModelSpec):
        self.spec = spec
        self.species = [s["name"] for s in data.get("species", [])]
        drug = [re.compile(p, re.IGNORECASE) for p in spec.drug_patterns]
        rdt = [re.compile(p, re.IGNORECASE) for p in spec.readout_patterns]
        self.nodes = [s for s in self.species
                      if not any(p.search(s) for p in drug + rdt)]
        self._norm2node = {_norm(n): n for n in self.nodes}
        # plural/singular tolerance + spec aliases, all normalized
        self._alias = {}
        for n in self.nodes:
            k = _norm(n)
            self._alias.setdefault(k.rstrip("S"), n)      # BCells -> BCELL
            self._alias.setdefault(k + "S", n)
        for a, n in (spec.aliases or {}).items():
            self._alias[_norm(a)] = n

        rules = data.get("rules", [])
        self.rule_rhs = self._rule_map(rules)
        self.edges = self._edges(rules, data.get("parameters", []))
        self.params = [QSPParam(p["name"], str(p.get("units") or "dimensionless"),
                                float(p["value"]))
                       for p in data.get("parameters", [])
                       if _isnum(p.get("value"))]
        self.readout_drivers = self._readout()

    # -- matching -------------------------------------------------------- #
    def canon(self, raw: str):
        n = _norm(raw)
        if n in self._norm2node:
            return self._norm2node[n]
        return None

    def resolve(self, raw: str):
        """Tolerant: exact node, then alias (plural/singular/spec), then 'cell(s)' strip."""
        c = self.canon(raw)
        if c:
            return c
        n = _norm(raw)
        if n in self._alias:
            return self._alias[n]
        for suf in ("CELLS", "CELL"):
            if n.endswith(suf) and len(n) > len(suf):
                b = n[: -len(suf)]
                if b in self._norm2node:
                    return self._norm2node[b]
                if b in self._alias:
                    return self._alias[b]
        return None

    # -- derivation ------------------------------------------------------ #
    def _rule_map(self, rules) -> dict:
        out = {}
        for r in rules or []:
            expr = r.get("rule", "") if isinstance(r, dict) else str(r)
            if "=" in expr:
                lhs, rhs = expr.split("=", 1)
                out[lhs.strip().split(".")[-1]] = rhs
        return out

    def _edges(self, rules, params) -> list[Edge]:
        # rule/param names use abbreviations (EndoInflux, MacroProlif) while species use
        # full names (Endothelial, Macrophages), so match tolerantly (aliases), not strict.
        out: dict[tuple, Edge] = {}
        for r in rules or []:
            expr = r.get("rule", "") if isinstance(r, dict) else str(r)
            m = _RULE_RE.match(expr or "")
            if not m:
                continue
            sign = _SIGN[m.group(1)]
            tgt, proc = _split_target_general(m.group(2), self.resolve)
            if tgt is None:
                continue
            for sm in _MM_SRC.finditer(m.group(3)):
                src = self.resolve(sm.group(1))
                if src and src != tgt:
                    out.setdefault((src, sign, tgt), Edge(src, sign, tgt, proc))
        for p in params or []:
            nm = p.get("name", "") or ""
            m = _NAME_RE.match(nm)
            if m:
                sign = _SIGN[m.group(1)]
                tgt, proc = _split_target_general(m.group(2), self.resolve)
                src = self.resolve(m.group(3))
            else:                                    # unprefixed effect-strength params
                m2 = _NAME_RE2.match(nm)             # 'IL6SecMacro_MaxbyTNFa' -> TNFa->IL6
                if not m2:
                    continue
                sign = 1
                tgt, proc = _split_target_general(m2.group(1), self.resolve)
                src = self.resolve(m2.group(2))
            if src and tgt and src != tgt:
                out.setdefault((src, sign, tgt), Edge(src, sign, tgt, proc))
        return list(out.values())

    def _readout(self) -> list[str]:
        drivers: set[str] = set()
        visited: set[str] = set()
        frontier = [t for t in self.spec.readout_targets]
        while frontier:
            name = frontier.pop()
            if name in visited:
                continue
            visited.add(name)
            rhs = self.rule_rhs.get(name)
            if rhs is None:
                continue
            for tok in _TOKEN.findall(rhs):
                c = self.canon(tok)
                if c:
                    drivers.add(c)
                elif tok in self.rule_rhs and tok not in visited:
                    frontier.append(tok)
        return sorted(drivers)

    # -- model-scoped scoring (uses this model's matcher, not RA's) ------ #
    def edges_from_proposal(self, items) -> list:
        out: dict[tuple, Edge] = {}
        for it in items or []:
            src = self.resolve(str(it.get("source", "")))
            tgt = self.resolve(str(it.get("target", "")))
            if src is None or tgt is None or src == tgt:
                continue
            s = it.get("sign", 1)
            if isinstance(s, str):
                s = -1 if s.lower() in ("-", "inhibit", "suppress", "anti", "down", "-1") else 1
            e = Edge(src, 1 if s >= 0 else -1, tgt)
            out[e.signed()] = e
        return list(out.values())

    def resolve_all(self, raw: str) -> set:
        """Resolve a free-text entry that may pack several nodes or annotations into one
        string ('B cell / plasma cell', 'Macrophage (synovial)', 'CD4 T cell (Th1, Th17)').
        Returns the set of nodes it references (empty if none)."""
        r = self.resolve(raw)
        if r:
            return {r}
        out = set()
        for part in re.split(r"[\/,;+]|\(|\)|\band\b", str(raw)):
            c = self.resolve(part.strip())
            if c:
                out.add(c)
        return out

    def score_node_set(self, proposed: list, truth: list) -> dict:
        """Generic precision/recall/F1 of proposed node names vs a truth node set. Handles
        compound/annotated free-text entries by splitting them before matching."""
        picks, junk = set(), set()
        for p in (proposed or []):
            got = self.resolve_all(p)
            if got:
                picks |= got
            else:
                junk.add(_norm(p))
        T = set(truth)
        hit = len(picks & T)
        n_prop = len(picks) + len(junk)
        prec = hit / n_prop if n_prop else 0.0
        rec = hit / len(T) if T else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        return {"n_truth": len(T), "hit": hit, "precision": round(prec, 3),
                "recall": round(rec, 3), "f1": round(f1, 3),
                "missed": sorted(T - picks),
                "extra": sorted(picks - T) + sorted(junk)}

    # -- sensitivity (needs the spec's external GSA list) ---------------- #
    def sensitivity_pool(self, n_distractors: int = 30) -> list[str]:
        """The GSA top params hidden among real distractor param names from this model."""
        top = list(self.spec.gsa_top)
        pnames = [p.name for p in self.params if p.name not in set(top)]
        # deterministic, evenly-spaced distractor sample (no RNG)
        if len(pnames) > n_distractors:
            step = len(pnames) / n_distractors
            pnames = [pnames[int(i * step)] for i in range(n_distractors)]
        return sorted(top + pnames)

    def score_sensitivity(self, ranked: list) -> dict:
        top = set(self.spec.gsa_top)
        # the task is 'rank the ~top-K most sensitive'; keep only the agent's top-K so
        # dumping the whole pool cannot trivially score recall 1.0.
        picks = [p for p in dict.fromkeys(ranked or [])][: max(1, len(top))]
        hit = [p for p in picks if p in top]
        pool_n = len(self.sensitivity_pool())
        prec = len(hit) / len(picks) if picks else 0.0
        rec = len(hit) / len(top) if top else 0.0
        rand = (len(picks) * len(top) / pool_n / len(top)) if (pool_n and top) else 0.0
        return {"n_picked": len(picks), "hit": len(hit),
                "precision": round(prec, 3), "recall": round(rec, 3),
                "random_baseline_recall": round(rand, 3), "beats_random": rec > rand,
                "missed_top": sorted(top - set(hit))}

    @classmethod
    def from_network_json(cls, path: str, spec: QSPModelSpec) -> "QSPModel":
        with open(path, encoding="utf-8") as fh:
            return cls(json.load(fh), spec)

    @classmethod
    def from_sbml(cls, path: str, spec: QSPModelSpec = None, name: str = "QSP model"):
        """Build directly from a standard SBML export - no MATLAB. If no spec is given the
        spec is heuristically inferred from the parsed structure."""
        from .sbml_import import sbml_to_network
        data = sbml_to_network(path)
        return cls(data, spec or infer_spec(data, name))

    @classmethod
    def inferred(cls, path: str, name: str = "QSP model") -> "QSPModel":
        """Build with a HEURISTICALLY inferred spec - no hand config. Best-effort: it can
        over-include drug-conjugate species that carry no generic token, and cannot supply
        the external GSA list (sensitivity is skipped). Validate before trusting."""
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return cls(data, infer_spec(data, name))


def _split_target_general(core: str, canon):
    """Strip a process suffix and canon the target, using the model's own matcher."""
    for p in ("SecMacro", "SecFLS", "Sec", "Prolif", "Influx", "Apop", "Death", "Prod"):
        if core.endswith(p):
            return canon(core[: -len(p)]), p
    return canon(core), ""


def _isnum(x) -> bool:
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


# Per-model config is DATA, not code: each project ships projects/<name>/spec.json. No
# model's specifics are hardcoded in the engine - this just loads them.
_PROJECTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             "..", "..", "projects"))


def spec_from_dict(d: dict) -> QSPModelSpec:
    return QSPModelSpec(
        name=d.get("name", "QSP model"),
        readout_targets=list(d.get("readout_targets", [])),
        drug_patterns=list(d.get("drug_patterns", [])),
        readout_patterns=list(d.get("readout_patterns", [])),
        aliases=dict(d.get("aliases", {})),
        readout_name=d.get("readout_name", "the disease-severity score"),
        gsa_top=list(d.get("gsa_top", [])))


def load_spec(name: str, projects_dir: str = None) -> QSPModelSpec:
    """Load a project's spec.json into a QSPModelSpec. `name` is the project folder."""
    base = projects_dir or _PROJECTS_DIR
    path = os.path.join(base, name, "spec.json")
    with open(path, encoding="utf-8") as fh:
        return spec_from_dict(json.load(fh))

_READOUT_TOKENS = ["DAS", "ACR", "SCORE", "ACTIVITY", "REMISSION", "RESPONSE",
                   "NONRESP", "DELTA", "PASI", "SLEDAI", "EULAR"]
# PK-compartment tokens. Deliberately NOT 'PLASMA' (collides with PlasmaCells) - the
# residual drug-conjugate species (e.g. TNFa_Ada) that carry no generic token are the
# known blind spot of heuristic inference, and are what an LLM classifier would catch.
_DRUG_TOKENS = ["DRUG", "DOSE", "CENTRAL", "PERIPHERAL", "AVAILABLE", "DEPOT", "GI",
                "SUBCUT"]


def _auto_aliases(rules: list, params: list, species: list) -> dict:
    """Recover naming abbreviations (Endo->Endothelial, Macro->Macrophages) by matching
    the tokens that appear in rule/param names against the species by normalized prefix."""
    snorm = {_norm(s): s for s in species}
    toks: set[str] = set()
    for r in rules or []:
        expr = r.get("rule", "") if isinstance(r, dict) else str(r)
        for m in re.finditer(r"(?:Pro|Anti|Hill)_([A-Za-z0-9]+?)_(?:effect|by([A-Za-z0-9]+))",
                             expr or ""):
            toks.add(m.group(1))
            if m.group(2):
                toks.add(m.group(2))
        for m in _MM_SRC.finditer(expr or ""):
            toks.add(m.group(1))
    for p in params or []:
        m = _NAME_RE.match(p.get("name", "") or "")
        if m:
            toks.add(m.group(2)); toks.add(m.group(3))
    out = {}
    for t in toks:
        # strip a trailing process word so 'FLSProlif' -> 'FLS'
        core = t
        for suf in ("SecMacro", "SecFLS", "Sec", "Prolif", "Influx", "Apop", "Death", "Prod"):
            if core.endswith(suf):
                core = core[: -len(suf)]
                break
        n = _norm(core)
        if not n or n in snorm:
            continue
        # prefix match to a species (Endo -> Endothelial); shortest species wins
        cands = sorted((s for k, s in snorm.items() if k.startswith(n)), key=len)
        if cands:
            out[n.lower()] = cands[0]
    return out


def infer_spec(data: dict, name: str = "QSP model") -> QSPModelSpec:
    """Heuristically infer a QSPModelSpec from a network.json dump - no hand config.
    Classifies species into readout/drug/biology by token, finds the readout rules, and
    auto-derives naming aliases. The one thing it cannot invent is the external GSA list."""
    species = [s["name"] for s in data.get("species", [])]
    rules = data.get("rules", [])
    rule_lhs = []
    for r in rules or []:
        expr = r.get("rule", "") if isinstance(r, dict) else str(r)
        if "=" in expr:
            rule_lhs.append(expr.split("=", 1)[0].strip().split(".")[-1])

    # readout rules: the composite score, matched by PREFIX (so 'ACR' hits 'ACR_Perc'
    # but not 'Macro'); the derivation walks intermediates like delta_ from there.
    readout_targets = [l for l in rule_lhs
                       if any(l.upper().startswith(t) for t in ["DAS", "ACR", "SCORE"])]

    def _pat(tok: str) -> str:
        # short/ambiguous tokens (ACR, DAS, GI) are anchored at a word boundary to avoid
        # matching inside a species name (ACR in Macro); longer tokens can be substrings.
        return rf"(^|_){tok}" if len(tok) <= 3 else tok

    drug_pat = sorted({_pat(t) for s in species for t in _DRUG_TOKENS if t in s.upper()})
    rdt_pat = sorted({_pat(t) for s in species for t in _READOUT_TOKENS
                      if re.search(rf"(^|_){t}" if len(t) <= 3 else t, s, re.IGNORECASE)})
    return QSPModelSpec(
        name=name, readout_targets=readout_targets or ["DAS28_CRP"],
        drug_patterns=drug_pat, readout_patterns=rdt_pat,
        aliases=_auto_aliases(rules, data.get("parameters", []), species))


# Add a new model = add a projects/<name>/ folder with a spec.json. Short aliases map to
# project folder names; no spec literal lives in code.
_SPEC_ALIASES = {"ra": "vantage_ra", "vantage_ra": "vantage_ra"}


def get_spec(name: str, projects_dir: str = None) -> QSPModelSpec:
    key = (name or "ra").lower()
    folder = _SPEC_ALIASES.get(key, key)
    base = projects_dir or _PROJECTS_DIR
    if not os.path.isdir(os.path.join(base, folder)):
        known = sorted(d for d in os.listdir(base)
                       if os.path.isdir(os.path.join(base, d))) if os.path.isdir(base) else []
        raise KeyError(f"unknown project '{name}'. Known: {known}. "
                       "Add a projects/<name>/spec.json.")
    return load_spec(folder, base)


# Backward-compat handle: the RA spec, now loaded from data rather than a code literal.
VANTAGE_RA_SPEC = load_spec("vantage_ra")

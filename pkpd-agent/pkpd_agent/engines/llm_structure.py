"""LLM structure extractor: read ANY model's naming, output the canonical schema.

This removes the brittle-naming limitation of the regex derivation (which only understands
this model's ``Pro_/Anti_/MM()`` convention). An LLM reads the model's species + rule
expressions and returns, in a fixed schema: which species are biology vs drug/readout, the
signed regulatory edges, and the readout drivers - whatever naming the author used.

The critical safety rule for a BENCHMARK: the extracted structure becomes the ANSWER KEY,
so an extraction error corrupts the key. Two guards, both here:
  1. The extractor LLM is SEPARATE from the benchmarked LLM (different call, different job).
  2. ``compare_to_model`` regresses the extractor against the deterministic parser on a
     model whose convention we DO understand (RA): the extractor must reproduce its edges /
     readout drivers before we trust it on a new-convention model. Trust the LLM to READ
     arbitrary naming; trust the regression to prove it read correctly.

The extractor is pluggable via ``call_fn(system, user) -> str`` so tests use a stub and
real runs use Claude. Extraction is low-ambiguity (reading a formula), which is why an LLM
does it reliably - but it is verified, never assumed.
"""

from __future__ import annotations

import json
import re


def _parse_json(text: str):
    """Pull the JSON object/array out of an LLM reply (tolerating code fences / prose)."""
    t = re.sub(r"^```[a-zA-Z]*\n?|```$", "", text.strip()).strip()
    try:
        return json.loads(t)                       # clean JSON (object or array) - common
    except json.JSONDecodeError:
        pass
    # else find the EARLIEST top-level bracket (not a nested one) and balance from there
    starts = sorted((t.find(o), o, c) for o, c in (("{", "}"), ("[", "]")) if t.find(o) >= 0)
    for i, opener, closer in starts:
        depth = 0
        for j in range(i, len(t)):
            if t[j] == opener:
                depth += 1
            elif t[j] == closer:
                depth -= 1
                if depth == 0:
                    try:
                        return json.loads(t[i:j + 1])
                    except json.JSONDecodeError:
                        break
    return json.loads(t)          # last resort (raises with a clear error)


def default_call(config):
    """A real Claude call: (system, user) -> text. Uses the configured model."""
    import anthropic
    client = anthropic.Anthropic()

    def call(system: str, user: str) -> str:
        resp = client.messages.create(
            model=config.model, max_tokens=8000, system=system,
            messages=[{"role": "user", "content": user}])
        return "".join(getattr(b, "text", "") for b in resp.content
                       if getattr(b, "type", None) == "text")
    return call


_SYS = ("You are a systems-biology model parser. You read a QSP model's species and rule "
        "expressions and extract its STRUCTURE into strict JSON. You do not invent biology "
        "- you report only what the equations encode. Map any abbreviations in rule names "
        "to the exact species names given. Output JSON only, no prose.")


def classify_species(species: list[str], call) -> dict:
    """Split species into biology / drug-PK / readout using the LLM (replaces the token
    heuristic). Returns {'biology': [...], 'drug': [...], 'readout': [...]}."""
    user = (
        "Classify each species as one of: 'biology' (a cell type or a soluble mediator that "
        "participates in the disease dynamics), 'drug' (a drug/PK compartment or drug-bound "
        "species), or 'readout' (a clinical score, response flag, or derived output). "
        'Return JSON {"biology": [...], "drug": [...], "readout": [...]} using the exact '
        f"names.\n\nSPECIES:\n{json.dumps(species)}")
    r = _parse_json(call(_SYS, user))
    return {k: list(r.get(k, [])) for k in ("biology", "drug", "readout")}


def extract_edges(bio_species: list[str], rules: list[str], call, batch: int = 25) -> list:
    """Signed regulatory edges from the rule expressions. Each edge {source, target, sign}
    uses biology species names; source drives the target's process, sign +1 activating /
    -1 inhibiting."""
    out, seen = [], set()
    for k in range(0, len(rules), batch):
        chunk = rules[k:k + batch]
        user = (
            "From each rule, extract the regulatory edges it encodes, using this EXACT edge "
            "convention (match it precisely):\n"
            "- A rule governs a TARGET species' PROCESS. The process is named in the rule's "
            "left-hand side: Prolif/Growth/Influx/Apop for a CELL, or Sec (secretion) for a "
            "MEDIATOR.\n"
            "- SECRETION rule for mediator C (e.g. 'C_Sec...' or 'CSecK...'): emit the edge "
            "from the PRODUCING CELL K (named in the LHS, e.g. 'IL6SecFLS' -> FLS) to C, AND "
            "an edge from each DRIVER species D appearing in the expression (the first arg of "
            "each MM()/Hill term) to C.\n"
            "- CELL-process rule for cell K (Prolif/Influx/Apop of K): emit an edge from each "
            "DRIVER D to K.\n"
            "- Prefer these cell<->mediator edges; do NOT invent cytokine->cytokine shortcuts "
            "unless a cytokine literally appears as the MM() driver of another mediator's "
            "secretion rule.\n"
            "- sign = +1 if the rule increases the target (Pro/promoting), -1 if it decreases "
            "it (Anti/inhibiting).\n"
            "Use ONLY these biology species names; map abbreviations in the rule name to "
            "them. Ignore drug/readout terms.\n"
            'Return JSON [{"source": "...", "target": "...", "sign": 1 or -1}, ...].\n\n'
            f"BIOLOGY SPECIES:\n{json.dumps(bio_species)}\n\nRULES:\n"
            + "\n".join(chunk))
        for e in _parse_json(call(_SYS, user)):
            src, tgt = e.get("source"), e.get("target")
            sgn = 1 if int(e.get("sign", 1)) >= 0 else -1
            if src in bio_species and tgt in bio_species and src != tgt:
                key = (src, sgn, tgt)
                if key not in seen:
                    seen.add(key)
                    out.append({"source": src, "sign": sgn, "target": tgt})
    return out


def extract_readout(bio_species: list[str], readout_rules: list[str], call) -> list:
    """Which biology species the clinical-readout rule(s) are a direct function of."""
    if not readout_rules:
        return []
    user = (
        "These rules define the clinical readout. List the biology species the readout is a "
        "direct function of (the variables it sums/combines). Use the exact species names.\n"
        'Return JSON {"drivers": [...]}.\n\nBIOLOGY SPECIES:\n'
        f"{json.dumps(bio_species)}\n\nREADOUT RULES:\n" + "\n".join(readout_rules))
    r = _parse_json(call(_SYS, user))
    return [d for d in r.get("drivers", []) if d in bio_species]


def _rule_strings(data: dict) -> list[str]:
    return [r.get("rule", "") if isinstance(r, dict) else str(r)
            for r in data.get("rules", [])]


def extract_structure(data: dict, call, readout_hint: list[str] = None) -> dict:
    """Full LLM extraction of the canonical schema from a network.json dict."""
    species = [s["name"] for s in data.get("species", [])]
    rules = _rule_strings(data)
    cls = classify_species(species, call)
    bio = cls["biology"]
    # readout rules: those whose LHS names a readout species, else the hinted targets
    readout_names = set(cls["readout"]) | set(readout_hint or [])
    readout_rules = [r for r in rules
                     if "=" in r and r.split("=", 1)[0].strip().split(".")[-1] in readout_names]
    return {
        "classification": cls,
        "nodes": bio,
        "edges": extract_edges(bio, rules, call),
        "readout_drivers": extract_readout(bio, readout_rules, call),
    }


def compare_to_model(result: dict, model) -> dict:
    """Regression: does the LLM extraction reproduce the deterministic QSPModel's keys?
    Reports precision/recall of edges and readout drivers, and node-set agreement."""
    def prf(pred: set, truth: set) -> dict:
        hit = len(pred & truth)
        p = hit / len(pred) if pred else 0.0
        r = hit / len(truth) if truth else 0.0
        f = 2 * p * r / (p + r) if (p + r) else 0.0
        return {"precision": round(p, 3), "recall": round(r, 3), "f1": round(f, 3),
                "hit": hit, "n_pred": len(pred), "n_truth": len(truth)}

    llm_e = {(e["source"], e["sign"], e["target"]) for e in result["edges"]}
    tru_e = {e.signed() for e in model.edges}
    llm_e_topo = {(s, t) for s, _, t in llm_e}
    tru_e_topo = {(s, t) for s, _, t in tru_e}
    return {
        "edges_signed": prf(llm_e, tru_e),
        "edges_topology": prf(llm_e_topo, tru_e_topo),
        "readout": prf(set(result["readout_drivers"]), set(model.readout_drivers)),
        "nodes": prf(set(result["nodes"]), set(model.nodes)),
        "edge_disagreements": {
            "llm_only": sorted(str(x) for x in (llm_e_topo - tru_e_topo))[:15],
            "deterministic_only": sorted(str(x) for x in (tru_e_topo - llm_e_topo))[:15]},
    }

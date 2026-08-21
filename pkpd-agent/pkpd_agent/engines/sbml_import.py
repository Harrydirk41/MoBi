"""Read a model from standard SBML - no MATLAB, no bespoke dumper.

SimBiology (and every other tool) can export SBML, so parsing SBML directly removes the
dependency on our MATLAB ``sb_network_json.m``. The one non-trivial step is converting
SBML's MathML (used for rule and kinetic-law expressions) back into the infix strings the
QSPModel derivation expects (``Pro_FLSProlif_effect = min(10, MM(TNFa,...)+...)``). This
handles the operators QSP models actually use (+ - * / ^, min/max, and user-function
application like MM(...)); anything exotic degrades gracefully.

    from pkpd_agent.engines.sbml_import import sbml_to_network
    data = sbml_to_network("model.xml")          # same dict shape as network.json
    model = QSPModel(data, spec)                  # ...then everything else is unchanged

Honest limits: SBML encodes units by reference to unit definitions, so the parameter
unit STRINGS may need a per-model unit-name map for the param benchmark's unit split; the
structural benchmarks (scope/topology/signs/readout) need no units and work as-is.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

_BINOPS = {"plus": "+", "minus": "-", "times": "*", "divide": "/", "power": "^"}


def _local(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


def mathml_to_infix(e) -> str:
    """Convert a MathML element to an infix expression string."""
    tag = _local(e.tag)
    if tag in ("cn", "ci", "csymbol"):
        return (e.text or "").strip()
    if tag == "apply":
        ch = list(e)
        if not ch:
            return ""
        op, args = ch[0], [mathml_to_infix(x) for x in ch[1:]]
        ot = _local(op.tag)
        if ot == "ci":                                   # user-function application: MM(...)
            return f"{(op.text or '').strip()}({','.join(args)})"
        if ot in _BINOPS:
            if ot == "minus" and len(args) == 1:
                return f"(-{args[0]})"
            return "(" + f" {_BINOPS[ot]} ".join(args) + ")"
        return f"{ot}({','.join(args)})"                 # min/max/root/exp/ln/...
    if tag in ("math", "lambda"):
        parts = [mathml_to_infix(x) for x in e if _local(x.tag) not in ("bvar",)]
        return parts[-1] if parts else ""
    # unknown wrapper: concatenate children
    return "".join(mathml_to_infix(x) for x in e)


def _find(parent, name):
    return [c for c in parent.iter() if _local(c.tag) == name]


def sbml_to_network(path: str) -> dict:
    """Parse an SBML file into the network.json dict shape (species/reactions/rules/params)."""
    root = ET.parse(path).getroot()
    model = next((c for c in root.iter() if _local(c.tag) == "model"), root)

    def _list(container_name, item_name):
        conts = [c for c in model.iter() if _local(c.tag) == container_name]
        out = []
        for cont in conts:
            out += [c for c in cont if _local(c.tag) == item_name]
        return out

    species = [{"name": s.get("id") or s.get("name")}
               for s in _list("listOfSpecies", "species")]
    parameters = [{"name": p.get("id") or p.get("name"),
                   "value": _num(p.get("value")),
                   "units": p.get("units") or "dimensionless"}
                  for p in _list("listOfParameters", "parameter")]

    rules = []
    for r in model.iter():
        t = _local(r.tag)
        if t in ("assignmentRule", "rateRule"):
            var = r.get("variable") or ""
            math = next((c for c in r if _local(c.tag) == "math"), None)
            if var and math is not None:
                rules.append({"type": "repeatedAssignment" if t == "assignmentRule"
                              else "rate", "rule": f"{var} = {mathml_to_infix(math)}"})

    reactions = []
    for rx in _list("listOfReactions", "reaction"):
        reac = [sr.get("species") for sr in rx.iter()
                if _local(sr.tag) == "speciesReference"
                and _local(_parent_tag(rx, sr)) == "listOfReactants"]
        prod = [sr.get("species") for sr in rx.iter()
                if _local(sr.tag) == "speciesReference"
                and _local(_parent_tag(rx, sr)) == "listOfProducts"]
        kl = next((c for c in rx.iter() if _local(c.tag) == "kineticLaw"), None)
        rate = ""
        if kl is not None:
            math = next((c for c in kl if _local(c.tag) == "math"), None)
            if math is not None:
                rate = mathml_to_infix(math)
        reactions.append({"name": rx.get("id") or rx.get("name") or "",
                          "reaction": " + ".join(reac) + " -> " + " + ".join(prod),
                          "rate": rate,
                          "reactants": [r for r in reac if r],
                          "products": [p for p in prod if p]})

    return {"species": species, "parameters": parameters, "rules": rules,
            "reactions": reactions}


def _parent_tag(root, child):
    for p in root.iter():
        for c in p:
            if c is child:
                return p.tag
    return ""


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")

"""Stage-1 model ASSEMBLY (steps 3+4+7): turn a topology + a rate-law motif per process +
known parameter values into a runnable SBML model - the piece that was never built, tested
end-to-end WITHOUT reading any paper (parameters are supplied, not extracted).

  * ``combined_effect``  - the library motif: a process rate modulated by its regulators as a
                           product of fold-change Hill terms (the convention this model uses).
  * ``build_subsystem``  - assemble a self-contained runnable subsystem (species + reactions +
                           rules) from edges + params, with the cytokine inputs clamped.
  * ``to_sbml``          - emit standard SBML (MathML kinetic laws / rules) that SimBiology can
                           import and that our own sbml_import round-trips - the assembly output.

General: no disease vocabulary; the motif and wiring come from the passed spec. Params are
inputs here (the "assume parameters known" path), so this measures the ASSEMBLY, not reading.
"""

from __future__ import annotations

import re
import xml.sax.saxutils as _sx

_TOKEN = re.compile(r"\s*(\d+\.?\d*(?:[eE][-+]?\d+)?|[A-Za-z_]\w*|[()+\-*/^,]|\S)")


def _tokenize(expr):
    pos, out = 0, []
    for m in _TOKEN.finditer(expr):
        out.append(m.group(1))
    return out


def infix_to_mathml(expr: str) -> str:
    """Convert an infix rate/rule expression to MathML (the inverse of sbml_import's reader).
    Handles + - * / ^, unary minus, function application f(a,b), numbers and identifiers - the
    operators QSP rate laws use. Precedence: ^ > * / > + -."""
    toks = _tokenize(expr)
    i = 0

    def peek():
        return toks[i] if i < len(toks) else None

    def eat(t=None):
        nonlocal i
        tok = toks[i]; i += 1
        return tok

    def prim():
        nonlocal i
        t = peek()
        if t == "(":
            eat(); e = add(); eat()  # ')'
            return e
        if t == "-":
            eat(); return ("apply", "minus", [prim()])
        if re.match(r"^\d", t):
            eat(); return ("cn", t)
        name = eat()
        if peek() == "(":                                # function application
            eat(); args = [add()]
            while peek() == ",":
                eat(); args.append(add())
            eat()  # ')'
            return ("fn", name, args)
        return ("ci", name)

    def power():
        left = prim()
        while peek() == "^":
            eat(); left = ("apply", "power", [left, prim()])
        return left

    def mul():
        left = power()
        while peek() in ("*", "/"):
            op = "times" if eat() == "*" else "divide"
            left = ("apply", op, [left, power()])
        return left

    def add():
        left = mul()
        while peek() in ("+", "-"):
            op = "plus" if eat() == "+" else "minus"
            left = ("apply", op, [left, mul()])
        return left

    def emit(node):
        k = node[0]
        if k == "cn":
            return f"<cn>{node[1]}</cn>"
        if k == "ci":
            return f"<ci>{_sx.escape(node[1])}</ci>"
        if k == "fn":
            inner = "".join(emit(a) for a in node[2])
            return f"<apply><ci>{_sx.escape(node[1])}</ci>{inner}</apply>"
        _, op, args = node
        return f"<apply><{op}/>" + "".join(emit(a) for a in args) + "</apply>"

    return emit(add())


def combined_effect(base_param: str, regulators: list) -> str:
    """The library motif: a process's rate = its baseline constant times, for each regulator,
    a fold-change Hill term ``(1 + (Max-1) * X/(K+X))``. ``regulators`` is
    [{species, max_param, k_param}]. Reproduces the model's 'combined regulators multiply the
    baseline rate' convention - the standard assembly, not this model's exact hand-tuning."""
    terms = [base_param]
    for r in regulators:
        x, mx, k = r["species"], r["max_param"], r["k_param"]
        terms.append(f"(1 + ({mx} - 1) * {x} / ({k} + {x}))")
    return " * ".join(terms)


def build_subsystem(cell: str, base_param: str, apop_param: str, regulators: list,
                    values: dict, clamp: dict) -> dict:
    """Assemble a self-contained runnable subsystem: one cell whose proliferation is modulated
    by its regulators (combined_effect motif) and balanced by first-order apoptosis, with the
    regulator cytokines clamped at ``clamp`` levels (boundary species). ``regulators`` is
    [{species, max_param, k_param}]; ``values`` maps every parameter name to its known value
    (the 'assume parameters known' path). Returns a spec for ``to_sbml``."""
    species = [{"name": cell, "initial": values.get(cell + "_init", 1e6)}]
    for r in regulators:
        species.append({"name": r["species"], "initial": clamp.get(r["species"], 0.0),
                        "boundary": True})
    params = [{"name": k, "value": v} for k, v in values.items() if not k.endswith("_init")]
    rate = combined_effect(base_param, regulators) + f" * {cell}"
    return {"name": cell + "_subsystem", "species": species, "parameters": params,
            "reactions": [
                {"id": cell + "_prolif", "reactants": [], "products": [cell], "rate": rate},
                {"id": cell + "_apop", "reactants": [cell], "products": [],
                 "rate": f"{apop_param} * {cell}"}],
            "rules": []}


def to_sbml(spec: dict) -> str:
    """Emit SBML L2v4 from a model spec: {species:[{name,initial,boundary?}],
    parameters:[{name,value}], reactions:[{id,reactants,products,rate}], rules:[{target,expr}]}.
    Rate laws and rules carry MathML. One default compartment. Round-trips through sbml_import."""
    def _species(s):
        bc = ' boundaryCondition="true"' if s.get("boundary") else ""
        return (f'<species id="{s["name"]}" compartment="c" '
                f'initialAmount="{s.get("initial", 0)}"{bc} hasOnlySubstanceUnits="true"/>')
    sp = "".join(_species(s) for s in spec.get("species", []))
    pr = "".join(f'<parameter id="{p["name"]}" value="{p["value"]}"/>'
                 for p in spec.get("parameters", []))
    rx = ""
    for r in spec.get("reactions", []):
        reac = "".join(f'<speciesReference species="{x}"/>' for x in r.get("reactants", []))
        prod = "".join(f'<speciesReference species="{x}"/>' for x in r.get("products", []))
        law = f'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML">' \
              f'{infix_to_mathml(r["rate"])}</math></kineticLaw>'
        rx += (f'<reaction id="{r["id"]}" reversible="false">'
               f'{f"<listOfReactants>{reac}</listOfReactants>" if reac else ""}'
               f'{f"<listOfProducts>{prod}</listOfProducts>" if prod else ""}'
               f'{law}</reaction>')
    ru = "".join(
        f'<assignmentRule variable="{r["target"]}">'
        f'<math xmlns="http://www.w3.org/1998/Math/MathML">{infix_to_mathml(r["expr"])}</math>'
        f'</assignmentRule>' for r in spec.get("rules", []))
    return ('<?xml version="1.0" encoding="UTF-8"?>'
            '<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">'
            f'<model id="{spec.get("name","assembled")}">'
            '<listOfCompartments><compartment id="c" size="1"/></listOfCompartments>'
            f'{f"<listOfSpecies>{sp}</listOfSpecies>" if sp else ""}'
            f'{f"<listOfParameters>{pr}</listOfParameters>" if pr else ""}'
            f'{f"<listOfRules>{ru}</listOfRules>" if ru else ""}'
            f'{f"<listOfReactions>{rx}</listOfReactions>" if rx else ""}'
            '</model></sbml>')

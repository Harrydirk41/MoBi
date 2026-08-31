"""A layered, provenance-tagged model specification - the single data structure the whole build
pipeline fills in, one layer per modelling decision.

The four modelling targets are orthogonal fields, and every element carries a ``source`` tag so the
honest-build distinctions (given vs agent-produced vs looked-up vs fitted vs data-needed) are
first-class data, not implicit:

    frame      (tier 1-2: human gives objective+acceptance, approves scope+scale)
    nodes      (1) which state variables exist            - agent enumerates / given
    edges      (2) topology: which regulates which, signed - agent proposes
    forms      (3) rate-law motif per process             - agent proposes
    constants  (4) values, each tagged literature|fit|prior + identifiable + needs_data

This module is pure: the AGENT decisions (chosen regulators, motif) are passed IN (produced by
llm_tasks/model_assembly via the ``call`` boundary), so the spec builder is deterministic and
testable. ``to_subsystem`` bridges a spec to the existing SBML assembly.
"""

from __future__ import annotations

from typing import Any

from pkpd_agent.engines import model_assembly as MA


# ─────────────────────────── provenance vocabulary ───────────────────────────
GIVEN, AGENT, LITERATURE, FIT, PRIOR = "given", "agent", "literature", "fit", "prior"


def _sign(direction: str | None, max_value: float | None) -> str:
    if max_value is not None:
        return "up" if max_value > 1 else "down"
    return direction or "up"


def build_spec(frame: dict, node: str, regulators: list, motif: dict,
               prov: dict, levels: dict, truth_regulators: "set | None" = None) -> dict:
    """Assemble the layered spec for one node from the agent's decisions + config data.

    ``frame``       : {objective, acceptance, scope, scale} (human tier; may be partial).
    ``regulators``  : the agent's chosen [{cytokine, direction, basis}] (layer 2, topology).
    ``motif``       : the agent's rate-law choice (layer 3).
    ``prov``/``levels`` : config data (param provenance, steady-state levels).
    ``truth_regulators`` : the model's own regulator set for this node, used ONLY to tag each edge
                           with a verify flag (never fed to the agent).
    """
    truth = set(truth_regulators or [])
    regs = [r for r in regulators if r.get("cytokine") in levels and r["cytokine"] != node]

    # ── layer 1: nodes (the built node + its regulators as boundary inputs) ──
    nodes = [{"id": node, "kind": "cytokine", "unit": "ng/mL", "role": "dynamic", "source": AGENT}]
    for r in regs:
        nodes.append({"id": r["cytokine"], "kind": "cytokine", "role": "input",
                      "source": GIVEN})

    # ── layer 2: edges (topology) ──
    edges, constants = [], []
    excess = []
    for r in regs:
        c = r["cytokine"]
        mx, from_lit = MA_max_provenance(prov, node, c)
        edges.append({"src": c, "dst": node, "process": "secretion",
                      "sign": _sign(r.get("direction"), mx), "source": AGENT,
                      "basis": r.get("basis"),
                      "verify": {"in_model": c in truth}})
        # ── layer 4: constants for this edge ──
        if mx is not None:
            constants.append({"param": f"M_{c}", "value": float(mx), "provenance": LITERATURE,
                              "ref": _ref(prov, node, c), "identifiable": None,
                              "needs_data": None, "load_bearing": None})
            mval = float(mx)
        else:
            mval = 1.5 if _sign(r.get("direction"), None) == "up" else 0.6
            constants.append({"param": f"M_{c}", "value": mval, "provenance": PRIOR,
                              "ref": None, "identifiable": None, "needs_data": "in_vitro_fold",
                              "load_bearing": None})
        constants.append({"param": f"K_{c}", "value": float(levels[c]), "provenance": FIT,
                          "ref": None, "identifiable": False, "needs_data": "dose_response",
                          "load_bearing": "pending"})
        excess.append((mval - 1.0) * 0.5)

    # ── layer 3: forms (rate-law motif) + a purely-computable verification ──
    responds = True
    if motif.get("combination") == "capped_sum" and motif.get("cap") not in (None, "", 0):
        responds = sum(excess) <= float(motif["cap"])       # cap binding at baseline -> unresponsive
    forms = [{"process": "secretion", "order": motif.get("proliferation_order"),
              "combination": motif.get("combination"), "cap": motif.get("cap"),
              "per_regulator": motif.get("per_regulator", "hill"),
              "reason": motif.get("reason"), "source": AGENT,
              "verify": {"steady_state_stable": motif.get("proliferation_order") == "zeroth",
                         "responds_to_single_knockdown": responds}}]

    # ── baseline rate (fitted, identifiable) + clearance (literature) ──
    clr = _clearance_param(prov, node)
    if clr:
        constants.append({"param": clr, "value": float(prov[clr]["value_from_reference"]),
                          "provenance": LITERATURE, "ref": _ref_name(prov, clr),
                          "identifiable": None, "needs_data": None, "load_bearing": None})
    constants.append({"param": f"kg_{node}", "value": None, "provenance": FIT, "ref": None,
                      "identifiable": True, "needs_data": None, "load_bearing": None})

    return {"frame": frame, "target": node, "nodes": nodes, "edges": edges,
            "forms": forms, "constants": constants}


def provenance_rollup(spec: dict) -> dict:
    """Count constants by provenance and flag the 'must-fit' / 'needs-data' subset - the honest
    four-quadrant view falls straight out of the tags."""
    roll = {LITERATURE: 0, FIT: 0, PRIOR: 0}
    needs_data = []
    for c in spec["constants"]:
        roll[c["provenance"]] = roll.get(c["provenance"], 0) + 1
        if c.get("needs_data"):
            needs_data.append((c["param"], c["needs_data"]))
    edges = spec["edges"]
    verified = sum(1 for e in edges if e["verify"].get("in_model"))
    return {"constants_by_source": roll, "needs_data": needs_data,
            "edges_total": len(edges), "edges_matching_model": verified}


def to_subsystem(spec: dict, prov: dict, levels: dict, clamp: "dict | None" = None) -> dict:
    """Bridge the spec to the existing SBML assembly: materialise kg (fit to the node's level under
    the chosen motif) and hand a build_subsystem spec ready for to_sbml."""
    node = spec["target"]
    reg_specs, values, excess = [], {}, []
    for e in spec["edges"]:
        c = e["src"]
        mx = next(x["value"] for x in spec["constants"] if x["param"] == f"M_{c}")
        reg_specs.append({"species": c, "max_param": f"M_{c}", "k_param": f"K_{c}"})
        values[f"M_{c}"] = mx; values[f"K_{c}"] = levels[c]
        excess.append((mx - 1.0) * 0.5)
    motif = {"proliferation_order": spec["forms"][0]["order"],
             "combination": spec["forms"][0]["combination"], "cap": spec["forms"][0]["cap"]}
    if motif["combination"] == "capped_sum":
        s = sum(excess); cap = motif["cap"]
        eff = min(cap, 1.0 + s) if cap not in (None, "", 0) else 1.0 + s
    else:
        eff = 1.0
        for e in excess:
            eff *= 1.0 + e
    clr = _clearance_param(prov, node)
    kcl = float(prov[clr]["value_from_reference"])
    kg = float(levels[node]) * kcl / eff
    values[f"kg_{node}"] = kg; values[clr] = kcl; values[f"{node}_init"] = 0.0
    if clamp is None:
        clamp = {e["src"]: levels[e["src"]] for e in spec["edges"]}
    return MA.build_subsystem(node, f"kg_{node}", clr, reg_specs, values, clamp=clamp, motif=motif)


# ─────────────────────────── small provenance helpers ───────────────────────────
def MA_max_provenance(prov, node, cyt):
    import re
    for n, p in prov.items():
        m = re.search(r"(?i)([A-Za-z0-9]+)Sec[A-Za-z0-9]*_Maxby([A-Za-z0-9]+)", n)
        if m and m.group(1) == node and m.group(2) == cyt:
            return p.get("value_from_reference"), p.get("from_literature")
    return None, None


def _ref(prov, node, cyt):
    import re
    for n, p in prov.items():
        m = re.search(r"(?i)([A-Za-z0-9]+)Sec[A-Za-z0-9]*_Maxby([A-Za-z0-9]+)", n)
        if m and m.group(1) == node and m.group(2) == cyt:
            return p.get("reference") or p.get("citation")
    return None


def _ref_name(prov, name):
    p = prov.get(name) or {}
    return p.get("reference") or p.get("citation")


def _clearance_param(prov, node):
    import re
    return next((n for n in prov if re.match(rf"(?i)(kcl|kd)_{node}\b", n)), None)

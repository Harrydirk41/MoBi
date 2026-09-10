"""Sealed sandbox for a continuously-running, leak-safe modeling agent (e.g. Claude Code).

The internal decision loop is leak-safe only because the agent sees tool RETURNS. A
filesystem-capable agent (Claude Code with Bash/Read) is not: it can `cat` the sibling
reference model ``json/<Compound>-Model.json``, the evaluation report, the ``answer_key/``,
and the held-out concentrations - i.e. the answers. Letting such an agent "work
continuously" therefore REQUIRES a physically sealed workspace.

``seal_case`` splits a case into two trees:

  * ``workspace/`` - answer-free. Contains ONLY what a modeler legitimately has: the
    blanked snapshot, the BUILDING-split clinical data, the physchem/biology inputs, and
    the task text. The agent runs here (cwd), so ``ls``/``grep``/``cat`` find nothing that
    reveals the target.
  * ``judge/``     - the answers: the reference snapshot, the held-out concentrations, the
    building/verification split, and the set of forbidden tokens (reference fitted numbers,
    held-out dataset names/concentrations). NEVER given to the agent; used only to grade and
    to drive the no-cheat verifier.

``verify_workspace`` is the seal check (planted-leak safe): the blanked snapshot must carry
NO fitted (ParameterIdentification-origin) value, no held-out dataset name/concentration may
appear in any workspace file, and no reference fitted number may appear (to significant
figures). It also scans an optional agent transcript for the same tokens. Run it after
sealing and after a session - a non-empty report means the wall has a hole.

Nothing here runs PK-Sim; it is pure stdlib file surgery, so it is fully testable offline.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from typing import Any, Iterator

from . import osp_split

# A value whose provenance is a parameter identification is a FITTED answer - it must never
# be visible to the agent. This is the same origin tag the dissolution leak carried.
_FITTED_SOURCE = "ParameterIdentification"


# --------------------------------------------------------------------------- #
# Walking a snapshot for parameters + their provenance
# --------------------------------------------------------------------------- #

def _iter_params(node: Any, path: str = "") -> Iterator[tuple[str, dict]]:
    """Yield (dotted-path, param-dict) for every dict in the snapshot that carries a numeric
    ``Value`` (these are the model's parameters, wherever they sit: compound, process,
    formulation, individual...). Recursive and structure-agnostic so nothing is missed."""
    if isinstance(node, dict):
        if isinstance(node.get("Value"), (int, float)):
            name = node.get("Name") or "?"
            yield (f"{path}/{name}", node)
        for k, v in node.items():
            yield from _iter_params(v, f"{path}/{k}")
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from _iter_params(v, f"{path}[{i}]")


def fitted_leaks(snap: dict) -> list[dict]:
    """Every parameter in ``snap`` whose value ORIGIN is a parameter identification, i.e. a
    fitted answer that a blanked snapshot must NOT still carry. An empty list means the
    snapshot is answer-free (for the fitted layer). This generalises the dissolution-leak
    check to the whole snapshot, so no future blanking miss slips through."""
    out = []
    for pth, p in _iter_params(snap):
        src = ((p.get("ValueOrigin") or {}).get("Source"))
        if src == _FITTED_SOURCE:
            out.append({"path": pth, "name": p.get("Name"),
                        "value": p.get("Value"), "unit": p.get("Unit"),
                        "origin": p.get("ValueOrigin")})
    return out


def _neutral_for(param: dict) -> float:
    """A benign naive value to replace a redacted fitted number with, so the sealed model
    still runs but carries no answer. Reference concentrations (µmol/l) fall back to the
    PK-Sim default of 1.0; everything else to 1.0 as well (the redaction is flagged in the
    manifest, so any non-obvious case is visible)."""
    return 1.0


def redact_fitted(snap: dict) -> tuple[dict, list[dict]]:
    """Return (snap, redactions) with EVERY fitted (ParameterIdentification) value removed, so
    the snapshot is safe to hand a filesystem-capable agent. Two shapes:

      * ``Simulations[].Parameters`` are path-keyed VALUE OVERRIDES baked into the built
        simulation (e.g. fitted gut-wall permeability in each gut Neighborhood). The PI-origin
        ones are DROPPED, so the built simulation uses the building-block value the agent
        controls instead of the fitted override.
      * a fitted value anywhere else (e.g. a fitted enzyme reference concentration in an
        ExpressionProfile) is NEUTRALISED to a naive default and its origin re-tagged, so the
        agent must determine it rather than read it.

    Operates on the passed dict in place and also returns it. Every change is recorded.
    """
    redactions: list[dict] = []

    for si, sim in enumerate(snap.get("Simulations") or []):
        ps = sim.get("Parameters")
        if isinstance(ps, list):
            kept = []
            for p in ps:
                if isinstance(p, dict) and \
                        ((p.get("ValueOrigin") or {}).get("Source")) == _FITTED_SOURCE:
                    redactions.append({"action": "dropped-override",
                                       "where": f"Simulations[{si}]/{sim.get('Name')}",
                                       "path": p.get("Path"), "value": p.get("Value")})
                else:
                    kept.append(p)
            sim["Parameters"] = kept

    def walk(node: Any, path: str) -> None:
        if isinstance(node, dict):
            if isinstance(node.get("Value"), (int, float)) and \
                    ((node.get("ValueOrigin") or {}).get("Source")) == _FITTED_SOURCE:
                old = node["Value"]
                node["Value"] = _neutral_for(node)
                node["ValueOrigin"] = {"Source": "Unknown",
                                       "Description": "sandbox redaction (fitted value removed)"}
                redactions.append({"action": "neutralized", "where": path,
                                   "unit": node.get("Unit"), "old": old, "new": node["Value"]})
            for k, v in node.items():
                if k == "Simulations":
                    continue                                # handled above
                walk(v, f"{path}/{k}")
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{path}[{i}]")

    for k, v in snap.items():
        if k == "Simulations":
            continue
        walk(v, f"/{k}")
    return snap, redactions


def strip_observed_blocks(snap: dict) -> dict:
    """Remove the snapshot blocks that carry observed CONCENTRATIONS and the saved parameter-
    identification SETUP - neither is load-bearing for the agent, and both leak answers:

      * ``ObservedData`` holds every observed curve, building AND held-out; the harness scores
        against the building data passed to the tools, not this block, and the building curves
        are in ``task.input.json`` anyway - so the whole block goes (name-agnostic, so a
        compact snapshot id that doesn't match the split's dataset name cannot slip through).
      * ``ParameterIdentifications`` is the saved PI setup, whose stored targets are fitted
        answers (e.g. Sildenafil's 1.37405); PK-Sim's own snap step already blanks it when
        pruning, so removing it changes nothing about the run.
      * each simulation's ``OutputMappings`` / ``ObservedData`` reference those datasets by
        name (goodness-of-fit overlay only, not needed to simulate).

    Operates in place; returns a count of what was removed."""
    removed = {"ObservedData": 0, "ParameterIdentifications": 0,
               "OutputMappings": 0, "sim_ObservedData": 0}
    for block in ("ObservedData", "ParameterIdentifications"):
        if isinstance(snap.get(block), list):
            removed[block] = len(snap[block])
            snap[block] = []
    for sim in snap.get("Simulations") or []:
        if isinstance(sim.get("OutputMappings"), list):
            removed["OutputMappings"] += len(sim["OutputMappings"])
            sim["OutputMappings"] = []
        if isinstance(sim.get("ObservedData"), list):
            removed["sim_ObservedData"] += len(sim["ObservedData"])
            sim["ObservedData"] = []
    return removed


def reference_fitted_values(ref_snap: dict) -> list[float]:
    """The numeric values the reference model FITTED (the answers), for the forbidden-token
    set the verifier scans the workspace for."""
    vals = []
    for _, p in _iter_params(ref_snap):
        if ((p.get("ValueOrigin") or {}).get("Source")) == _FITTED_SOURCE:
            v = p.get("Value")
            if isinstance(v, (int, float)):
                vals.append(float(v))
    return vals


# --------------------------------------------------------------------------- #
# Numeric tokens (for the leak scan)
# --------------------------------------------------------------------------- #

def _num_tokens(v: float) -> set[str]:
    """Distinctive string forms of a fitted number, at 4-6 significant figures. Only values
    with >=4 significant digits produce tokens, so round quantities (1, 2, 30, 0.5) - which
    are NOT distinctive of a fit and would false-positive on placeholders/constants - are
    ignored. A genuinely fitted value like 17.5448 yields '17.54'/'17.545'/'17.5448'."""
    if not isinstance(v, (int, float)) or v == 0:
        return set()
    toks: set[str] = set()
    for sig in (4, 5, 6):
        s = f"{v:.{sig}g}"
        digits = sum(c.isdigit() for c in s.split("e")[0].lstrip("-").replace(".", "").lstrip("0"))
        if digits >= 4:
            toks.add(s)
    return toks


def _scan_text(text: str, tokens: set[str]) -> list[str]:
    """Return the tokens that appear as a contiguous numeric run in ``text``."""
    hits = []
    for t in tokens:
        # a token is a number; require it not be immediately flanked by more digits, so
        # '17.545' does not spuriously match inside '917.5455'
        pat = r"(?<![\d.])" + re.escape(t) + r"(?![\d])"
        if re.search(pat, text):
            hits.append(t)
    return hits


# --------------------------------------------------------------------------- #
# Sealing a case
# --------------------------------------------------------------------------- #

def _workspace_input(inp: dict, build_datasets: set[str], has_split: bool) -> dict:
    """The input.json the agent gets: identical to the source EXCEPT clinical_observed_data
    is reduced to the BUILDING datasets (held-out concentrations removed). Study DESIGNS are
    kept - a modeler predicting held-out studies legitimately knows their dose/route."""
    out = json.loads(json.dumps(inp))                       # deep copy
    gd = out.get("given_data") or {}
    obs = gd.get("clinical_observed_data") or []
    if has_split:
        gd["clinical_observed_data"] = [o for o in obs if o.get("dataset") in build_datasets]
    out["given_data"] = gd
    return out


def seal_case(blanked_path: str, input_path: str, out_dir: str,
              reference_path: str | None = None, task_name: str | None = None) -> dict:
    """Build ``out_dir/workspace`` (answer-free) and ``out_dir/judge`` (answers) from a case.

    Returns a manifest dict (also written to ``judge/manifest.json``). The workspace is what a
    continuous agent is pointed at; the judge tree is for grading + verification only.
    """
    with open(blanked_path, encoding="utf-8") as fh:
        blanked = json.load(fh)
    with open(input_path, encoding="utf-8") as fh:
        inp = json.load(fh)

    observed = ((inp.get("given_data") or {}).get("clinical_observed_data")) or []
    split = osp_split.split_studies(blanked, observed)
    has_split = bool(split.get("verification"))
    build_set = set(split.get("building") or [])
    verif_set = set(split.get("verification") or [])
    build_obs = [o for o in observed if o.get("dataset") in build_set] if has_split else observed
    held_obs = [o for o in observed if o.get("dataset") in verif_set]

    ws = os.path.join(out_dir, "workspace")
    jd = os.path.join(out_dir, "judge")
    os.makedirs(ws, exist_ok=True)
    os.makedirs(jd, exist_ok=True)

    # ---- workspace/ : answer-free -------------------------------------- #
    # The blanked snapshot still bakes fitted values into its Simulations block (and, in a
    # few models, ExpressionProfile reference concentrations) - readable by a filesystem
    # agent. Redact them from the sealed copy so the workspace carries no answer.
    redacted, redactions = redact_fitted(json.loads(json.dumps(blanked)))
    stripped = strip_observed_blocks(redacted)
    with open(os.path.join(ws, "model.blanked.json"), "w", encoding="utf-8") as fh:
        json.dump(redacted, fh, ensure_ascii=False, indent=1)
    ws_input = _workspace_input(inp, build_set, has_split)
    with open(os.path.join(ws, "task.input.json"), "w", encoding="utf-8") as fh:
        json.dump(ws_input, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(ws, "TASK.md"), "w", encoding="utf-8") as fh:
        fh.write(_task_md(inp, task_name or _case_name(input_path), split))

    # ---- judge/ : the answers ------------------------------------------ #
    if reference_path and os.path.exists(reference_path):
        shutil.copyfile(reference_path, os.path.join(jd, "reference.json"))
    with open(os.path.join(jd, "heldout.observed.json"), "w", encoding="utf-8") as fh:
        json.dump(held_obs, fh, ensure_ascii=False, indent=1)
    with open(os.path.join(jd, "split.json"), "w", encoding="utf-8") as fh:
        json.dump(split, fh, ensure_ascii=False, indent=1)

    ref_vals: list[float] = []
    if reference_path and os.path.exists(reference_path):
        with open(reference_path, encoding="utf-8") as fh:
            ref_vals = reference_fitted_values(json.load(fh))
    # held-out concentration TOKENS that don't also occur in building data - only these are
    # distinctive of the held-out curves. Exclude at the token (significant-figure) level, not
    # exact float, so a building value that merely rounds to the same token as a held-out one
    # does not false-positive on the building concentrations the agent legitimately has.
    build_tokens: set[str] = set()
    for v in _heldout_conc_values(build_obs):
        build_tokens |= _num_tokens(v)
    held_tokens: set[str] = set()
    for v in _heldout_conc_values(held_obs):
        held_tokens |= _num_tokens(v)
    forbidden = {
        "heldout_datasets": sorted(verif_set),
        "reference_fitted_values": ref_vals,
        "heldout_conc_tokens": sorted(held_tokens - build_tokens),
    }
    with open(os.path.join(jd, "forbidden.json"), "w", encoding="utf-8") as fh:
        json.dump(forbidden, fh, ensure_ascii=False, indent=1)

    manifest = {
        "task": task_name or _case_name(input_path),
        "source": {"blanked": os.path.abspath(blanked_path),
                   "input": os.path.abspath(input_path),
                   "reference": os.path.abspath(reference_path) if reference_path else None},
        "split": {"method": split.get("method"), "n_studies": split.get("n_studies"),
                  "n_building": len(build_obs), "n_verification": len(held_obs),
                  "held_out_studies": split.get("held_out_studies")},
        "workspace": os.path.abspath(ws),
        "judge": os.path.abspath(jd),
        "source_blank_leaks": len(fitted_leaks(blanked)),   # fitted values still in the SOURCE blanked snapshot
        "redactions": len(redactions),                      # fitted values removed for the sealed copy
        "stripped_blocks": stripped,                        # observed/PI blocks removed from the sealed copy
        "sealed_leaks": len(fitted_leaks(redacted)),        # must be 0 - the sealed copy is answer-free
    }
    with open(os.path.join(jd, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=1)
    return manifest


def _heldout_conc_values(held_obs: list[dict]) -> list[float]:
    """Distinctive concentration values from the held-out data - if any appears in a
    workspace file, a held-out curve leaked."""
    out: list[float] = []
    for o in held_obs:
        for v in (o.get("conc_mg_L") or []):
            if isinstance(v, (int, float)):
                out.append(float(v))
    return out


def _case_name(input_path: str) -> str:
    base = os.path.basename(input_path)
    for suf in (".input.json", "-Model.input.json", ".json"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base


def _task_md(inp: dict, name: str, split: dict) -> str:
    obj = inp.get("objective") or "Build the whole-body PBPK model for this compound."
    bg = (inp.get("background") or {}).get("description") or ""
    n_b = len(split.get("building") or [])
    n_v = len(split.get("verification") or [])
    return (
        f"# PBPK modeling task: {name}\n\n"
        f"{obj}\n\n"
        f"## Background\n{bg}\n\n"
        "## What is in this workspace\n"
        "- `model.blanked.json` - the whole-body PBPK snapshot. The physiology is fixed; the\n"
        "  drug parameters marked as placeholders are unknown and must be determined.\n"
        "- `task.input.json` - the objective, known biology, literature physicochemical\n"
        "  priors, and the BUILDING clinical datasets you fit to.\n"
        f"- Building datasets: {n_b}. Held-out (verification) studies exist ({n_v} datasets)\n"
        "  and are NOT in this workspace - you are graded on how well the model you build\n"
        "  PREDICTS them. Do not try to find them; they are intentionally absent.\n\n"
        "## The tools\n"
        "Drive the model with `osp_agent_cli.py` (run from this directory):\n"
        "```\n"
        "python -m examples.osp_agent_cli inspect      # objective, biology, priors, data, model\n"
        "python -m examples.osp_agent_cli options      # editable params, legal methods, addable mechanisms\n"
        "python -m examples.osp_agent_cli sweep   --estimate '{\"Lipophilicity\":[0,5]}'   # choose distribution method\n"
        "python -m examples.osp_agent_cli optimize --estimate '{\"<param>\":[lo,hi]}' --structure '{...}'\n"
        "```\n"
        "You decide the modeling JUDGMENT (structure, which parameters to fit vs fix); the\n"
        "optimizer does the fitting. Iterate until the fit is good and the parameters are\n"
        "identifiable, then stop. Grading against the held-out data happens afterwards.\n"
    )


# --------------------------------------------------------------------------- #
# Verifying the seal
# --------------------------------------------------------------------------- #

def _workspace_files(workspace_dir: str) -> list[str]:
    out = []
    for root, _dirs, files in os.walk(workspace_dir):
        for f in files:
            if f == ".osp_session.json":            # the agent's own scratch state - not a leak source
                continue
            out.append(os.path.join(root, f))
    return out


def verify_workspace(workspace_dir: str, judge_dir: str | None = None,
                     transcript_path: str | None = None) -> dict:
    """Assert the workspace (and an optional agent transcript) leaks no answer. Returns
    ``{"ok": bool, "leaks": [...], "info": [...], "checks": {...}}``. A block-level leak -> ok
    False; ``info`` items are reported but do not fail the seal.

    Block-level (the answers):
      1. the blanked snapshot carries NO fitted (ParameterIdentification) value;
      2. no reference FITTED number appears (to significant figures);
      3. no held-out CONCENTRATION value appears.
    Info-level (not the answer, reported for transparency):
      * a held-out dataset NAME appears - a study's identity/design (author, route, dose,
        cohort) is legitimately knowable by a modeler predicting it; only its concentration
        curve is secret, and that is checked at block level.
      (2-3 and the info scan need judge_dir; 1 is structural and always runs.)
    """
    leaks: list[dict] = []
    info: list[dict] = []
    checks: dict[str, Any] = {}

    # (1) structural: blanked snapshot must be fitted-value-free.
    ws_snap = os.path.join(workspace_dir, "model.blanked.json")
    if os.path.exists(ws_snap):
        with open(ws_snap, encoding="utf-8") as fh:
            snap = json.load(fh)
        fl = fitted_leaks(snap)
        checks["blanked_snapshot_fitted_values"] = len(fl)
        for x in fl:
            leaks.append({"kind": "fitted-value-in-snapshot", "where": "model.blanked.json",
                          "detail": x})
    else:
        checks["blanked_snapshot_present"] = False
        leaks.append({"kind": "missing", "where": "model.blanked.json",
                      "detail": "workspace has no model.blanked.json"})

    # gather workspace text once for the token scans.
    ws_texts: dict[str, str] = {}
    for f in _workspace_files(workspace_dir):
        try:
            with open(f, encoding="utf-8") as fh:
                ws_texts[os.path.relpath(f, workspace_dir)] = fh.read()
        except (UnicodeDecodeError, OSError):
            continue
    if transcript_path and os.path.exists(transcript_path):
        try:
            with open(transcript_path, encoding="utf-8") as fh:
                ws_texts[f"<transcript:{os.path.basename(transcript_path)}>"] = fh.read()
        except (UnicodeDecodeError, OSError):
            pass

    if judge_dir and os.path.exists(os.path.join(judge_dir, "forbidden.json")):
        with open(os.path.join(judge_dir, "forbidden.json"), encoding="utf-8") as fh:
            forb = json.load(fh)

        # held-out dataset names (INFO - study identity/design is legitimately knowable)
        held_names = [d for d in (forb.get("heldout_datasets") or []) if d]
        seen_names: set[str] = set()
        for rel, text in ws_texts.items():
            for name in held_names:
                if name and name in text and name not in seen_names:
                    seen_names.add(name)
                    info.append({"kind": "heldout-dataset-name", "where": rel, "detail": name})
        checks["heldout_dataset_names_scanned"] = len(held_names)

        # (2) reference fitted numbers (BLOCK)
        ref_tokens: set[str] = set()
        for v in (forb.get("reference_fitted_values") or []):
            ref_tokens |= _num_tokens(v)
        for rel, text in ws_texts.items():
            for hit in _scan_text(text, ref_tokens):
                leaks.append({"kind": "reference-fitted-value", "where": rel, "detail": hit})
        checks["reference_value_tokens_scanned"] = len(ref_tokens)

        # (3) held-out concentrations (BLOCK) - tokens precomputed to exclude any shared
        # with building data, so a building concentration cannot false-positive.
        conc_tokens = set(forb.get("heldout_conc_tokens") or [])
        for rel, text in ws_texts.items():
            for hit in _scan_text(text, conc_tokens):
                leaks.append({"kind": "heldout-concentration", "where": rel, "detail": hit})
        checks["heldout_conc_tokens_scanned"] = len(conc_tokens)
    else:
        checks["judge_dir"] = "absent - only structural check (1) ran"

    return {"ok": not leaks, "leaks": leaks, "info": info, "checks": checks}

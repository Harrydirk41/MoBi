"""Tests for the sealed sandbox + no-cheat verifier (pkpd_agent.engines.osp_sandbox).

The sandbox is the wall that lets a filesystem-capable agent (Claude Code) run continuously
without reaching the reference model or held-out data. These tests pin the invariants that
make the wall trustworthy: fitted values and held-out concentrations never survive into the
sealed workspace, and the verifier catches them if they do.
"""

from __future__ import annotations

import glob
import json
import os

import pytest

from pkpd_agent.engines import osp_sandbox as S

_LIB = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                    "OSP-PBPK-Model-Library"))


def _pi(v, unit=None):
    return {"Value": v, "Unit": unit, "ValueOrigin": {"Source": "ParameterIdentification"}}


def _synthetic_snapshot() -> dict:
    return {
        "Compounds": [{
            "Name": "Drug",
            "Parameters": [{"Name": "Molecular weight", "Value": 300.0}],
            "Processes": [{"Molecule": "CYP3A4", "InternalName": "Meta",
                           "Parameters": [dict(Name="kcat", **_pi(17.5448))]}],
        }],
        "Formulations": [{"Name": "Tab", "Parameters": [dict(Name="Dissolution time (50% dissolved)",
                                                             **_pi(1.7958))]}],
        "ExpressionProfiles": [{"Molecule": "CYP3A4",
                                "Parameters": [_pi(3.4842, "µmol/l")]}],
        "Simulations": [{
            "Name": "PO_1mg",
            "Parameters": [{"Path": "Neigh|P", **_pi(0.00015477)},
                           {"Path": "Infusion", "Value": 1.0}],   # non-PI, must survive
            "OutputMappings": [{"ObservedData": "StudyA", "Path": "x"}],
            "ObservedData": ["StudyA"],
        }],
        "ObservedData": [{"Name": "StudyA",
                          "Columns": [{"Values": [0.111, 0.222], "Unit": "mg/l"}]}],
        "ParameterIdentifications": [{"IdentificationParameters":
                                      [{"Parameters": [{"Value": 1.37405}]}]}],
    }


# --------------------------------------------------------------------------- #
# leak detection
# --------------------------------------------------------------------------- #

def test_fitted_leaks_detects_pi_values_everywhere():
    leaks = S.fitted_leaks(_synthetic_snapshot())
    paths = " ".join(l["path"] for l in leaks)
    # a fitted value in the compound process, the formulation, the expression profile,
    # and the simulation overrides must all be found
    assert any("Processes" in l["path"] for l in leaks)
    assert any("Formulations" in l["path"] for l in leaks)
    assert any("ExpressionProfiles" in l["path"] for l in leaks)
    assert any("Simulations" in l["path"] for l in leaks)
    assert "Compounds" in paths


def test_fitted_leaks_empty_when_clean():
    snap = {"Compounds": [{"Parameters": [{"Name": "MW", "Value": 300.0}]}]}
    assert S.fitted_leaks(snap) == []


# --------------------------------------------------------------------------- #
# redaction
# --------------------------------------------------------------------------- #

def test_redact_fitted_removes_every_fitted_value():
    snap, redactions = S.redact_fitted(_synthetic_snapshot())
    assert redactions                                   # something was redacted
    assert S.fitted_leaks(snap) == []                   # and nothing fitted remains

    # simulation PI override dropped, non-PI sim param kept
    sim_params = snap["Simulations"][0]["Parameters"]
    paths = [p.get("Path") for p in sim_params]
    assert "Neigh|P" not in paths                        # the fitted override is gone
    assert "Infusion" in paths                           # the structural param survives

    # a fitted value elsewhere (expression profile) is neutralised + retagged, not left as-is
    ep_val = snap["ExpressionProfiles"][0]["Parameters"][0]
    assert ep_val["Value"] != 3.4842
    assert ep_val["ValueOrigin"]["Source"] != "ParameterIdentification"


# --------------------------------------------------------------------------- #
# stripping observed / PI blocks
# --------------------------------------------------------------------------- #

def test_strip_observed_blocks_empties_all_observed_and_pi():
    snap = _synthetic_snapshot()
    removed = S.strip_observed_blocks(snap)
    assert snap["ObservedData"] == []
    assert snap["ParameterIdentifications"] == []
    assert snap["Simulations"][0]["OutputMappings"] == []
    assert snap["Simulations"][0]["ObservedData"] == []
    assert removed["ObservedData"] == 1 and removed["ParameterIdentifications"] == 1


# --------------------------------------------------------------------------- #
# numeric tokens
# --------------------------------------------------------------------------- #

def test_num_tokens_distinctive_only():
    assert S._num_tokens(17.5448)                        # a fitted-looking value -> tokens
    assert "17.54" in S._num_tokens(17.5448)
    assert S._num_tokens(2.0) == set()                   # a round value -> no token
    assert S._num_tokens(30.0) == set()
    assert S._num_tokens(0.0) == set()


# --------------------------------------------------------------------------- #
# seal_case
# --------------------------------------------------------------------------- #

def _write_case(tmp_path):
    blanked = _synthetic_snapshot()
    # make it a plausible blanked snapshot: compound/formulation already blanked, only the
    # simulation + expression-profile fitted values remain (the real leak surface).
    reference = _synthetic_snapshot()
    observed = {"given_data": {"clinical_observed_data": [
        {"dataset": "StudyA", "study": "A", "route": "PO", "dose": "1 mg",
         "conc_mg_L": [0.111, 0.222]}]},
        "objective": "Build it."}
    bp = tmp_path / "blanked.json"; bp.write_text(json.dumps(blanked))
    ip = tmp_path / "input.json"; ip.write_text(json.dumps(observed))
    rp = tmp_path / "reference.json"; rp.write_text(json.dumps(reference))
    return str(bp), str(ip), str(rp)


def test_seal_case_produces_answer_free_workspace(tmp_path):
    bp, ip, rp = _write_case(tmp_path)
    out = tmp_path / "sb"
    manifest = S.seal_case(bp, ip, str(out), reference_path=rp, task_name="Drug")

    assert manifest["sealed_leaks"] == 0
    ws = os.path.join(str(out), "workspace")
    jd = os.path.join(str(out), "judge")
    assert os.path.exists(os.path.join(ws, "model.blanked.json"))
    assert os.path.exists(os.path.join(ws, "task.input.json"))
    assert os.path.exists(os.path.join(jd, "reference.json"))
    assert os.path.exists(os.path.join(jd, "forbidden.json"))

    assert os.path.exists(os.path.join(ws, "PROMPT.md"))  # ready-to-paste CC task
    assert os.path.exists(os.path.join(ws, "TASK.md"))
    assert os.path.exists(os.path.join(jd, "heldout.observed.json"))

    sealed = json.load(open(os.path.join(ws, "model.blanked.json")))
    assert S.fitted_leaks(sealed) == []                  # no fitted value in the sealed copy
    assert sealed.get("ObservedData") == []              # observed block stripped
    assert sealed.get("ParameterIdentifications") == []  # PI setup stripped

    # the reference fitted values live only in judge/, never in the workspace
    ws_text = open(os.path.join(ws, "model.blanked.json")).read()
    assert "17.5448" not in ws_text
    assert "1.37405" not in ws_text


def test_strip_expression_profiles_removes_block():
    snap = {"ExpressionProfiles": [{"Molecule": "CYP3A4"}, {"Molecule": "P-gp"}],
            "Individuals": [{"ExpressionProfiles": [{"Molecule": "CYP3A4"}]}],
            "Simulations": [{"ExpressionProfiles": [{"Molecule": "CYP3A4"}]}]}
    n = S.strip_expression_profiles(snap)
    assert n == 4                                        # 2 top-level + 1 individual + 1 sim
    assert snap["ExpressionProfiles"] == []
    assert snap["Individuals"][0]["ExpressionProfiles"] == []
    assert snap["Simulations"][0]["ExpressionProfiles"] == []


def test_label_scrub_blanks_enzyme_names_in_individual_labels():
    snap = {"Individuals": [{"Name": "Japanese (P-gp modified, CYP3A4 36 h)"}],
            "Simulations": [{"Individual": "Japanese (P-gp modified, CYP3A4 36 h)"}],
            "SimulationClassifications": [{"Name": "CYP2D6_EMs_PMs"}]}
    n = S._scrub_mechanism_names(snap, ["CYP3A4", "P-gp", "CYP2D6"])
    assert n >= 3
    # the simulation's individual reference stays consistent with the renamed individual
    assert snap["Individuals"][0]["Name"] == snap["Simulations"][0]["Individual"]
    assert "CYP3A4" not in snap["Individuals"][0]["Name"]
    assert "CYP2D6" not in snap["SimulationClassifications"][0]["Name"]


def test_hard_seal_strips_enzyme_identity_but_normal_keeps_it(tmp_path):
    # A hard (structure-blind) snapshot: processes already stripped, but the
    # ExpressionProfiles still NAME the clearing enzyme - which for many compounds IS
    # the answer. Sealing a hard case must remove that block from the workspace so a
    # filesystem agent can't read the mechanism off the raw snapshot.
    hard_snap = {
        "Compounds": [{"Name": "Drug",
                       "Parameters": [{"Name": "Molecular weight", "Value": 300.0}],
                       "Processes": []}],                    # hard: structure stripped
        "ExpressionProfiles": [{"Molecule": "CYP3A4", "Parameters": []}],
        "Simulations": [{"Name": "PO_1mg", "Parameters": [], "Compounds": []}],
    }
    inp = {"given_data": {"clinical_observed_data": [
        {"dataset": "StudyA", "study": "A", "route": "PO", "dose": "1 mg",
         "conc_mg_L": [0.1, 0.05]}]},
        "background": {"candidate_clearance_molecules": [
            {"molecule": "CYP3A4"}, {"molecule": "CYP2D6"}, {"molecule": "UGT1A1"}]}}

    # hard: name must be gone from the workspace snapshot
    bp = tmp_path / "Drug-Model.hard_blanked.json"; bp.write_text(json.dumps(hard_snap))
    ip = tmp_path / "Drug-Model.hard.input.json"; ip.write_text(json.dumps(inp))
    mf = S.seal_case(str(bp), str(ip), str(tmp_path / "hard_out"), task_name="Drug")
    assert mf["hard"] is True
    assert mf["stripped_expression_profiles"] == 1
    ws_snap = tmp_path / "hard_out" / "workspace" / "model.blanked.json"
    assert json.load(open(ws_snap)).get("ExpressionProfiles") == []
    assert "CYP3A4" not in ws_snap.read_text()               # the giveaway is closed

    # the seal-time self-check confirms no mechanism name survived into the snapshot,
    # and the forbidden set records the enzyme names for the standalone verifier.
    assert mf["hard_mechanism_leaks"] == 0
    forb = json.load(open(tmp_path / "hard_out" / "judge" / "forbidden.json"))
    assert "CYP3A4" in forb["mechanism_molecules"]

    # normal (non-hard): the expressed panel legitimately stays (it is the given skeleton)
    bp2 = tmp_path / "Drug-Model.blanked.json"; bp2.write_text(json.dumps(hard_snap))
    ip2 = tmp_path / "Drug-Model.input.json"; ip2.write_text(json.dumps(inp))
    mf2 = S.seal_case(str(bp2), str(ip2), str(tmp_path / "norm_out"), task_name="Drug")
    assert mf2["hard"] is False
    assert mf2["stripped_expression_profiles"] == 0
    ws_snap2 = tmp_path / "norm_out" / "workspace" / "model.blanked.json"
    assert json.load(open(ws_snap2)).get("ExpressionProfiles")   # kept
    # normal mode has no mechanism-molecule forbidden set (the enzyme is a given)
    forb2 = json.load(open(tmp_path / "norm_out" / "judge" / "forbidden.json"))
    assert forb2["mechanism_molecules"] == []


def test_verifier_catches_enzyme_name_in_hard_snapshot(tmp_path):
    # regression guard: if a hard workspace snapshot still NAMES the answer enzyme
    # (e.g. the strip regressed), the verifier must flag it as a block-level leak -
    # but the SAME name in task.input.json's candidate pool must NOT trip it.
    ws = tmp_path / "workspace"; ws.mkdir()
    jd = tmp_path / "judge"; jd.mkdir()
    (ws / "model.blanked.json").write_text(json.dumps(
        {"Compounds": [{"Name": "Drug"}],
         "ExpressionProfiles": [{"Molecule": "CYP3A4"}]}))        # LEAK: names the enzyme
    (ws / "task.input.json").write_text(json.dumps(
        {"background": {"candidate_clearance_molecules": [{"molecule": "CYP3A4"}]}}))  # legit
    (jd / "forbidden.json").write_text(json.dumps(
        {"heldout_datasets": [], "reference_fitted_values": [],
         "heldout_conc_tokens": [], "mechanism_molecules": ["CYP3A4"]}))
    res = S.verify_workspace(str(ws), judge_dir=str(jd))
    assert res["ok"] is False
    kinds = {(l["kind"], l["where"]) for l in res["leaks"]}
    assert ("mechanism-molecule-in-snapshot", "model.blanked.json") in kinds
    # exactly one hit (the snapshot), NOT the pool in task.input.json
    assert sum(1 for l in res["leaks"] if l["kind"] == "mechanism-molecule-in-snapshot") == 1

    # clean snapshot (enzyme only in the pool) passes the mechanism check
    (ws / "model.blanked.json").write_text(json.dumps({"Compounds": [{"Name": "Drug"}]}))
    res2 = S.verify_workspace(str(ws), judge_dir=str(jd))
    assert not [l for l in res2["leaks"] if l["kind"] == "mechanism-molecule-in-snapshot"]


# --------------------------------------------------------------------------- #
# verify_workspace
# --------------------------------------------------------------------------- #

def _clean_workspace(tmp_path):
    """A minimal sealed workspace + judge with a clean snapshot and a forbidden set."""
    ws = tmp_path / "workspace"; ws.mkdir()
    jd = tmp_path / "judge"; jd.mkdir()
    (ws / "model.blanked.json").write_text(json.dumps(
        {"Compounds": [{"Parameters": [{"Name": "MW", "Value": 300.0}]}]}))
    (ws / "task.input.json").write_text(json.dumps({"given_data": {}}))
    (jd / "forbidden.json").write_text(json.dumps({
        "heldout_datasets": ["Friedman 1988 - Control"],
        "reference_fitted_values": [17.5448180963],
        "heldout_conc_tokens": ["0.001209"],
    }))
    return str(ws), str(jd)


def test_verify_workspace_passes_when_clean(tmp_path):
    ws, jd = _clean_workspace(tmp_path)
    rep = S.verify_workspace(ws, judge_dir=jd)
    assert rep["ok"] is True
    assert rep["leaks"] == []


def test_verify_workspace_flags_reference_value_leak(tmp_path):
    ws, jd = _clean_workspace(tmp_path)
    # plant the reference fitted value into a workspace file
    with open(os.path.join(ws, "notes.txt"), "w") as fh:
        fh.write("kcat came out around 17.5448 units")
    rep = S.verify_workspace(ws, judge_dir=jd)
    assert rep["ok"] is False
    assert any(l["kind"] == "reference-fitted-value" for l in rep["leaks"])


def test_verify_workspace_flags_heldout_concentration(tmp_path):
    ws, jd = _clean_workspace(tmp_path)
    with open(os.path.join(ws, "task.input.json"), "w") as fh:
        fh.write(json.dumps({"leaked": [0.001209]}))
    rep = S.verify_workspace(ws, judge_dir=jd)
    assert rep["ok"] is False
    assert any(l["kind"] == "heldout-concentration" for l in rep["leaks"])


def test_verify_workspace_heldout_name_is_info_not_leak(tmp_path):
    ws, jd = _clean_workspace(tmp_path)
    with open(os.path.join(ws, "task.input.json"), "w") as fh:
        fh.write(json.dumps({"predict": "Friedman 1988 - Control"}))
    rep = S.verify_workspace(ws, judge_dir=jd)
    assert rep["ok"] is True                             # a study name is not the answer
    assert any(i["kind"] == "heldout-dataset-name" for i in rep["info"])


def test_verify_workspace_flags_snapshot_fitted_value(tmp_path):
    ws, jd = _clean_workspace(tmp_path)
    # a fitted value that survived into the snapshot itself
    (tmp_path / "workspace" / "model.blanked.json").write_text(json.dumps(
        {"Compounds": [{"Parameters": [dict(Name="kcat", **_pi(5.0))]}]}))
    rep = S.verify_workspace(ws, judge_dir=jd)
    assert rep["ok"] is False
    assert any(l["kind"] == "fitted-value-in-snapshot" for l in rep["leaks"])


# --------------------------------------------------------------------------- #
# integration on a real library case (skipped if the library isn't present)
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not glob.glob(os.path.join(_LIB, "Triazolam", "benchmark", "*blanked*.json")),
                    reason="OSP model library not available")
def test_real_triazolam_seal_is_answer_free(tmp_path):
    blanked = glob.glob(os.path.join(_LIB, "Triazolam", "benchmark", "*Model.blanked.json"))[0]
    inp = os.path.join(_LIB, "Triazolam", "json_input", "Triazolam-Model.input.json")
    ref = glob.glob(os.path.join(_LIB, "Triazolam", "json", "*-Model.json"))
    manifest = S.seal_case(blanked, inp, str(tmp_path / "sb"),
                           reference_path=ref[0] if ref else None, task_name="Triazolam")
    assert manifest["sealed_leaks"] == 0
    assert manifest["source_blank_leaks"] > 0            # the source really did leak
    rep = S.verify_workspace(os.path.join(str(tmp_path / "sb"), "workspace"),
                             judge_dir=os.path.join(str(tmp_path / "sb"), "judge"))
    assert rep["ok"] is True, rep["leaks"]

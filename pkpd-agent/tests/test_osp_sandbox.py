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

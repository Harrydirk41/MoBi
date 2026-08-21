"""Conversational onboarding tools: registration, the build->validate->set->save loop.

Stub call_fn + dict network, and a temp out_dir so no API / MATLAB / real project folder
is touched. Exercises the whole path an onboarding agent drives.
"""

import json
import os
import tempfile
import unittest

from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.onboarding_loop_tools import register_onboarding_tools, _set_path


_NETWORK = {
    "name": "Test QSP",
    "species": [{"name": s} for s in ("ACR20", "MTX_NonResp", "TCZ_ACR20", "DAS28_CRP")],
    "parameters": [{"name": p} for p in ("F_IL6", "F_TNFa", "KD_TCZ")],
}


def _stub_call(system, user):
    # structure classify / edges / readout, tasks role classify / readout, config build
    if "role in a two-line trial" in user:
        return ('{"first_line_flags":["ACR20"],"subgroup_flag":"MTX_NonResp",'
                '"second_line_flags":["TCZ_ACR20"],"severity_states":["DAS28_CRP"]}')
    if "disease_drivers" in user and "Assign each parameter" in user:
        return '{"disease_drivers":["F_IL6"],"druggable":["F_IL6"],"calibratable":["KD_TCZ"]}'
    if '"biology"' in user or "Classify each species" in user:
        return '{"biology":["DAS28_CRP"],"drug":[],"readout":["DAS28_CRP"]}'
    if "regulatory edges" in user:
        return "[]"
    if "readout is a direct function" in user:
        return '{"drivers":["DAS28_CRP"]}'
    if "tasks.json object" in user:      # the config builder
        return ('{"vpop_target":{"mean":5,"sd":1,"band":[3,8]},'
                '"run_columns":{"first_line":{"ACR20":"ACR20"},'
                '"severity":{"baseline":"DAS28_base","readout":"DAS28_read"}}}')
    return "{}"


class _Session:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestSetPath(unittest.TestCase):
    def test_nested_set(self):
        d = {}
        _set_path(d, "vpop_target.band", [6, 20])
        self.assertEqual(d["vpop_target"]["band"], [6, 20])

    def test_through_nonobject_raises(self):
        with self.assertRaises(ValueError):
            _set_path({"a": 5}, "a.b", 1)


class TestOnboardingLoop(unittest.TestCase):
    def _reg(self, out_dir):
        reg = ToolRegistry()
        register_onboarding_tools(reg, None, {
            "network": _NETWORK, "description": "RA-like model, IL-6 driven",
            "name": "test_proj", "call": _stub_call, "out_dir": out_dir})
        return reg

    def test_registers(self):
        reg = self._reg("/tmp")
        for n in ("onboard_inspect", "onboard_build", "onboard_set", "onboard_save"):
            self.assertIn(n, reg)

    def test_inspect_reports_counts(self):
        res = self._reg("/tmp").dispatch("onboard_inspect", {}, _Session())
        self.assertEqual(res.data["n_parameters"], 3)
        self.assertTrue(res.data["have_description"])

    def test_build_validate_save_flow(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = self._reg(tmp)
            sess = _Session()
            b = reg.dispatch("onboard_build", {}, sess)
            self.assertTrue(b.ok)
            self.assertIsNotNone(sess.get("onboard_tasks"))
            # save should write files (no ERRORS expected for this stub config)
            sv = reg.dispatch("onboard_save", {}, sess)
            self.assertTrue(sv.ok, sv.message)
            self.assertTrue(os.path.isfile(os.path.join(tmp, "test_proj", "tasks.json")))
            self.assertTrue(os.path.isfile(os.path.join(tmp, "test_proj", "spec.json")))

    def test_save_blocked_on_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = self._reg(tmp)
            sess = _Session()
            reg.dispatch("onboard_build", {}, sess)
            # inject an error: a vpop_driver that is not a real parameter
            reg.dispatch("onboard_set", {"path": "vpop_drivers.F_GHOST", "value": {}}, sess)
            sv = reg.dispatch("onboard_save", {}, sess)
            self.assertFalse(sv.ok)
            self.assertFalse(os.path.isdir(os.path.join(tmp, "test_proj")))

    def test_set_fixes_and_unblocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg = self._reg(tmp)
            sess = _Session()
            reg.dispatch("onboard_build", {}, sess)
            reg.dispatch("onboard_set", {"path": "vpop_target.band", "value": [6, 20]}, sess)
            self.assertEqual(sess.get("onboard_tasks")["vpop_target"]["band"], [6, 20])

    def test_set_before_build_errors(self):
        res = self._reg("/tmp").dispatch("onboard_set", {"path": "x", "value": 1}, _Session())
        self.assertFalse(res.ok)


if __name__ == "__main__":
    unittest.main()

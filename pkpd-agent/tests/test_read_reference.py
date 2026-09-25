"""The reference library is handed to osp_inspect as a bare INDEX; to use an
analogue's structure the agent must OPEN it with osp_read_reference. This forces
real reading instead of skimming a dumped digest."""

import unittest

from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools, _reference_index
from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession

_LIB = {"mode": "all/structural", "note": "leave-one-out.",
        "models": [
            {"compound": "Sufentanil",
             "calculation_methods": ["Distribution: Schmitt"],
             "processes": [{"molecule": "CYP3A4", "type": "metabolization_first_order"}],
             "estimated_parameters": ["Lipophilicity", "Intrinsic clearance@CYP3A4"]},
            {"compound": "Midazolam",
             "calculation_methods": ["Distribution: Rodgers and Rowland"],
             "processes": [{"molecule": "CYP3A4", "type": "metabolization_first_order"}],
             "estimated_parameters": ["Lipophilicity"]},
        ]}


class TestReferenceIndex(unittest.TestCase):
    def test_index_hides_structure_keeps_names_and_molecules(self):
        idx = _reference_index(_LIB)
        self.assertEqual({m["compound"] for m in idx["available"]},
                         {"Sufentanil", "Midazolam"})
        # the index exposes involved molecules but NOT methods / parameter names
        suf = next(m for m in idx["available"] if m["compound"] == "Sufentanil")
        self.assertIn("CYP3A4", suf["involves"])
        self.assertNotIn("calculation_methods", suf)
        self.assertIn("osp_read_reference", idx["note"])

    def test_index_none_when_no_library(self):
        self.assertIsNone(_reference_index(None))


class TestReadReference(unittest.TestCase):
    def _tool(self, context_reports=None):
        reg = ToolRegistry()
        register_osp_loop_tools(reg, AgentConfig(mock=False), {
            "cli": object(), "snapshot_path": "x", "observed": [],
            "input": {"reference_library": _LIB},
            "context_reports": context_reports or {}})
        return reg.get("osp_read_reference").handler

    def test_open_analogue_returns_full_structure(self):
        r = self._tool()({"compound": "Sufentanil"}, ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        self.assertEqual(r.data["processes"][0]["molecule"], "CYP3A4")
        self.assertIn("Intrinsic clearance@CYP3A4", r.data["estimated_parameters"])
        self.assertIn("Reference model: Sufentanil", r.data["reference_markdown"])

    def test_open_analogue_serves_full_report_when_present(self):
        import tempfile, os as _os
        d = tempfile.mkdtemp()
        rp = _os.path.join(d, "Sufentanil.md")
        open(rp, "w").write("# Sufentanil model\nUses Schmitt distribution, CYP3A4.")
        r = self._tool({"Sufentanil": rp})({"compound": "Sufentanil"}, ModelingSession(goal="g"))
        self.assertIn("Schmitt distribution", r.data["report_markdown"])

    def test_unknown_analogue_lists_available(self):
        r = self._tool()({"compound": "Nope"}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)
        self.assertIn("Sufentanil", r.message)

    def test_no_library_is_denovo_error(self):
        reg = ToolRegistry()
        register_osp_loop_tools(reg, AgentConfig(mock=False), {
            "cli": object(), "snapshot_path": "x", "observed": [], "input": {}})
        r = reg.get("osp_read_reference").handler({"compound": "X"}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)
        self.assertIn("de-novo", r.message)


if __name__ == "__main__":
    unittest.main()

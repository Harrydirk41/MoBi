"""Context-data-lake tools: the agent builds its own givens from raw materials.

These verify the honest self-extraction path (read report/data -> record givens
-> the givens feed the downstream fix-vs-fit split) against the pbpk-realworld
tree, and that recording mutates the shared task input in place.
"""
import glob
import os
import unittest

from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.context_tools import register_context_tools

_RW = os.path.join(os.path.dirname(__file__), "..", "..", "pbpk-realworld")
_TRZ = os.path.join(_RW, "Triazolam", "test")


def _tools(inp):
    reg = ToolRegistry()
    register_context_tools(reg, AgentConfig(mock=False), {
        "report_path": (glob.glob(os.path.join(_TRZ, "report", "*.md")) or [None])[0],
        "data_dir": os.path.join(_TRZ, "data"),
        "input": inp})
    return {n: reg.get(n).handler for n in
            ("osp_read_report", "osp_list_studies", "osp_read_study", "osp_record_givens")}


@unittest.skipUnless(os.path.isdir(_TRZ), "pbpk-realworld not generated")
class TestContextTools(unittest.TestCase):
    def test_report_has_givens_but_not_the_answer(self):
        t = _tools({})
        r = t["osp_read_report"]({}, ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        md = r.data["report_markdown"]
        self.assertIn("CYP3A4", md)                 # enzyme identity is a given
        self.assertIn("unbound", md.lower())        # fu is a given
        self.assertNotIn("GMFE", md)                # the grade is redacted
        self.assertNotIn("Rodgers and Rowland", md)  # the chosen method is redacted

    def test_section_filter(self):
        t = _tools({})
        r = t["osp_read_report"]({"section": "in vitro"}, ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        self.assertIn("Molecular weight", r.data["report_markdown"])

    def test_list_and_read_studies(self):
        t = _tools({})
        r = t["osp_list_studies"]({}, ModelingSession(goal="g"))
        self.assertTrue(r.ok and r.data["studies"])
        s0 = r.data["studies"][0]
        self.assertIn(s0["route"], ("PO", "IV", ""))
        rd = t["osp_read_study"]({"study": s0["study"]}, ModelingSession(goal="g"))
        self.assertTrue(rd.ok and len(rd.data["points"]) >= 2)

    def test_record_givens_mutates_task_input(self):
        inp = {}
        t = _tools(inp)
        sess = ModelingSession(goal="g")
        r = t["osp_record_givens"]({"givens": [
            {"parameter": "Lipophilicity", "value": 2.42, "unit": "log", "source": "PMID x",
             "provenance": "given"},
            {"parameter": "Effective clearance", "provenance": "judged"}]}, sess)
        self.assertTrue(r.ok)
        # written into the SHARED input so options()/optimize() read the agent's givens
        lit = inp["given_data"]["literature_physicochemical"]
        self.assertEqual({e["parameter"] for e in lit},
                         {"Lipophilicity", "Effective clearance"})
        self.assertEqual(sess.get("extracted_givens"), lit)
        # a bad provenance tag is coerced, not rejected
        self.assertTrue(all(e["provenance"] in ("given", "judged") for e in lit))

    def test_record_givens_rejects_empty(self):
        t = _tools({})
        r = t["osp_record_givens"]({"givens": []}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)


class TestContextToolsNoMaterials(unittest.TestCase):
    def test_graceful_when_nothing_handed(self):
        reg = ToolRegistry()
        register_context_tools(reg, AgentConfig(mock=False),
                               {"report_path": None, "data_dir": None, "input": {}})
        r = reg.get("osp_read_report").handler({}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)
        r = reg.get("osp_list_studies").handler({}, ModelingSession(goal="g"))
        self.assertFalse(r.ok)


class TestSelfExtractNoSplitLeak(unittest.TestCase):
    """In self-extract mode osp_inspect must NOT hand the agent the reference
    model's fix-vs-fit split (a chosen method) — the agent decides it itself."""
    def _inspect(self, self_extract):
        import glob
        from pkpd_agent.tools.registry import ToolRegistry
        from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools
        snaps = glob.glob(os.path.join(
            os.path.dirname(__file__), "..", "..", "OSP-PBPK-Model-Library",
            "*", "benchmark", "*-Model.blanked.json"))
        if not snaps:
            self.skipTest("OSP library not present")
        reg = ToolRegistry()
        register_osp_loop_tools(reg, AgentConfig(mock=False), {
            "cli": object(), "snapshot_path": snaps[0], "observed": [],
            "input": {}, "self_extract": self_extract})
        return reg.get("osp_inspect").handler({}, ModelingSession(goal="g")).data

    def test_split_withheld_in_self_extract(self):
        d = self._inspect(True)
        self.assertNotIn("parameters_to_determine", d)
        self.assertNotIn("given_input_parameters", d)
        self.assertIn("fix_vs_fit_note", d)
        self.assertFalse(any("value_status" in p for p in d["current_model"]["parameters"]))

    def test_split_kept_in_normal_mode(self):
        d = self._inspect(False)
        self.assertIn("parameters_to_determine", d)
        self.assertIn("given_input_parameters", d)


if __name__ == "__main__":
    unittest.main()

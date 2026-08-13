import unittest

from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession
from pkpd_agent.tools import build_default_registry


class TestRegistry(unittest.TestCase):
    def setUp(self):
        self.cfg = AgentConfig(mock=True)
        self.reg = build_default_registry(self.cfg)
        self.session = ModelingSession(goal="test")

    def test_all_engines_contribute_tools(self):
        names = self.reg.names()
        self.assertIn("pharmpy_fit", names)
        self.assertIn("osp_simulate", names)
        self.assertIn("nca_analyze", names)

    def test_anthropic_schema_shape(self):
        schema = self.reg.to_anthropic_schema()
        self.assertTrue(schema)
        for tool in schema:
            self.assertIn("name", tool)
            self.assertIn("description", tool)
            self.assertIn("input_schema", tool)
            self.assertEqual(tool["input_schema"]["type"], "object")

    def test_dispatch_nca_computes_auc(self):
        res = self.reg.dispatch(
            "nca_analyze",
            {"times": [0, 1, 2, 3], "concentrations": [0, 10, 5, 2.5]},
            self.session,
        )
        self.assertTrue(res.ok)
        self.assertEqual(res.data["c_max"], 10)
        self.assertGreater(res.data["auc_trapezoidal"], 0)

    def test_dispatch_unknown_tool_is_error_not_exception(self):
        res = self.reg.dispatch("does_not_exist", {}, self.session)
        self.assertFalse(res.ok)
        self.assertIn("unknown tool", res.message)

    def test_handler_exception_becomes_error_result(self):
        # nca with mismatched lengths raises inside the engine -> caught
        res = self.reg.dispatch(
            "nca_analyze", {"times": [0, 1], "concentrations": [1]}, self.session
        )
        self.assertFalse(res.ok)


if __name__ == "__main__":
    unittest.main()

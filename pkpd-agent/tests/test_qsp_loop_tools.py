"""Model-agnostic loop tools exercised on the RA fixture model (no LLM)."""

import unittest

from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.qsp_topology_loop_tools import register_qsp_topology_loop_tools
from pkpd_agent.tools import qsp_loop_tools as Q
from tests.test_qsp_model import _model


class _S:
    def __init__(self):
        self._d = {}

    def get(self, k, d=None):
        return self._d.get(k, d)

    def put(self, k, v):
        self._d[k] = v


class TestGeneralLoops(unittest.TestCase):
    def setUp(self):
        self.m = _model()

    def _reg(self, fn):
        r = ToolRegistry()
        fn(r, None, {"model": self.m})
        return r

    def test_topology(self):
        r = self._reg(register_qsp_topology_loop_tools)
        s = _S()
        r.dispatch("network_propose", {"edges": [
            {"source": "TNFa", "target": "FLS", "sign": 1}]}, s)
        res = r.dispatch("network_finalize", {}, s)
        self.assertTrue(res.ok)
        self.assertGreaterEqual(res.data["topology"]["hit"], 1)

    def test_scope(self):
        r = self._reg(Q.register_qsp_scope_loop_tools)
        s = _S()
        r.dispatch("scope_propose", {"nodes": ["macrophage", "FLS", "aspirin"]}, s)
        res = r.dispatch("scope_finalize", {}, s)
        self.assertEqual(res.data["hit"], 2)          # macrophage->Macrophages, FLS

    def test_signs(self):
        r = self._reg(Q.register_qsp_sign_loop_tools)
        s = _S()
        r.dispatch("sign_predict", {"edges": [
            {"source": "TNFa", "target": "FLS", "sign": 1},
            {"source": "TGFb", "target": "endothelial", "sign": -1}]}, s)
        res = r.dispatch("sign_finalize", {}, s)
        self.assertTrue(res.ok)
        self.assertIn("accuracy", res.data)

    def test_readout(self):
        r = self._reg(Q.register_qsp_readout_loop_tools)
        s = _S()
        r.dispatch("readout_propose", {"nodes": ["macrophage", "FLS", "IL-6"]}, s)
        res = r.dispatch("readout_finalize", {}, s)
        self.assertEqual(res.data["hit"], 2)          # Macrophages, FLS are drivers; IL6 not

    def test_params(self):
        r = self._reg(Q.register_qsp_params_loop_tools)
        s = _S()
        r.dispatch("param_estimate", {"predictions": [
            {"name": "kd_FLS_Baseline", "value": 0.1}]}, s)
        res = r.dispatch("param_finalize", {}, s)
        self.assertTrue(res.ok)

    def test_sensitivity(self):
        r = self._reg(Q.register_qsp_sensitivity_loop_tools)
        s = _S()
        r.dispatch("sens_rank", {"ranked": self.m.spec.gsa_top}, s)
        res = r.dispatch("sens_finalize", {}, s)
        self.assertEqual(res.data["recall"], 1.0)
        self.assertTrue(res.data["beats_random"])


if __name__ == "__main__":
    unittest.main()

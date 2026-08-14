"""The measured-quantity principle, enforced (not merely advised).

  * measured physical CONSTANTS (MW, pKa, reference pH) can never be estimated;
  * measured-but-refinable quantities (fraction unbound, solubility) may be
    estimated only WITHIN their measured range - requested bounds are clamped;
  * freely-estimated parameters are untouched.
"""

import unittest

from pkpd_agent.engines import osp_catalog as C
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools
from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession
from pkpd_agent.engines import osp_optimize as OO


class TestTierClassification(unittest.TestCase):
    def test_constants(self):
        for n in ["Molecular weight", "Reference pH", "pKa (base)"]:
            self.assertEqual(C.param_tier(n), "constant", n)

    def test_measured_soft(self):
        for n in ["Fraction unbound (plasma, reference value)",
                  "Solubility at reference pH"]:
            self.assertEqual(C.param_tier(n), "measured_soft", n)

    def test_estimate(self):
        for n in ["Lipophilicity", "Intrinsic clearance", "Permeability"]:
            self.assertEqual(C.param_tier(n), "estimate", n)

    def test_measured_range_from_literature(self):
        lit = [{"parameter": "Fraction unbound in plasma",
                "reported_range_percent": [8.6, 12.0]}]
        self.assertEqual(
            C.measured_range("Fraction unbound (plasma, reference value)", lit),
            (0.086, 0.12))
        # a freely-estimated parameter has no measured anchor
        self.assertIsNone(C.measured_range("Lipophilicity", lit))


class TestOptimizeEnforcement(unittest.TestCase):
    def setUp(self):
        self.inp = {"given_data": {"literature_physicochemical": [
            {"parameter": "Fraction unbound in plasma",
             "reported_range_percent": [8.6, 12.0]}]}}
        self.captured = {}

        def fake_run(cli, snap, observed, estimate, **kw):
            self.captured["estimate"] = estimate
            return {"ok": True,
                    "optimized": {k: sum(v) / 2 for k, v in estimate.items()},
                    "fit": {"gmfe": 1.5}, "by_route": {}, "worst_datasets": [],
                    "params_at_bound": [], "sensitivity": {}, "n_evals": 10,
                    "fit_simulations": ["s1"]}
        self._orig = OO.run_optimization
        OO.run_optimization = fake_run
        self.reg = ToolRegistry()
        register_osp_loop_tools(self.reg, AgentConfig(mock=False), {
            "cli": object(), "snapshot_path": "x", "observed": [],
            "input": self.inp})
        self.opt = self.reg.get("osp_optimize").handler

    def tearDown(self):
        OO.run_optimization = self._orig

    def test_constant_cannot_be_estimated(self):
        r = self.opt({"estimate": {"Molecular weight": [400, 430]}},
                     ModelingSession(goal="g"))
        self.assertFalse(r.ok)
        self.assertIn("cannot be estimated", r.message)
        self.assertNotIn("estimate", self.captured)   # optimizer never ran

    def test_measured_soft_bounds_clamped_to_measured_range(self):
        r = self.opt({"estimate": {
            "Fraction unbound (plasma, reference value)": [0.001, 1.0],
            "Intrinsic clearance": [0.01, 3.0]}}, ModelingSession(goal="g"))
        self.assertTrue(r.ok)
        fu = self.captured["estimate"]["Fraction unbound (plasma, reference value)"]
        self.assertEqual(fu, [0.086, 0.12])                    # clamped
        self.assertEqual(self.captured["estimate"]["Intrinsic clearance"],
                         [0.01, 3.0])                           # untouched
        self.assertTrue(r.data.get("measured_constraints"))

    def test_in_range_bounds_are_not_widened(self):
        # a request already inside the measured range is kept as-is
        r = self.opt({"estimate": {
            "Fraction unbound (plasma, reference value)": [0.09, 0.11]}},
            ModelingSession(goal="g"))
        self.assertEqual(
            self.captured["estimate"]["Fraction unbound (plasma, reference value)"],
            [0.09, 0.11])


if __name__ == "__main__":
    unittest.main()

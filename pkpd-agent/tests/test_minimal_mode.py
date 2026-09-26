"""Minimal / purist mode: trust the LLM's trial-and-error, drop the modeling-scaffold.

These lock in WHAT minimal mode removes (the method sweep, the optimizer's directive
recommendations, staged-fit method_guidance, and the soft biological-plausibility
opinions) while keeping the tools and the code-enforced invariants.
"""

import unittest

from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools
from pkpd_agent.config import AgentConfig
from pkpd_agent.state import ModelingSession
from pkpd_agent.engines import osp_optimize as OO, osp_score


def _reg(minimal):
    reg = ToolRegistry()
    register_osp_loop_tools(reg, AgentConfig(mock=False), {
        "cli": object(), "snapshot_path": "x", "observed": [], "input": {},
        "minimal": minimal})
    return reg


class TestMinimalDropsScaffold(unittest.TestCase):
    def test_sweep_tool_absent_in_minimal(self):
        self.assertIn("osp_sweep_methods", _reg(False).names())   # full mode has it
        self.assertNotIn("osp_sweep_methods", _reg(True).names())  # minimal drops it

    def test_core_tools_still_present_in_minimal(self):
        names = _reg(True).names()
        for t in ("osp_inspect", "osp_options", "osp_try_model", "osp_optimize"):
            self.assertIn(t, names)

    def test_optimize_drops_recommendations_in_minimal(self):
        series = [{"study": "S", "observed": [[1, 0.1]], "simulated": [[1, 0.11]]}]

        def fake_run(cli, snap, observed, estimate, **kw):
            return {"ok": True, "optimized": {k: sum(v) / 2 for k, v in estimate.items()},
                    "fit": {"gmfe": 1.4}, "by_route": {}, "worst_datasets": [],
                    "params_at_bound": [], "sensitivity": {"Intrinsic clearance": {"relative": 0.02}},
                    "recommendations": [{"parameter": "Intrinsic clearance", "issue": "weak",
                                         "action": "FIX it", "severity": "high"}],
                    "n_evals": 9, "fit_simulations": ["s1"], "series": series}
        orig = OO.run_optimization
        OO.run_optimization = fake_run
        try:
            est = {"estimate": {"Intrinsic clearance": [0.01, 3.0]}}
            r_full = _reg(False).get("osp_optimize").handler(est, ModelingSession(goal="g"))
            r_min = _reg(True).get("osp_optimize").handler(est, ModelingSession(goal="g"))
            self.assertTrue(r_full.data.get("recommendations"))   # full mode: directive advice
            self.assertFalse(r_min.data.get("recommendations"))   # minimal: none
            self.assertNotIn("NEXT:", r_min.message)              # and no NEXT: nudge in the message
            # the raw sensitivity numbers are still there in both — the model reasons itself
            self.assertIn("sensitivity", r_min.data)
        finally:
            OO.run_optimization = orig


class TestPlausibilityHardOnly(unittest.TestCase):
    def test_hard_only_keeps_impossibility_drops_opinions(self):
        params = [{"parameter": "Lipophilicity", "value": 9.0},                 # soft: out of [-2,7]
                  {"parameter": "Intrinsic clearance", "value": 3.0, "unit": "l/min"},  # soft: > blood flow
                  {"parameter": "Fraction unbound", "value": 1.5}]              # HARD: fu > 1 impossible
        full = {f["parameter"] for f in osp_score.plausibility(params)}
        hard = {f["parameter"] for f in osp_score.plausibility(params, hard_only=True)}
        self.assertEqual(full, {"Lipophilicity", "Intrinsic clearance", "Fraction unbound"})
        self.assertEqual(hard, {"Fraction unbound"})              # only the impossible one survives

    def test_hard_only_still_flags_negative_clearance(self):
        f = osp_score.plausibility([{"parameter": "Intrinsic clearance", "value": -1}], hard_only=True)
        self.assertTrue(f)                                        # <=0 is impossible, always flagged


class TestMinimalSystemPrompt(unittest.TestCase):
    def test_minimal_loads_the_minimal_playbook(self):
        import importlib.util, pathlib, sys
        p = pathlib.Path(__file__).resolve().parents[1] / "examples" / "run_llm_build.py"
        spec = importlib.util.spec_from_file_location("_rlb", p)
        m = importlib.util.module_from_spec(spec)
        sys.modules["_rlb"] = m
        spec.loader.exec_module(m)
        full = m._system_prompt(1.6, minimal=False)
        mini = m._system_prompt(1.6, minimal=True)
        self.assertIn("STOPPING RULE", mini)                     # the explicit stop condition
        self.assertIn("biologically reasonable", mini)
        self.assertNotEqual(full, mini)
        # the full playbook prescribes a staged/sweep workflow; the minimal one does not
        self.assertNotIn("osp_sweep_methods", mini)


if __name__ == "__main__":
    unittest.main()

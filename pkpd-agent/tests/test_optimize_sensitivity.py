"""Local-sensitivity ranking (identifiability evidence) - pure logic, no PK-Sim.

A parameter the objective is flat in must come out near 0 (weakly identified);
a parameter the objective responds to strongly must come out near 1.
"""

import math
import unittest

from pkpd_agent.engines.osp_optimize import _local_sensitivity


class TestLocalSensitivity(unittest.TestCase):
    def setUp(self):
        # synthetic objective in log10 space: strongly depends on p_strong,
        # ignores p_weak entirely. eval_at returns a token; log_sse computes the
        # objective from the parameter values baked into that token.
        self.names = ["p_strong", "p_weak"]
        # optimum at 1.0 for both (log10 = 0)
        self.optimized = {"p_strong": 1.0, "p_weak": 1.0}
        self.best = [0.0, 0.0]          # log10 of the optimum
        self.lx = [-3.0, -3.0]
        self.hx = [3.0, 3.0]

    def _eval_at(self, values, sims):
        return dict(values), {"ok": True}   # pass the values straight through

    def _log_sse(self, observed_sub, predicted):
        # predicted IS the values dict; objective = (log10 p_strong)^2, flat in p_weak
        ps = predicted["p_strong"]
        return (math.log10(ps) ** 2, 1)

    def test_flat_parameter_ranks_near_zero(self):
        s = _local_sensitivity(self._eval_at, self._log_sse, [], ["s1"],
                               self.optimized, self.names, self.best,
                               self.lx, self.hx)
        self.assertEqual(s["p_strong"]["relative"], 1.0)   # most influential
        self.assertLess(s["p_weak"]["relative"], 0.05)     # data ignores it
        self.assertEqual(s["p_weak"]["obj_change"], 0.0)


class TestCollinearity(unittest.TestCase):
    """Two parameters can each be individually influential yet only their SUM is
    identifiable (they trade off). One-at-a-time sensitivity misses this; the
    pairwise collinearity signal must catch it."""

    def setUp(self):
        self.names = ["a", "b"]
        self.optimized = {"a": 1.0, "b": 1.0}
        self.best = [0.0, 0.0]
        self.lx = [-3.0, -3.0]
        self.hx = [3.0, 3.0]

    def _eval_at(self, values, sims):
        return dict(values), {"ok": True}

    def test_sum_identifiable_split_is_collinear(self):
        # objective depends only on (log a + log b): moving a up and b down by
        # equal amounts leaves it unchanged -> a flat trade-off direction.
        def log_sse(obs, pred):
            return ((math.log10(pred["a"]) + math.log10(pred["b"])) ** 2, 1)
        s = _local_sensitivity(self._eval_at, log_sse, [], ["s1"],
                               self.optimized, self.names, self.best,
                               self.lx, self.hx)
        # both look influential one-at-a-time
        self.assertGreater(s["a"]["relative"], 0.5)
        self.assertGreater(s["b"]["relative"], 0.5)
        # but they are strongly collinear (opposed move is ~flat)
        self.assertGreater(s["a"]["collinearity"], 0.8)
        self.assertEqual(s["a"]["collinear_with"], "b")

    def test_independent_parameters_not_collinear(self):
        # separable objective: no flat trade-off direction between a and b.
        def log_sse(obs, pred):
            return (math.log10(pred["a"]) ** 2 + math.log10(pred["b"]) ** 2, 1)
        s = _local_sensitivity(self._eval_at, log_sse, [], ["s1"],
                               self.optimized, self.names, self.best,
                               self.lx, self.hx)
        self.assertLess(s["a"]["collinearity"], 0.5)


if __name__ == "__main__":
    unittest.main()

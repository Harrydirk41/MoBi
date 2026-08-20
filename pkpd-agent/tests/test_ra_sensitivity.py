"""Sensitivity-ranking benchmark: pool, scoring, baseline, loop tools (no LLM)."""

import unittest

from pkpd_agent.engines import ra_sensitivity as SEN
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_sensitivity_loop_tools import register_ra_sensitivity_loop_tools


class TestPool(unittest.TestCase):
    def test_pool_contains_top20_and_distractors(self):
        p = set(SEN.pool())
        self.assertTrue(set(SEN.GSA_TOP20) <= p)
        self.assertEqual(len(p), len(SEN.GSA_TOP20) + len(SEN.DISTRACTORS))

    def test_pool_is_deterministic(self):
        self.assertEqual(SEN.pool(), SEN.pool())


class TestScore(unittest.TestCase):
    def test_perfect(self):
        s = SEN.score_sensitivity(SEN.GSA_TOP20)
        self.assertEqual(s["recall"], 1.0)
        self.assertEqual(s["precision"], 1.0)
        self.assertTrue(s["beats_random"])

    def test_all_distractors(self):
        s = SEN.score_sensitivity(SEN.DISTRACTORS[:20])
        self.assertEqual(s["hit"], 0)
        self.assertFalse(s["beats_random"])

    def test_random_baseline_recall(self):
        # picking 20 of a ~50 pool blind -> expected recall 20/pool_size
        s = SEN.score_sensitivity(SEN.GSA_TOP20)     # 20 picks
        self.assertAlmostEqual(s["random_baseline_recall"], 20 / s["pool_size"], places=2)

    def test_spearman_perfect_order(self):
        s = SEN.score_sensitivity(SEN.GSA_TOP20)
        self.assertEqual(s["spearman_on_hits"], 1.0)

    def test_spearman_reversed(self):
        s = SEN.score_sensitivity(list(reversed(SEN.GSA_TOP20)))
        self.assertEqual(s["spearman_on_hits"], -1.0)

    def test_dedupe(self):
        s = SEN.score_sensitivity(["kg_FLS_Baseline", "kg_FLS_Baseline"])
        self.assertEqual(s["n_picked"], 1)


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestLoopTools(unittest.TestCase):
    def _reg(self):
        reg = ToolRegistry()
        register_ra_sensitivity_loop_tools(reg, None, {})
        return reg

    def test_registers(self):
        reg = self._reg()
        for t in ("sens_inspect", "sens_rank", "sens_finalize"):
            self.assertIn(t, reg)

    def test_inspect_gives_pool_no_key(self):
        res = self._reg().dispatch("sens_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        self.assertIn("candidate_pool", res.data)
        # the ranked GSA order must not leak - pool is alphabetical
        self.assertEqual(res.data["candidate_pool"], sorted(res.data["candidate_pool"]))

    def test_rank_and_finalize(self):
        reg = self._reg()
        sess = _FakeSession()
        reg.dispatch("sens_rank", {"ranked": SEN.GSA_TOP20}, sess)
        res = reg.dispatch("sens_finalize", {}, sess)
        self.assertTrue(res.ok)
        self.assertEqual(sess.get("sens_final")["recall"], 1.0)

    def test_finalize_requires_ranking(self):
        self.assertFalse(self._reg().dispatch("sens_finalize", {}, _FakeSession()).ok)


if __name__ == "__main__":
    unittest.main()

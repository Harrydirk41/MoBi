"""Stage-2 scope-selection: scoring and loop-tool registration (synthetic, no LLM)."""

import unittest

from pkpd_agent.engines import ra_scope as S
from pkpd_agent.tools.registry import ToolRegistry
from pkpd_agent.tools.ra_scope_loop_tools import register_ra_scope_loop_tools


class TestScoreScope(unittest.TestCase):
    def test_perfect(self):
        sc = S.score_scope(S.MODEL_NODES)
        self.assertEqual((sc.precision, sc.recall, sc.f1), (1.0, 1.0, 1.0))
        self.assertEqual(sc.missed, 0)
        self.assertEqual(sc.extra, 0)

    def test_aliases_match(self):
        sc = S.score_scope(["macrophage", "TNF-a", "IL-6", "IFNgamma"])
        self.assertEqual(sc.hit, 4)
        self.assertEqual(sc.extra, 0)

    def test_over_inclusion_hits_precision(self):
        sc = S.score_scope(S.MODEL_NODES + ["IL-2", "IL8", "NK cells"])
        self.assertEqual(sc.recall, 1.0)
        self.assertLess(sc.precision, 1.0)
        self.assertEqual(sc.extra, 3)

    def test_known_excluded_flagged(self):
        sc = S.score_scope(["IL-2", "IL-18", "gibberishXYZ"])
        self.assertIn("IL2", sc.extra_known_mediators)
        self.assertIn("IL18", sc.extra_known_mediators)
        self.assertNotIn("GIBBERISHXYZ", sc.extra_known_mediators)

    def test_partial_recall(self):
        sc = S.score_scope(S.MODEL_CELLS)          # cells only, no cytokines
        self.assertEqual(sc.hit, len(S.MODEL_CELLS))
        self.assertGreater(sc.missed, 0)
        self.assertEqual(sc.precision, 1.0)        # all correct, just incomplete

    def test_empty(self):
        sc = S.score_scope([])
        self.assertEqual((sc.hit, sc.precision, sc.recall), (0, 0.0, 0.0))

    def test_dedupes_extras(self):
        sc = S.score_scope(["IL-2", "IL2", "il 2"])   # same mediator, 3 spellings
        self.assertEqual(sc.extra, 1)


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
        register_ra_scope_loop_tools(reg, None, {})
        return reg

    def test_registers(self):
        reg = self._reg()
        for t in ("scope_inspect", "scope_propose", "scope_finalize"):
            self.assertIn(t, reg)

    def test_inspect_no_key(self):
        res = self._reg().dispatch("scope_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        blob = str(res.data).lower()
        self.assertNotIn("macro", blob)          # must not leak the cast
        self.assertNotIn("tnfa", blob)

    def test_propose_and_finalize(self):
        reg = self._reg()
        sess = _FakeSession()
        reg.dispatch("scope_propose", {"nodes": S.MODEL_NODES}, sess)
        res = reg.dispatch("scope_finalize", {}, sess)
        self.assertTrue(res.ok)
        self.assertEqual(sess.get("scope_final")["f1"], 1.0)

    def test_finalize_requires_proposal(self):
        self.assertFalse(self._reg().dispatch("scope_finalize", {}, _FakeSession()).ok)


if __name__ == "__main__":
    unittest.main()

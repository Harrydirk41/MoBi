"""Benchmark run-order priority: most-informative tasks first (adult > pediatric >
DDI, mechanism-diverse within adult), so a partial run covers what most informs the
scoreboard. The order is an informed default; these tests pin its structural rules."""

import unittest

from pkpd_agent.bench.priority import ADULT_IMPORTANCE, order_tasks, task_sort_key


class TestPriority(unittest.TestCase):
    def test_kind_groups_adult_then_pediatric_then_ddi(self):
        tasks = [{"dir": "X", "kind": "ddi"}, {"dir": "X", "kind": "pediatric"},
                 {"dir": "X", "kind": "adult"}]
        assert [t["kind"] for t in order_tasks(tasks)] == ["adult", "pediatric", "ddi"]

    def test_adult_ordered_by_importance_list(self):
        # Midazolam (index 0) must precede Digoxin, whatever the input order.
        tasks = [{"dir": "Digoxin", "kind": "adult"}, {"dir": "Midazolam", "kind": "adult"}]
        assert [t["dir"] for t in order_tasks(tasks)] == ["Midazolam", "Digoxin"]

    def test_untrustworthy_reference_deprioritised(self):
        # Propofol (TCI arms not reproduced) is last among the listed adults.
        assert ADULT_IMPORTANCE[-1] == "Propofol"
        tasks = [{"dir": "Propofol", "kind": "adult"}, {"dir": "Alfentanil", "kind": "adult"}]
        assert [t["dir"] for t in order_tasks(tasks)] == ["Alfentanil", "Propofol"]

    def test_unlisted_compound_sorts_after_listed_then_alphabetical(self):
        tasks = [{"dir": "Zzz", "kind": "adult"}, {"dir": "Aaa", "kind": "adult"},
                 {"dir": "Midazolam", "kind": "adult"}]
        # Midazolam (listed) first; then the two unlisted alphabetically
        assert [t["dir"] for t in order_tasks(tasks)] == ["Midazolam", "Aaa", "Zzz"]

    def test_order_is_stable_and_nonmutating(self):
        tasks = [{"dir": "Midazolam", "kind": "adult"}, {"dir": "Digoxin", "kind": "adult"}]
        before = list(tasks)
        order_tasks(tasks)
        assert tasks == before                       # input not mutated

    def test_sort_key_shape(self):
        k = task_sort_key("Midazolam", "adult")
        assert k == (0, 0, "Midazolam")

    def test_adult_run_order_excludes_untrustworthy_reference(self):
        from pkpd_agent.bench.priority import adult_run_order, UNTRUSTWORTHY_ADULT_REFERENCE
        ro = adult_run_order()
        assert "Propofol" in UNTRUSTWORTHY_ADULT_REFERENCE
        assert "Propofol" not in ro                       # excluded from the adult default
        assert set(ro) == set(ADULT_IMPORTANCE) - UNTRUSTWORTHY_ADULT_REFERENCE
        assert ro[0] == "Midazolam"                       # still importance-ordered

    def test_order_names_sorts_subset_and_keeps_unexcluded(self):
        from pkpd_agent.bench.priority import order_names
        # pediatric set: Propofol is NOT excluded here (only its ADULT ref is untrusted)
        out = order_names(["Propofol", "Vancomycin", "Raltegravir"])
        assert out == ["Vancomycin", "Raltegravir", "Propofol"]   # by importance index


class TestScoreboardUsesSharedOrder(unittest.TestCase):
    def test_scoreboard_orders_come_from_priority(self):
        import examples.run_osp_scoreboard as SB
        from pkpd_agent.bench import priority
        self.assertEqual(SB.DEFAULT_ORDER, priority.adult_run_order())
        self.assertNotIn("Propofol", SB.DEFAULT_ORDER)
        self.assertIn("Propofol", SB.PEDIATRIC_ORDER)     # pediatric keeps it


if __name__ == "__main__":
    unittest.main()

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


if __name__ == "__main__":
    unittest.main()

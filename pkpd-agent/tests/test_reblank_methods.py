"""Re-blanking the distribution/permeability method: a blanked snapshot must start from PK-Sim's
neutral default, not the reference model's method (which the value-blanking left in place, handing the
agent the answer's structural choice for free)."""

import unittest

from examples.reblank_methods import reblank_snapshot
from examples.verify_benchmark import _method_leaks


def _snap(part, perm="PK-Sim Standard"):
    return {"Compounds": [{"CalculationMethods": [
        f"Cellular partition coefficient method - {part}",
        f"Cellular permeability - {perm}"]}]}


class TestReblank(unittest.TestCase):
    def test_resets_non_default_partition(self):
        snap = _snap("Berezhkovskiy")
        changes = reblank_snapshot(snap)
        self.assertTrue(any(k == "partition" for k, *_ in changes))
        self.assertEqual(snap["Compounds"][0]["CalculationMethods"][0],
                         "Cellular partition coefficient method - PK-Sim Standard")

    def test_resets_non_default_permeability(self):
        snap = _snap("PK-Sim Standard", perm="Charge dependent Schmitt")
        reblank_snapshot(snap)
        self.assertEqual(snap["Compounds"][0]["CalculationMethods"][1],
                         "Cellular permeability - PK-Sim Standard")

    def test_already_neutral_is_untouched(self):
        snap = _snap("PK-Sim Standard")
        self.assertEqual(reblank_snapshot(snap), [])


class TestMethodLeakCheck(unittest.TestCase):
    def test_flags_non_default_method(self):
        self.assertTrue(_method_leaks(_snap("Schmitt")))          # leaked -> flagged
        self.assertTrue(_method_leaks(_snap("PK-Sim Standard", perm="Charge dependent Schmitt")))

    def test_neutral_default_is_clean(self):
        self.assertEqual(_method_leaks(_snap("PK-Sim Standard")), [])

    def test_reblank_then_no_leak(self):
        snap = _snap("Rodgers and Rowland", perm="Charge dependent Schmitt")
        reblank_snapshot(snap)
        self.assertEqual(_method_leaks(snap), [])                 # re-blank clears the leak


if __name__ == "__main__":
    unittest.main()

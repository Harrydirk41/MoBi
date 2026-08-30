"""L1 identifiability: the FLS steady-state target pins only the baseline rate; the shape
parameters (K's) are unidentified - any K reproduces the target once the rate re-absorbs it."""

import unittest
from examples.run_qsp_l1_calibrate import effect


class TestL1Identifiability(unittest.TestCase):
    def test_different_K_same_target_via_rate(self):
        maxes = {"IL6": 2.0, "TNFa": 1.5}
        levels = {"IL6": 289.0, "TNFa": 0.9}
        kd, target = 0.145, 2.28e7
        reproduced = []
        for kfac in (0.1, 1.0, 10.0):
            ks = {c: levels[c] * kfac for c in maxes}
            eff = effect(maxes, levels, ks)
            kg = target * kd / eff                     # 1-D fit soaks up any effect
            reproduced.append(kg / kd * eff)
        # every K guess reproduces the SAME target -> K's are not identified by it
        for r in reproduced:
            self.assertAlmostEqual(r, target, delta=target * 1e-6)

    def test_effect_is_capped(self):
        maxes = {c: 10.0 for c in "abcde"}
        levels = {c: 100.0 for c in "abcde"}
        ks = {c: 1e-9 for c in "abcde"}               # huge sum, must cap
        self.assertEqual(effect(maxes, levels, ks, cap=10.0), 11.0)


if __name__ == "__main__":
    unittest.main()

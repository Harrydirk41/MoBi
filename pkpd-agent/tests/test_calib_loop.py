"""Closed-loop calibration brain: cytokine parsing, the 1-D isolating solve, action picking.

All pure - the 1-D solver is tested against a synthetic monotonic response, the action picker
against a stub LLM. The perturbation simulations that feed the real solver live in the runner.
"""

import unittest

from pkpd_agent.engines import calib_loop as CL


class TestRegulatorCytokine(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(CL.regulator_cytokine("FLSProlif_MaxbyIL6"), "IL6")
        self.assertEqual(CL.regulator_cytokine("FLSProlif_MaxbyTNFa"), "TNFa")

    def test_no_match(self):
        self.assertEqual(CL.regulator_cytokine("kg_FLS_Baseline"), "")


class TestSolve1d(unittest.TestCase):
    def test_recovers_monotonic_param(self):
        # response(m) = 100 * m^1.3 (monotonic); target generated from true m=2.0
        true = 2.0
        resp = lambda m: 100 * m ** 1.3
        target = resp(true)
        got = CL.solve_1d(resp, target, lo=0.1, hi=20.0)
        self.assertAlmostEqual(got, true, delta=0.02)

    def test_widens_bracket_if_needed(self):
        true = 15.0
        resp = lambda m: m
        got = CL.solve_1d(resp, resp(true), lo=0.1, hi=5.0)   # target above hi -> widen
        self.assertAlmostEqual(got, true, delta=0.1)


class TestChooseExperiment(unittest.TestCase):
    def test_picks_valid_cytokine(self):
        call = lambda s, u: '{"action":"experiment","cytokine":"IL6","reason":"isolate it"}'
        a = CL.choose_experiment({"remaining": ["IL6", "TNFa"], "pinned": []}, call)
        self.assertEqual(a, {"action": "experiment", "cytokine": "IL6", "reason": "isolate it"})

    def test_invalid_pick_stops(self):
        call = lambda s, u: '{"action":"experiment","cytokine":"Nonexistent"}'
        a = CL.choose_experiment({"remaining": ["IL6"], "pinned": []}, call)
        self.assertEqual(a["action"], "stop")

    def test_stop(self):
        call = lambda s, u: '{"action":"stop"}'
        self.assertEqual(CL.choose_experiment({"remaining": ["IL6"]}, call)["action"], "stop")


if __name__ == "__main__":
    unittest.main()

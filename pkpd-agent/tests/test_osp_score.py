"""Tests for osp_score: simulation<->observed mapping, GMFE, bias, plausibility."""

import unittest

from pkpd_agent.engines.osp_cli import PredictedProfile
from pkpd_agent.engines import osp_score


def _obs(dataset, study, route, dose, t, c):
    return {"dataset": dataset, "study": study, "route": route, "dose": dose,
            "time_h": t, "conc_mg_L": c}


class TestMapping(unittest.TestCase):
    def test_dose_unit_equivalence_and_fanout(self):
        observed = [
            _obs("Ferrier healthy", "Ferrier 1985", "IV", "0.05 mg/kg", [1, 2], [10, 5]),
            _obs("Ferrier cirrhosis", "Ferrier 1985", "IV", "0.05 mg/kg", [1, 2], [12, 6]),
            # observed dose in µg/kg must match a simulation named in mg/kg
            _obs("Kharasch iv", "Kharasch 1997", "IV", "20 µg/kg", [1, 2], [1.0, 0.5]),
        ]
        profiles = [
            PredictedProfile("Ferrier 1985, Alfentanil iv 0.05 mg_kg", "Ferrier 1985",
                             "IV", "0.05 mg/kg", [1, 2], [10, 5]),
            PredictedProfile("Kharasch 1997, Alfentanil iv 0.02 mg_kg", "Kharasch 1997",
                             "IV", "0.02 mg/kg", [1, 2], [1.0, 0.5]),
        ]
        mapped, unmatched = osp_score.map_predictions(profiles, observed)
        self.assertEqual(unmatched, [])
        self.assertEqual(len(mapped), 3)                 # both Ferrier arms + Kharasch
        # 20 µg/kg matched the 0.02 mg/kg simulation
        k = [m for m in mapped if m["dataset"] == "Kharasch iv"][0]
        self.assertIn("Kharasch 1997", k["_from_simulation"])


    def test_formulation_and_food_disambiguate_same_study_dose(self):
        """Four arms of one study share (study, route, dose) and differ only in
        formulation + food. The scorer MUST pair each observation with the matching
        simulation (fed<->fed, tablet<->tablet), not collapse them onto one - else a
        fed arm is scored against a fasted simulation and the GMFE inflates."""
        observed = [
            _obs("Shah2006.Tiz.8mg.po.sd.Tab_Fed", "Shah 2006", "po", "8 mg", [1, 2], [.1, .05]),
            _obs("Shah2006.Tiz.8mg.po.sd.Cap_Fasted", "Shah 2006", "po", "8 mg", [1, 2], [.1, .05]),
        ]
        profiles = [
            PredictedProfile("Tizanidine 8mg po tablet fed", "Shah 2006", "po", "8 mg", [1, 2], [.1, .05]),
            PredictedProfile("Tizanidine 8mg po tablet fasted", "Shah 2006", "po", "8 mg", [1, 2], [9, 9]),
            PredictedProfile("Tizanidine 8mg po capsule fed", "Shah 2006", "po", "8 mg", [1, 2], [9, 9]),
            PredictedProfile("Tizanidine 8mg po capsule fasted", "Shah 2006", "po", "8 mg", [1, 2], [.1, .05]),
        ]
        mapped, unmatched = osp_score.map_predictions(profiles, observed)
        self.assertEqual(unmatched, [])
        by = {m["dataset"]: m["_from_simulation"] for m in mapped}
        self.assertEqual(by["Shah2006.Tiz.8mg.po.sd.Tab_Fed"], "Tizanidine 8mg po tablet fed")
        self.assertEqual(by["Shah2006.Tiz.8mg.po.sd.Cap_Fasted"], "Tizanidine 8mg po capsule fasted")


class TestMetrics(unittest.TestCase):
    def test_perfect_fit(self):
        observed = [_obs("D", "S", "IV", "1 mg", [1, 2, 3], [10, 5, 2.5])]
        pred = [{"dataset": "D", "time_h": [1, 2, 3], "pred_conc_mg_L": [10, 5, 2.5]}]
        r = osp_score.score_fit(observed, pred)
        self.assertAlmostEqual(r["overall"]["gmfe"], 1.0, places=6)
        self.assertAlmostEqual(r["overall"]["bias"], 1.0, places=6)

    def test_bias_direction(self):
        # prediction 2x observed -> over-prediction, bias 2.0, gmfe 2.0
        observed = [_obs("D", "S", "IV", "1 mg", [1, 2, 3], [10, 5, 2.5])]
        pred = [{"dataset": "D", "time_h": [1, 2, 3], "pred_conc_mg_L": [20, 10, 5]}]
        r = osp_score.score_fit(observed, pred)
        self.assertAlmostEqual(r["overall"]["bias"], 2.0, places=6)
        self.assertAlmostEqual(r["overall"]["gmfe"], 2.0, places=6)

    def test_interpolation_onto_observed_times(self):
        # dense prediction, sparse observation at t=1.5 (interp between 1 and 2)
        observed = [_obs("D", "S", "IV", "1 mg", [1.5], [7.5])]
        pred = [{"dataset": "D", "time_h": [1, 2], "pred_conc_mg_L": [10, 5]}]
        r = osp_score.score_fit(observed, pred)
        self.assertAlmostEqual(r["overall"]["gmfe"], 1.0, places=6)  # interp -> 7.5


class TestPlausibility(unittest.TestCase):
    def test_flags(self):
        flags = osp_score.plausibility([
            {"parameter": "Fraction unbound", "value": 1.4},
            {"parameter": "Intrinsic clearance", "value": 5.0, "unit": "l/min"},
            {"parameter": "Lipophilicity", "value": 2.0},
        ])
        msgs = {f["parameter"] for f in flags}
        self.assertIn("Fraction unbound", msgs)
        self.assertIn("Intrinsic clearance", msgs)
        self.assertNotIn("Lipophilicity", msgs)          # 2.0 is fine


if __name__ == "__main__":
    unittest.main()

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


    def test_bare_number_dose_still_constrains(self):
        """A unit-less observed dose (e.g. 300.0) must still constrain matching - it
        is mg. Otherwise every dose arm collapses onto the first same-route sim."""
        observed = [_obs("Pringle.100", "Pringle 1986", "po", 100.0, [1, 2], [.1, .05]),
                    _obs("Pringle.300", "Pringle 1986", "po", 300.0, [1, 2], [.3, .15])]
        profiles = [PredictedProfile("Mex 100 mg po", None, "po", "100 mg", [1, 2], [.1, .05]),
                    PredictedProfile("Mex 300 mg po", None, "po", "300 mg", [1, 2], [.3, .15])]
        by = {m["dataset"]: m["_from_simulation"]
              for m in osp_score.map_predictions(profiles, observed)[0]}
        self.assertEqual(by["Pringle.100"], "Mex 100 mg po")
        self.assertEqual(by["Pringle.300"], "Mex 300 mg po")

    def test_infusion_duration_and_phenotype_disambiguate(self):
        observed = [
            _obs("Wagner.iv.bolus", "Wagner 1981", "iv", "0.5 mg", [1, 2], [.1, .05]),
            _obs("Wagner.iv.3h", "Wagner 1981", "iv", "0.5 mg", [1, 2], [.9, .9]),
            _obs("Labbe.PM", "Labbe 2000", "po", "83.31 mg", [1, 2], [.1, .05]),
        ]
        profiles = [
            PredictedProfile("Digoxin iv, 0.5 mg, Bolus", "Wagner 1981", "iv", "0.5 mg", [1, 2], [.1, .05]),
            PredictedProfile("Digoxin iv, 0.5 mg, 3 h", "Wagner 1981", "iv", "0.5 mg", [1, 2], [.9, .9]),
            PredictedProfile("Mex 83.31 mg po bid EM", "Labbe 2000", "po", "83.31 mg", [1, 2], [9, 9]),
            PredictedProfile("Mex 83.31 mg po bid PM", "Labbe 2000", "po", "83.31 mg", [1, 2], [.1, .05]),
        ]
        by = {m["dataset"]: m["_from_simulation"]
              for m in osp_score.map_predictions(profiles, observed)[0]}
        self.assertEqual(by["Wagner.iv.bolus"], "Digoxin iv, 0.5 mg, Bolus")
        self.assertEqual(by["Wagner.iv.3h"], "Digoxin iv, 0.5 mg, 3 h")
        self.assertEqual(by["Labbe.PM"], "Mex 83.31 mg po bid PM")


    def test_study_is_hard_when_both_named(self):
        """A simulation named only by a study (no parseable route/dose, e.g.
        'Bovill 1984') must pair by STUDY - and an observation must NOT pile onto a
        different study's simulation just because the route matched."""
        observed = [_obs("Bovill", "Bovill 1984", "IV", "5 mg/kg", [1, 2], [.1, .05]),
                    _obs("Taverne", "Taverne 1992", "iv", "150 µg", [1, 2], [.2, .1])]
        profiles = [PredictedProfile("Bovill 1984", "Bovill 1984", None, None, [1, 2], [.1, .05]),
                    PredictedProfile("Taverne 1992", "Taverne 1992", None, None, [1, 2], [.2, .1]),
                    PredictedProfile("Willsie 2014 adults IV", "Willsie 2014", "IV", None, [1, 2], [9, 9])]
        by = {m["dataset"]: m["_from_simulation"]
              for m in osp_score.map_predictions(profiles, observed)[0]}
        self.assertEqual(by["Bovill"], "Bovill 1984")       # not the IV-named Willsie sim
        self.assertEqual(by["Taverne"], "Taverne 1992")

    def test_compact_dose_and_hyphen_study_in_names(self):
        """Dose parsed from a COMPACT sim name ('Healy1987_0.5gPer6h') and a study id
        with a hyphen ('Gepts-1987') so 0.5 g does not match the 1 g arm and a Schnider
        arm does not cross-match a Gepts simulation."""
        observed = [_obs("Healy.0.5g", "Healy 1987", "IV", "0.5 g", [1, 2], [.1, .05]),
                    _obs("Schnider.6", "Schnider 1998", "IV", "6 mg/kg", [1, 2], [.3, .1])]
        profiles = [PredictedProfile("Healy1987_0.5gPer6h", None, None, None, [1, 2], [.1, .05]),
                    PredictedProfile("Healy1987_1gBID", None, None, None, [1, 2], [9, 9]),
                    PredictedProfile("Gepts-1987 6mg_kg", None, None, None, [1, 2], [9, 9])]
        by = {m["dataset"]: m["_from_simulation"]
              for m in osp_score.map_predictions(profiles, observed)[0]}
        self.assertEqual(by["Healy.0.5g"], "Healy1987_0.5gPer6h")   # not the 1 g sim
        self.assertNotIn("Schnider.6", by)                          # not cross-matched to Gepts


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

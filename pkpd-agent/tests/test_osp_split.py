"""Building/verification study split for held-out predictive grading."""

import unittest

from pkpd_agent.engines import osp_split as SP


def _obs(dataset, study, route):
    return {"dataset": dataset, "study": study, "route": route, "dose": "1 mg"}


def _snap(sim_names, links, classifications=None):
    sims = []
    for nm in sim_names:
        sims.append({"Name": nm,
                     "OutputMappings": [{"ObservedData": ds}
                                        for ds, s in links.items() if s == nm]})
    snap = {"Simulations": sims}
    if classifications:
        snap["SimulationClassifications"] = classifications
    return snap


class TestSplit(unittest.TestCase):
    def test_holds_out_whole_studies_iv_stays_in_building(self):
        # 6 PO studies + 1 IV study; IV must stay in building, ~1/3 of PO held out
        obs, links, sims = [], {}, []
        for i, st in enumerate(["Aa 2001", "Bb 2002", "Cc 2003", "Dd 2004", "Ee 2005", "Ff 2006"]):
            ds = f"{st}.po"
            obs.append(_obs(ds, st, "po")); links[ds] = f"sim{i}"; sims.append(f"sim{i}")
        obs.append(_obs("Gg 2007.iv", "Gg 2007", "IV")); links["Gg 2007.iv"] = "simIV"; sims.append("simIV")
        sp = SP.split_studies(_snap(sims, links), obs)
        self.assertEqual(sp["method"], "study-holdout")
        self.assertIn("Gg 2007.iv", sp["building"])       # IV anchor stays in building
        self.assertTrue(sp["verification"])               # something held out
        # no study is split across building and verification
        held = set(sp["held_out_studies"])
        for ds in sp["verification"]:
            self.assertIn(SP.osp_score._study_token(ds), held)
        for ds in sp["building"]:
            self.assertNotIn(SP.osp_score._study_token(ds), held)

    def test_routeless_dataset_stays_in_building_never_held_out(self):
        # a linked dataset with no route/dose metadata is a shaky held-out target
        # (the blind grade may get no curve — the alfentanil Kharasch2012 rows), but
        # it IS simulatable for building, so it stays in building and is never held out.
        obs, links, sims = [], {}, []
        for i, st in enumerate(["Aa 2001", "Bb 2002", "Cc 2003", "Dd 2004", "Ee 2005"]):
            ds = f"{st}.po"; obs.append(_obs(ds, st, "po"))
            links[ds] = f"sim{i}"; sims.append(f"sim{i}")
        # a route-less study, linked to a simulation but with route/dose = None
        for nm in ("Kharasch2012_alone_a", "Kharasch2012_alone_b"):
            obs.append({"dataset": nm, "study": "Kharasch2012", "route": None, "dose": None})
            links[nm] = nm; sims.append(nm)
        sp = SP.split_studies(_snap(sims, links), obs)
        self.assertIn("Kharasch2012_alone_a", sp["building"])      # kept for building
        self.assertNotIn("Kharasch2012_alone_a", sp["verification"])  # never held out
        self.assertNotIn("Kharasch2012", sp["held_out_studies"])

    def test_iv_only_compound_still_gets_holdout(self):
        # every study carries an IV arm (no PO-only study), as for alfentanil. The
        # split must still hold out a whole IV-containing study while keeping >=1 IV
        # study in building, so the compound is graded blind instead of not at all.
        obs, links, sims = [], {}, []
        studies = ["Aa 2001", "Bb 2002", "Cc 2003", "Dd 2004", "Ee 2005", "Ff 2006"]
        for i, st in enumerate(studies):
            for route in ("IV", "po"):                 # each study has both arms
                ds = f"{st}.{route}"
                obs.append(_obs(ds, st, route))
                links[ds] = f"sim{i}"; sims.append(f"sim{i}")
        sp = SP.split_studies(_snap(sims, links), obs)
        self.assertEqual(sp["method"], "study-holdout-iv")
        self.assertTrue(sp["verification"])                       # graded blind
        self.assertTrue(sp["held_out_studies"])
        # at least one IV study remains in building to anchor disposition
        held = set(sp["held_out_studies"])
        self.assertLess(len(held), len(studies))
        self.assertGreaterEqual(len(sp["building"]), len(sp["verification"]))

    def test_building_stays_majority(self):
        # a single data-rich study must not dominate the held-out set
        obs, links, sims = [], {}, []
        for st in ["Aa 2001", "Bb 2002", "Cc 2003", "Dd 2004", "Ee 2005"]:
            for j in range(2):
                ds = f"{st}.arm{j}"; obs.append(_obs(ds, st, "po"))
                links[ds] = f"sim_{st}"; sims.append(f"sim_{st}")
        # one giant study
        for j in range(20):
            ds = f"Big 2009.arm{j}"; obs.append(_obs(ds, "Big 2009", "po"))
            links[ds] = "sim_big"; sims.append("sim_big")
        sp = SP.split_studies(_snap(sims, links), obs)
        self.assertGreaterEqual(len(sp["building"]), len(sp["verification"]))

    def test_osp_verification_tag_is_used(self):
        obs = [_obs("A.po", "Aa 2001", "po"), _obs("B.po", "Bb 2002", "po")]
        links = {"A.po": "simA", "B.po": "simB"}
        snap = _snap(["simA", "simB"], links,
                     classifications=[{"Name": "model verification", "Classifiables": ["simB"]}])
        sp = SP.split_studies(snap, obs)
        self.assertEqual(sp["method"], "osp-tag")
        self.assertEqual(sp["verification"], ["B.po"])
        self.assertEqual(sp["building"], ["A.po"])

    def test_too_few_studies_no_holdout(self):
        obs = [_obs("A.iv", "Aa 2001", "IV"), _obs("B.iv", "Bb 2002", "IV")]
        links = {"A.iv": "simA", "B.iv": "simB"}
        sp = SP.split_studies(_snap(["simA", "simB"], links), obs)
        self.assertEqual(sp["method"], "none-too-few")
        self.assertEqual(sp["verification"], [])
        self.assertEqual(len(sp["building"]), 2)

    def test_unmapped_datasets_excluded(self):
        # a dataset the model maps to no simulation is in neither split
        obs = [_obs(f"S{i} 200{i}.po", f"S{i} 200{i}", "po") for i in range(6)]
        obs.append(_obs("Ghost 2099.po", "Ghost 2099", "po"))   # not linked
        links = {o["dataset"]: f"sim{i}" for i, o in enumerate(obs[:6])}
        sp = SP.split_studies(_snap([f"sim{i}" for i in range(6)], links), obs)
        allsplit = set(sp["building"]) | set(sp["verification"])
        self.assertNotIn("Ghost 2099.po", allsplit)


if __name__ == "__main__":
    unittest.main()

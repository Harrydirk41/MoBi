"""propose_vpop_set: general, model-agnostic Vpop-driver selection from a full param list.

The engine must select the varied set by REASONING over parameter names, with no disease
vocabulary baked in - so these tests drive it on a NON-RA (tumor-shaped) parameter list
and a stub 'analyst' LLM, and assert the plumbing (batch filtering, aggregation, no
hallucinated names) rather than any RA specifics.
"""

import json
import unittest

from pkpd_agent.engines import llm_tasks as LT


class TestProposeVpopSet(unittest.TestCase):
    def test_general_selection_no_disease_vocab(self):
        # a tumor-immunology-shaped parameter set - nothing RA about it
        params = [
            {"name": "kg_TumorCell", "value": 0.1},
            {"name": "kmig_TCell", "value": 0.2},
            {"name": "ksec_VEGF", "value": 0.3},
            {"name": "F_IFNg", "value": 4.0},
            {"name": "CL_drugX", "value": 1.0},      # PK clearance -> fixed
            {"name": "KD_drugX", "value": 1e-9},     # binding constant -> fixed
        ]

        def stub(system, user):
            # a competent analyst: vary growth/migration/secretion/amplification, drop PK/KD
            return json.dumps({
                "selected": ["kg_TumorCell", "kmig_TCell", "ksec_VEGF", "F_IFNg"],
                "categories": {"kg_TumorCell": "cell growth",
                               "kmig_TCell": "cell migration",
                               "ksec_VEGF": "mediator secretion",
                               "F_IFNg": "pathway amplification"},
                "rationale": "varied growth/migration/secretion/amplification; excluded PK+binding"})

        r = LT.propose_vpop_set(params, stub)
        self.assertEqual(set(r["selected"]),
                         {"kg_TumorCell", "kmig_TCell", "ksec_VEGF", "F_IFNg"})
        self.assertEqual(r["n_candidates"], 6)
        self.assertEqual(r["n_selected"], 4)
        self.assertEqual(r["categories"]["ksec_VEGF"], "mediator secretion")
        self.assertNotIn("CL_drugX", r["selected"])   # fixed constants excluded
        self.assertNotIn("KD_drugX", r["selected"])
        self.assertTrue(r["rationale"])

    def test_drops_hallucinated_names(self):
        params = [{"name": "a"}, {"name": "b"}]

        def stub(system, user):
            return json.dumps({"selected": ["a", "NOT_A_REAL_PARAM"],
                               "categories": {"NOT_A_REAL_PARAM": "x"}, "rationale": ""})

        r = LT.propose_vpop_set(params, stub)
        self.assertEqual(r["selected"], ["a"])          # name not in the list is dropped
        self.assertNotIn("NOT_A_REAL_PARAM", r["categories"])

    def test_aggregates_across_batches(self):
        params = [{"name": f"p{i}"} for i in range(5)]

        def stub(system, user):
            # each batch: select whatever p-names are in it (echo the first one)
            payload = json.loads(user.split("PARAMETERS:\n", 1)[1])
            picks = [p["name"] for p in payload]
            return json.dumps({"selected": picks, "categories": {}, "rationale": "all"})

        r = LT.propose_vpop_set(params, stub, batch=2)  # forces 3 batches
        self.assertEqual(r["n_selected"], 5)
        self.assertEqual(set(r["selected"]), {f"p{i}" for i in range(5)})

    def test_survives_a_bad_reply(self):
        params = [{"name": "a"}, {"name": "b"}]

        def stub(system, user):
            return "not json at all"

        r = LT.propose_vpop_set(params, stub)     # must not raise
        self.assertEqual(r["selected"], [])
        self.assertEqual(r["n_candidates"], 2)


if __name__ == "__main__":
    unittest.main()

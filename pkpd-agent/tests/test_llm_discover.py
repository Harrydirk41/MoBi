"""Structure-discovery brain: candidate-set construction, symptom-to-edge ranking, scoring.
Pure - a stub call stands in for the LLM."""

import unittest

from pkpd_agent.engines import llm_discover as D


class TestCandidateSet(unittest.TestCase):
    def test_includes_true_edge_and_distractors(self):
        pool = [("A", "B"), ("C", "D"), ("E", "F"), ("G", "H")]
        cand = D.candidate_set(("IL12", "IL6"), pool, n_distractors=2, seed=1)
        self.assertIn(("IL12", "IL6"), cand)
        self.assertEqual(len(cand), 3)                 # true + 2 distractors
        self.assertNotIn(("IL12", "IL6"), pool)        # pool untouched by dedup


class TestMaxbyEdges(unittest.TestCase):
    def test_parses_edges_and_knobs(self):
        species = {"FLS", "IL6", "IL1b", "TNFa"}
        names = ["FLSProlif_MaxbyIL6", "IL6SecFLS_MaxbyIL1b", "kg_FLS_Baseline",
                 "FLSProlif_MaxbyTNFa"]
        e = D.maxby_edges(names, species)
        self.assertEqual(e[("IL6", "FLS")], "FLSProlif_MaxbyIL6")   # IL6 -> FLS proliferation
        self.assertEqual(e[("IL1b", "IL6")], "IL6SecFLS_MaxbyIL1b")  # IL1b -> IL6 secretion
        self.assertEqual(e[("TNFa", "FLS")], "FLSProlif_MaxbyTNFa")
        self.assertNotIn(("FLS", "FLS"), e)            # kg_FLS_Baseline has no Maxby


class TestRankCandidates(unittest.TestCase):
    def test_llm_ranking_respected_and_completed(self):
        cands = [("IL12", "IL6"), ("A", "B"), ("C", "D")]
        call = lambda s, u: ('{"ranking": [{"src":"IL12","dst":"IL6"},{"src":"A","dst":"B"}], '
                             '"reason":"IL6 low -> missing activator"}')
        r = D.rank_candidates([{"species": "IL6", "direction": "too low"}], cands, call)
        self.assertEqual(r[0], ("IL12", "IL6"))
        self.assertIn(("C", "D"), r)                   # dropped one appended back
        self.assertEqual(len(r), 3)

    def test_rank_of_true(self):
        ranking = [("A", "B"), ("IL12", "IL6"), ("C", "D")]
        self.assertEqual(D.rank_of_true(ranking, ("IL12", "IL6")), 2)
        self.assertEqual(D.rank_of_true(ranking, ("X", "Y")), 4)   # absent -> len+1


if __name__ == "__main__":
    unittest.main()

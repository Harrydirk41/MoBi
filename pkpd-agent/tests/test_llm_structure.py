"""LLM structure extractor: parsing, aggregation, and the regression comparison.

Uses a STUB call_fn (canned JSON) so the extraction/aggregation/compare logic is tested
without an API. The real LLM is swapped in only in run_llm_extract.
"""

import unittest

from pkpd_agent.engines import llm_structure as LS
from pkpd_agent.engines.qsp_model import QSPModel, get_spec

VANTAGE_RA_SPEC = get_spec("vantage_ra")
from tests.test_qsp_model import _DAS28, _BIO, _NONBIO


def _fixture_model():
    data = {
        "species": [{"name": n} for n in _BIO + _NONBIO],
        "rules": [{"rule": _DAS28},
                  {"rule": "Pro_FLSProlif_effect = min(10,MM(TNFa,a)+MM(IL6,a))"}],
        "parameters": [],
    }
    return QSPModel(data, VANTAGE_RA_SPEC)


class TestJsonParse(unittest.TestCase):
    def test_strips_fences(self):
        self.assertEqual(LS._parse_json('```json\n[{"a":1}]\n```'), [{"a": 1}])

    def test_finds_object_in_prose(self):
        self.assertEqual(LS._parse_json('here: {"x": 2} done'), {"x": 2})


class TestExtraction(unittest.TestCase):
    def test_classify(self):
        call = lambda s, u: '{"biology":["FLS","TNFa"],"drug":["MTXDrug_Central"],' \
                            '"readout":["DAS28_CRP"]}'
        c = LS.classify_species(["FLS", "TNFa", "MTXDrug_Central", "DAS28_CRP"], call)
        self.assertEqual(c["biology"], ["FLS", "TNFa"])
        self.assertEqual(c["readout"], ["DAS28_CRP"])

    def test_extract_edges_filters_to_bio(self):
        call = lambda s, u: '[{"source":"TNFa","target":"FLS","sign":1},' \
                            '{"source":"TNFa","target":"FLS","sign":1},' \
                            '{"source":"Aspirin","target":"FLS","sign":1}]'
        e = LS.extract_edges(["TNFa", "FLS", "IL6"], ["rule1"], call)
        self.assertEqual(len(e), 1)                       # deduped, Aspirin dropped
        self.assertEqual(e[0], {"source": "TNFa", "sign": 1, "target": "FLS"})

    def test_extract_readout(self):
        call = lambda s, u: '{"drivers":["FLS","Macrophages","NotASpecies"]}'
        d = LS.extract_readout(["FLS", "Macrophages"], ["DAS28_CRP = FLS+Macrophages"], call)
        self.assertEqual(set(d), {"FLS", "Macrophages"})


class TestRegression(unittest.TestCase):
    def test_perfect_extraction_scores_high(self):
        m = _fixture_model()
        # a "perfect" LLM returns exactly the deterministic edges/readout/nodes
        result = {
            "nodes": list(m.nodes),
            "edges": [{"source": e.source, "sign": e.sign, "target": e.target}
                      for e in m.edges],
            "readout_drivers": list(m.readout_drivers),
        }
        cmp = LS.compare_to_model(result, m)
        self.assertEqual(cmp["edges_signed"]["recall"], 1.0)
        self.assertEqual(cmp["readout"]["recall"], 1.0)
        self.assertEqual(cmp["nodes"]["recall"], 1.0)

    def test_missing_edges_lower_recall(self):
        m = _fixture_model()
        result = {"nodes": list(m.nodes), "edges": [], "readout_drivers": []}
        cmp = LS.compare_to_model(result, m)
        self.assertEqual(cmp["edges_signed"]["recall"], 0.0)
        self.assertTrue(cmp["edge_disagreements"]["deterministic_only"])


if __name__ == "__main__":
    unittest.main()

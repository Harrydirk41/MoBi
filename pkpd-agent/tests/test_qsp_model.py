"""Model-agnostic QSPModel derivation: nodes/edges/params/readout from a network.json dump.

Uses a fixture built from the REAL Vantage RA structure (species list, sample rules, the
actual DAS28 formula, param names/units) to prove the derivation reproduces the hardcoded
RA answer keys - the backward-compat regression that validates Problem A on this model.
"""

import unittest

from pkpd_agent.engines.qsp_model import QSPModel, VANTAGE_RA_SPEC, get_spec

_DAS28 = ("Synovium.DAS28_CRP = 2*FLS^2.5/(1.3E7^2.5+FLS^2.5)+0.5*Endothelial^2.5/"
          "(4.2e7^2.5+Endothelial^2.5)+1.5*Th1^2.5/(4e6^2.5+Th1^2.5)+0.5*Th17^2.5/"
          "(1e5^2.5+Th17^2.5)-0.5*Treg^2.5/(2E6^2.5+Treg^2.5)+0.5*CTL^2.5/"
          "(1.3E5^2.5+CTL^2.5)+1.5*BCells^2.5/(3E6^2.5+BCells^2.5)+1*PlasmaCells^2.5/"
          "(1.8E6^2.5+PlasmaCells^2.5)+2.5*Macrophages^2.5/(2.2e7^2.5+Macrophages^2.5)")

_BIO = ["FLS", "TNFa", "Endothelial", "IL1b", "IL6", "VEGF", "Macrophages", "BCells",
        "Th1", "Treg", "CTL", "IL17", "IL10", "TGFb", "IFNg", "AutoAb", "GMCSF", "MIP3",
        "MCP1", "CAM", "BAFF", "Th17", "PlasmaCells", "IL12", "IL23", "RANTES"]
_NONBIO = ["TNFa_Ada", "DAS28_CRP", "DAS28_BL", "ACR_Perc", "ACR20", "MTX_NonResp",
           "MTXDrug_Central", "AdaDrug_Central", "TCZDrug_Synovium", "ana_Dose"]


def _model():
    data = {
        "species": [{"name": n} for n in _BIO + _NONBIO],
        "rules": [
            {"rule": _DAS28},
            {"rule": "ACR_Perc = 100*delta_DAS28_CRP/DAS28_BL"},
            {"rule": "delta_DAS28_CRP = DAS28_BL-DAS28_CRP"},
            {"rule": "Pro_FLSProlif_effect = min(10,MM(TNFa,a,b,c)+MM(IL6,a,b,c))"},
            {"rule": "Anti_EndoInflux_effect = min(0.9,MM(TGFb,a,b,c))"},
            {"rule": "Pro_MacroProlif_effect = min(10,MM(GMCSF,a,b,c)+MM(TNFa,a,b,c))"},
        ],
        "parameters": [
            {"name": "kd_FLS_Baseline", "units": "1/day", "value": 0.1},
            {"name": "FLSProlif_MaxbyTNFa", "units": "dimensionless", "value": 1.5},
            {"name": "IL6SecFLS_MaxbyIL1b", "units": "nanogram/(molecule*day)", "value": 1.0},
            {"name": "KD_TCZ", "units": "M", "value": 2.5e-12},
            {"name": "Pro_IL6Sec_byMacro_effect", "units": "dimensionless", "value": 1},
        ],
    }
    return QSPModel(data, VANTAGE_RA_SPEC)


class TestNodeDerivation(unittest.TestCase):
    def test_excludes_drugs_and_readouts(self):
        m = _model()
        self.assertEqual(len(m.nodes), 26)
        self.assertIn("Macrophages", m.nodes)
        for bad in ("DAS28_CRP", "MTXDrug_Central", "ACR20", "TNFa_Ada"):
            self.assertNotIn(bad, m.nodes)


class TestReadoutDerivation(unittest.TestCase):
    def test_nine_cell_drivers(self):
        m = _model()
        self.assertEqual(set(m.readout_drivers),
                         {"FLS", "Endothelial", "Th1", "Th17", "Treg", "CTL", "BCells",
                          "PlasmaCells", "Macrophages"})


class TestEdgeDerivation(unittest.TestCase):
    def test_rule_and_abbreviation_edges(self):
        m = _model()
        sset = {(e.source, e.sign, e.target) for e in m.edges}
        self.assertIn(("TNFa", 1, "FLS"), sset)
        self.assertIn(("TGFb", -1, "Endothelial"), sset)     # abbreviation Endo->Endothelial
        self.assertIn(("GMCSF", 1, "Macrophages"), sset)     # abbreviation Macro->Macrophages
        self.assertIn(("Macrophages", 1, "IL6"), sset)       # from the _by param name


class TestParamDerivation(unittest.TestCase):
    def test_unit_split(self):
        m = _model()
        self.assertEqual(sum(p.physiological() for p in m.params), 2)   # 1/day, M
        self.assertEqual(sum(p.model_scaling() for p in m.params), 1)   # ng/molecule/day
        self.assertEqual(sum(p.dimensionless() for p in m.params), 2)


class TestMatcher(unittest.TestCase):
    def test_free_text_resolves(self):
        m = _model()
        self.assertEqual(m.resolve("macrophage"), "Macrophages")
        self.assertEqual(m.resolve("Th1 cell"), "Th1")
        self.assertEqual(m.resolve("IL-6"), "IL6")
        self.assertEqual(m.resolve("CCL2"), "MCP1")
        self.assertIsNone(m.resolve("aspirin"))

    def test_score_node_set(self):
        m = _model()
        s = m.score_node_set(["macrophage", "FLS", "aspirin"], ["Macrophages", "FLS", "IL6"])
        self.assertEqual(s["hit"], 2)
        self.assertIn("ASPIRIN", s["extra"])

    def test_resolve_all_splits_compound(self):
        m = _model()
        self.assertEqual(m.resolve_all("B cell / plasma cell"), {"BCells", "PlasmaCells"})
        self.assertEqual(m.resolve_all("Macrophage (synovial)"), {"Macrophages"})
        self.assertEqual(m.resolve_all("CD4 T cell (Th1, Th17)"), {"Th1", "Th17"})
        self.assertEqual(m.resolve_all("osteoclast"), set())

    def test_score_node_set_credits_compound(self):
        m = _model()
        s = m.score_node_set(["B cell / plasma cell", "Macrophage (synovial)"],
                             ["BCells", "PlasmaCells", "Macrophages"])
        self.assertEqual(s["hit"], 3)                 # all three, from 2 compound entries


class TestInferSpec(unittest.TestCase):
    def _data(self):
        return {
            "species": [{"name": n} for n in _BIO + _NONBIO],
            "rules": [{"rule": _DAS28},
                      {"rule": "ACR_Perc = 100*delta_DAS28_CRP/DAS28_BL"},
                      {"rule": "delta_DAS28_CRP = DAS28_BL-DAS28_CRP"},
                      {"rule": "Pro_FLSProlif_effect = min(10,MM(TNFa,a,b,c)+MM(IL6,a,b,c))"},
                      {"rule": "Anti_EndoInflux_effect = min(0.9,MM(TGFb,a,b,c))"}],
            "parameters": [{"name": "kd_FLS_Baseline", "units": "1/day", "value": 0.1}],
        }

    def test_inferred_recovers_nodes_and_readout(self):
        from pkpd_agent.engines.qsp_model import infer_spec
        data = self._data()
        m = QSPModel(data, infer_spec(data, "auto"))
        # all 26 biological nodes recovered (may add a drug-conjugate false positive)
        for bio in _BIO:
            self.assertIn(bio, m.nodes)
        self.assertEqual(set(m.readout_drivers),
                         {"FLS", "Endothelial", "Th1", "Th17", "Treg", "CTL", "BCells",
                          "PlasmaCells", "Macrophages"})

    def test_short_token_not_matched_midword(self):
        # 'ACR' must not exclude 'Macrophages' (M-acr-o)
        from pkpd_agent.engines.qsp_model import infer_spec
        data = self._data()
        m = QSPModel(data, infer_spec(data, "auto"))
        self.assertIn("Macrophages", m.nodes)


class TestSpecs(unittest.TestCase):
    def test_get_spec(self):
        # spec is loaded from projects/vantage_ra/spec.json (data, not a code literal)
        self.assertEqual(get_spec("ra"), VANTAGE_RA_SPEC)
        self.assertEqual(get_spec("ra").name, "Vantage RA")
        self.assertEqual(len(get_spec("ra").gsa_top), 20)
        with self.assertRaises(KeyError):
            get_spec("nonexistent_model")


if __name__ == "__main__":
    unittest.main()

"""Agent-parameter Vpop: severity-knob selection (ksec + non-marginal kprolif), response frac,
and Vpop xlsx writing. Pure - the MATLAB engine parts are exercised only on the user's machine."""

import os
import tempfile
import unittest

from examples import run_qsp_agent_vpop as V


class TestSelectSeverity(unittest.TestCase):
    def test_keeps_ksec_and_nonmarginal_kprolif_only(self):
        params = [{"name": "ksec_IL6", "value": 1600.0},
                  {"name": "kprolif_Th1", "value": 0.39},
                  {"name": "kprolif_FLS", "value": 0.044},        # FLS marginal -> excluded
                  {"name": "M_FLS_prolif_IL6", "value": 2.0},     # edge strength, not a severity knob
                  {"name": "kcl_IL6", "value": 5.5},              # clearance, not a knob
                  {"name": "ksec_TNFa", "value": 0.0}]            # non-positive -> skipped
        spec, chosen = V.select_severity_params(params, {"FLS", "Macrophages"}, 2.0)
        self.assertEqual(chosen, ["ksec_IL6", "kprolif_Th1"])
        self.assertIn("ksec_IL6,800,3200,log", spec)

    def test_span_is_log_uniform_around_nominal(self):
        spec, _ = V.select_severity_params([{"name": "ksec_X", "value": 10.0}], set(), 4.0)
        self.assertEqual(spec, "ksec_X,2.5,40,log")


class TestFrac(unittest.TestCase):
    def test_fraction_and_mask(self):
        self.assertEqual(V._frac({"ACR20": [1, 0, 1, 1, float("nan")]}, "ACR20"), (75.0, 4))
        self.assertEqual(V._frac({"ACR20": [1, 0, 1], "M": [1, 1, 0]}, "ACR20", mask="M"), (50.0, 2))
        self.assertEqual(V._frac({"ACR20": []}, "ACR20"), (None, 0))


class TestWriteXlsx(unittest.TestCase):
    def test_header_then_patient_rows(self):
        p = os.path.join(tempfile.gettempdir(), "test_vp.xlsx")
        V._write_vpop_xlsx(p, ["a", "b"], [{"a": 1, "b": 2}, {"a": 3, "b": 4}])
        import openpyxl
        ws = openpyxl.load_workbook(p).active
        got = [[c.value for c in r] for r in ws.iter_rows()]
        self.assertEqual(got, [["a", "b"], [1, 2], [3, 4]])


class TestMarginalFromProject(unittest.TestCase):
    def test_ra_marginal_cells(self):
        self.assertEqual(V._marginal_cells("ra"), {"FLS", "Macrophages", "PlasmaCells"})


if __name__ == "__main__":
    unittest.main()

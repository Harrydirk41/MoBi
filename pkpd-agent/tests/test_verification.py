import unittest

from pkpd_agent.verification import run_gates


class TestVerification(unittest.TestCase):
    def test_failed_minimization_blocks(self):
        data = {"minimization_successful": False, "condition_number": 10}
        findings = run_gates("pharmpy_fit", data, None)
        self.assertTrue(any(f.level == "block" and f.gate == "fit_convergence"
                            for f in findings))

    def test_high_condition_number_warns(self):
        data = {"minimization_successful": True, "condition_number": 5000}
        findings = run_gates("pharmpy_fit", data, None)
        self.assertTrue(any(f.gate == "fit_conditioning" and f.level == "warn"
                            for f in findings))

    def test_high_rse_warns(self):
        data = {"minimization_successful": True,
                "relative_standard_errors": {"KA": 0.8}}
        findings = run_gates("pharmpy_fit", data, None)
        self.assertTrue(any(f.gate == "parameter_precision" for f in findings))

    def test_negative_concentration_blocks(self):
        data = {"all_values_finite": True, "min_concentration": -0.1,
                "mass_balance_residual": 0.0}
        findings = run_gates("osp_simulate", data, None)
        self.assertTrue(any(f.level == "block" and f.gate == "physical_sanity"
                            for f in findings))

    def test_mass_balance_residual_warns(self):
        data = {"all_values_finite": True, "min_concentration": 0.0,
                "mass_balance_residual": 1e-2}
        findings = run_gates("osp_simulate", data, None)
        self.assertTrue(any(f.gate == "mass_balance" for f in findings))

    def test_clean_fit_has_no_findings(self):
        data = {"minimization_successful": True, "condition_number": 20,
                "relative_standard_errors": {"CL": 0.1, "V": 0.12}}
        self.assertEqual(run_gates("pharmpy_fit", data, None), [])


if __name__ == "__main__":
    unittest.main()

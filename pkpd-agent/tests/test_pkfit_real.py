"""Tests for the REAL pkfit engine and the decision policy.

Skipped automatically if numpy/scipy are not installed, so the stdlib-only
core test suite still runs everywhere.
"""

import unittest

try:
    import numpy  # noqa: F401
    import scipy  # noqa: F401
    HAVE_SCI = True
except Exception:  # noqa: BLE001
    HAVE_SCI = False


@unittest.skipUnless(HAVE_SCI, "numpy/scipy required for the real engine")
class TestPKFitEngine(unittest.TestCase):
    def setUp(self):
        from pkpd_agent.engines.pkfit import PKFitEngine, simulate_dataset
        self.eng = PKFitEngine()
        self.data, self.truth = simulate_dataset(seed=7)

    def test_recovers_known_truth(self):
        fit = self.eng.fit(self.data, "1cpt_oral")
        self.assertTrue(fit["minimization_successful"])
        est = fit["parameter_estimates"]
        # within ~25% of the simulation truth (naive-pooled, small sample)
        self.assertAlmostEqual(est["CL"], self.truth["CL_ref"], delta=0.25 * self.truth["CL_ref"])
        self.assertAlmostEqual(est["V"], self.truth["V_ref"], delta=0.25 * self.truth["V_ref"])
        self.assertAlmostEqual(est["Ka"], self.truth["Ka"], delta=0.25 * self.truth["Ka"])

    def test_reports_precision_and_conditioning(self):
        fit = self.eng.fit(self.data, "1cpt_oral")
        self.assertIn("relative_standard_errors", fit)
        self.assertTrue(fit["relative_standard_errors"])
        self.assertLess(fit["condition_number"], 1e3)

    def test_covariate_lowers_ofv_significantly(self):
        base = self.eng.fit(self.data, "1cpt_oral")
        cov = self.eng.fit(self.data, "1cpt_oral",
                           covariate={"param": "CL", "cov": "WT", "ref": 70})
        self.assertGreater(base["ofv"] - cov["ofv"], 3.84)  # LRT significant

    def test_vpc_runs_and_covers(self):
        fit = self.eng.fit(self.data, "1cpt_oral")
        vpc = self.eng.vpc(self.data, fit, n_sim=200)
        self.assertGreater(vpc["pct_observations_within_90_pi"], 70)


@unittest.skipUnless(HAVE_SCI, "numpy/scipy required for the real engine")
class TestPolicyOnRealEngine(unittest.TestCase):
    def test_full_workflow_picks_1cpt_and_keeps_covariate(self):
        from pkpd_agent.config import AgentConfig
        from pkpd_agent.loop import DecisionLoop
        from pkpd_agent.policies import PharmacometricPolicy

        loop = DecisionLoop(config=AgentConfig(mock=False),
                            policy=PharmacometricPolicy())
        session = loop.run("Build a popPK model.")

        self.assertTrue(session.finished)
        # 2-compartment fit is degenerate here -> must have raised a BLOCK
        self.assertTrue(any(o.blocked for o in session.observations))
        # final model is the 1-compartment model with the WT covariate kept
        fit_ids = [o.content.get("fit_id") for o in session.observations
                   if o.tool == "pkfit_fit"]
        self.assertIn("1cpt_oral+WT_on_CL", fit_ids)


if __name__ == "__main__":
    unittest.main()

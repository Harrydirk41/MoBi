import unittest

from pkpd_agent.doctor import run_probes, summary


class TestDoctor(unittest.TestCase):
    def test_probes_never_raise_and_report(self):
        probes = run_probes(rscript=None)
        self.assertTrue(probes)
        # every probe has the required fields and a boolean status
        for p in probes:
            self.assertIsInstance(p.ok, bool)
            self.assertTrue(p.name)
            self.assertTrue(p.unlocks)

    def test_summary_renders(self):
        text = summary(run_probes(None))
        self.assertIn("engine health check", text)
        self.assertIn("checks passing", text)

    def test_missing_rscript_is_reported_not_raised(self):
        probes = run_probes(rscript="/definitely/not/a/real/rscript")
        rprobe = [p for p in probes if p.name == "R:Rscript"][0]
        self.assertFalse(rprobe.ok)


if __name__ == "__main__":
    unittest.main()

"""Tests for the headless PK-Sim CLI engine (parsing/matching only; no PK-Sim).

The subprocess calls need PKSim.CLI.exe on Windows and are exercised by
examples/osp_run.py on the user's machine. Here we lock the pure logic: CSV
parsing + unit conversion, simulation-name parsing, and the mapping of
simulations to observed datasets used by the runner.
"""

import csv
import os
import tempfile
import unittest

from pkpd_agent.engines.osp_cli import OSPCli

CLI = OSPCli(pksim_cli_path="/nonexistent")

HEADER = ["IndividualId", "Time [min]",
          "Organism|PeripheralVenousBlood|Alfentanil|Plasma (Peripheral Venous Blood) [µmol/l]",
          "Organism|VenousBlood|Plasma|Alfentanil|Concentration in container [µmol/l]",
          "Organism|Lumen|Alfentanil|Fraction of oral drug mass absorbed into mucosa"]


class TestResultParsing(unittest.TestCase):
    def _write(self, d, rows):
        simdir = os.path.join(d, "Ferrier 1985, Alfentanil iv 0.05 mg_kg")
        os.makedirs(simdir)
        p = os.path.join(simdir, "Ferrier 1985, Alfentanil iv 0.05 mg_kg-Results.csv")
        with open(p, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(HEADER)
            for t, c in rows:
                w.writerow([0, t, c, c * 0.9, 0])
        return d

    def test_unit_and_time_conversion(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(d, [(0, 0), (3, 0.6376957), (6, 0.4858887)])
            prof = CLI._parse_results(d, mol_weight=416.52)[0]
            # 3 min -> 0.05 h ; 0.6376957 µmol/L * 416.52/1000 = 0.26561 mg/L
            self.assertAlmostEqual(prof.time_h[1], 0.05, places=6)
            self.assertAlmostEqual(prof.conc_mg_L[1], 0.6376957 * 416.52 / 1000, places=6)

    def test_picks_peripheral_venous_plasma_column(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(d, [(0, 0), (3, 1.0)])
            idx = CLI._pick_conc_column(HEADER)
            self.assertIn("PeripheralVenousBlood", HEADER[idx])
            self.assertIn("Plasma (Peripheral Venous Blood)", HEADER[idx])

    def test_population_mean_collapse(self):
        # two IndividualIds at the same time -> mean
        with tempfile.TemporaryDirectory() as d:
            simdir = os.path.join(d, "Pop, Alfentanil iv 1 mg")
            os.makedirs(simdir)
            p = os.path.join(simdir, "Pop, Alfentanil iv 1 mg-Results.csv")
            with open(p, "w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(HEADER)
                w.writerow([0, 3, 1.0, 0.9, 0])
                w.writerow([1, 3, 3.0, 2.7, 0])       # mean at t=3 is 2.0
            prof = CLI._parse_results(d, mol_weight=1000.0)[0]
            self.assertAlmostEqual(prof.conc_mg_L[0], 2.0, places=6)  # 2.0 µmol/L*1 = 2 mg/L


class TestSimNameParsing(unittest.TestCase):
    def test_names(self):
        cases = {
            "Ferrier 1985, Alfentanil iv 0.05 mg_kg": ("Ferrier 1985", "IV", "0.05 mg/kg"),
            "Kharasch 2011b, Alfentanil PO 4 mg": ("Kharasch 2011b", "PO", "4 mg"),
            "Meistelman 1987, 20µg_kg, adult male individual":
                ("Meistelman 1987", None, "20 µg/kg"),
        }
        for name, expect in cases.items():
            self.assertEqual(CLI._parse_sim_name(name), expect)

    def test_route_across_underscore(self):
        # route token separated by underscores (Kharasch2012_..._IV)
        _, route, _ = CLI._parse_sim_name("Kharasch2012_Alfentanil_alone_IV")
        self.assertEqual(route, "IV")


if __name__ == "__main__":
    unittest.main()

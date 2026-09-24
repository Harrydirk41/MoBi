"""PBPK Copilot server: the local web app over the agent.

The agent run itself needs PKSim.CLI + an Anthropic key (not in CI), so these
tests cover the wiring that does not: model listing, the leave-one-out library
endpoint, the UI page, and that a run with no key streams a clean error rather
than crashing.
"""
import os
import unittest

try:
    from fastapi.testclient import TestClient
    _HAVE_FASTAPI = True
except Exception:                                       # noqa: BLE001
    _HAVE_FASTAPI = False

_LIB = os.path.join(os.path.dirname(__file__), "..", "..", "OSP-PBPK-Model-Library")


@unittest.skipUnless(_HAVE_FASTAPI, "fastapi not installed")
@unittest.skipUnless(os.path.isdir(_LIB), "OSP library not present")
class TestCopilotServer(unittest.TestCase):
    def setUp(self):
        import sys
        pkg = os.path.join(os.path.dirname(__file__), "..")
        for p in (pkg, os.path.join(pkg, "examples")):
            if p not in sys.path:
                sys.path.insert(0, p)
        from copilot import server
        self.c = TestClient(server.create_app())

    def test_index_serves_ui(self):
        r = self.c.get("/")
        self.assertEqual(r.status_code, 200)
        self.assertIn("PBPK Copilot", r.text)

    def test_models_listed(self):
        r = self.c.get("/api/models")
        self.assertEqual(r.status_code, 200)
        comps = [m["compound"] for m in r.json()]
        self.assertIn("Triazolam", comps)
        self.assertIn("Digoxin", comps)

    def test_library_is_leave_one_out(self):
        r = self.c.get("/api/library?compound=Triazolam&same_type=true")
        self.assertEqual(r.status_code, 200)
        names = [m["compound"] for m in r.json().get("models", [])]
        self.assertNotIn("Triazolam", names)       # never its own model
        self.assertIn("Midazolam", names)          # a CYP3A4 analogue

    def test_rw_projects_and_series(self):
        r = self.c.get("/api/rw/projects")
        self.assertEqual(r.status_code, 200)
        projs = r.json()
        if not projs:
            self.skipTest("pbpk-realworld not generated")
        self.assertIn("Triazolam", projs)
        s = self.c.get("/api/rw/series?project=Triazolam&kind=test").json()["series"]
        self.assertTrue(s and all(len(e["points"]) >= 2 for e in s))

    def test_rw_report_is_redacted(self):
        r = self.c.get("/api/rw/report?project=Triazolam&kind=test")
        self.assertEqual(r.status_code, 200)
        md = r.json()["markdown"]
        if not md:
            self.skipTest("pbpk-realworld not generated")
        self.assertNotIn("GMFE", md)               # no held-out grade
        self.assertNotIn("Rodgers and Rowland", md)  # no chosen method leaked

    def test_rw_context_report_is_full_disclosure(self):
        # the support-context viewer opens OTHER projects' complete reports
        r = self.c.get("/api/rw/report?project=Midazolam&kind=context")
        self.assertEqual(r.status_code, 200)
        md = r.json()["markdown"]
        if not md:
            self.skipTest("pbpk-realworld not generated")
        self.assertIn("GMFE", md)                  # full model report, not redacted
        self.assertTrue(self.c.get(
            "/api/rw/series?project=Midazolam&kind=context").json()["series"])

    def test_rw_topology_parses_structure(self):
        r = self.c.get("/api/rw/topology?project=Tizanidine")
        if r.status_code == 404:
            self.skipTest("snapshot not present")
        t = r.json()
        self.assertTrue(t.get("distribution"))          # partition method
        self.assertTrue(t.get("renal"))                 # has GFR
        self.assertTrue(t.get("oral"))                  # tablet formulation
        self.assertIn("CYP1A2", [m["enzyme"] for m in t["metabolism"]])
        self.assertEqual(self.c.get("/api/rw/topology?project=NoSuch").status_code, 404)

    def test_rw_simulate_degrades_gracefully(self):
        # without PKSim.CLI the reference-fit run returns a clean error, never 500
        import os as _os
        from copilot import server
        prev = _os.environ.pop("PKPD_PKSIM_CLI", None)
        server._REFFIT_CACHE.pop("Alfentanil", None)
        try:
            r = self.c.get("/api/rw/simulate?project=Alfentanil")
            self.assertEqual(r.status_code, 200)
            j = r.json()
            # either it ran (CLI present) with series, or it reported why it couldn't
            self.assertIn("ok", j)
            if not j["ok"]:
                self.assertTrue(j.get("error"))
            miss = self.c.get("/api/rw/simulate?project=NoSuchProject").json()
            self.assertFalse(miss["ok"])
        finally:
            if prev is not None:
                _os.environ["PKPD_PKSIM_CLI"] = prev

    def test_rw_figure_serves_png_and_blocks_traversal(self):
        p = ("images/006_section_results-and-discussion/008_section_diagnostics-plots/"
             "2_gof_plot_predictedVsObserved.png")
        r = self.c.get("/api/rw/figure", params={"project": "Alfentanil", "path": p})
        if r.status_code == 404:
            self.skipTest("OSP figure not present")
        self.assertEqual(r.headers.get("content-type"), "image/png")
        self.assertEqual(r.content[:4], b"\x89PNG")
        bad = self.c.get("/api/rw/figure",
                         params={"project": "Alfentanil", "path": "../../../etc/passwd"})
        self.assertEqual(bad.status_code, 404)

    def test_context_projects_narrows_library(self):
        # an explicit allowlist restricts the reference library the run would see
        from pkpd_agent.engines import reference_library as RL
        from copilot import server
        f = server._find("Triazolam")
        allow = ["Midazolam", "Alprazolam"]
        lib = RL.library_for_snapshot(f["snapshot"], same_type=False, full=False, only=allow)
        got = {m["compound"] for m in lib["models"]}
        self.assertTrue(got.issubset(set(allow)))
        self.assertNotIn("Triazolam", got)

    def test_deliverable_packaging_and_zip(self):
        import io as _io
        import shutil
        import zipfile
        from copilot import server
        edits = {"calculation_methods": {"partition": "Schmitt"},
                 "estimate": {"Lipophilicity": [-2, 3]},
                 "parameters": {"Lipophilicity": 2.1, "CLspec": 4.4}}
        givens = [{"parameter": "Fraction unbound", "value": 0.17, "unit": "",
                   "source": "PMID x", "provenance": "given"}]
        split = {"building": ["A"], "held_out_studies": ["B"]}
        obs = [{"dataset": "A IV", "route": "IV", "dose": "1 mg",
                "time_h": [0.5, 1.0], "conc_mg_L": [0.01, 0.008]}]
        out = server._write_deliverable("Triazolam", "x", edits, 1.42, givens, split, obs)
        try:
            names = set(os.listdir(out))
            self.assertLessEqual({"README.md", "parameter_table.csv",
                                  "adopted_model.json", "observed_curves.csv"}, names)
            table = open(os.path.join(out, "parameter_table.csv"), encoding="utf-8").read()
            self.assertIn("given", table)          # extracted given
            self.assertIn("fitted", table)         # optimizer fit, tagged distinctly
            r = self.c.get("/api/deliverable?compound=Triazolam")
            self.assertEqual(r.status_code, 200)
            z = zipfile.ZipFile(_io.BytesIO(r.content))
            self.assertTrue(any(n.endswith("parameter_table.csv") for n in z.namelist()))
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def test_deliverable_missing_is_404(self):
        r = self.c.get("/api/deliverable?compound=NoSuchCompound")
        self.assertEqual(r.status_code, 404)

    def test_access_token_gates_when_set(self):
        import os as _os
        from copilot import server
        _os.environ["PBPK_COPILOT_TOKEN"] = "sekret123"
        try:
            c = TestClient(server.create_app())
            self.assertEqual(c.get("/api/models").status_code, 401)     # no token
            self.assertEqual(c.get("/api/models?token=nope").status_code, 401)
            r = c.get("/?token=sekret123")
            self.assertEqual(r.status_code, 200)
            self.assertIn("pbpk_token", r.headers.get("set-cookie", ""))
            self.assertEqual(c.get("/api/models").status_code, 200)     # via cookie
        finally:
            _os.environ.pop("PBPK_COPILOT_TOKEN", None)

    def test_obs_detail_surfaces_what_agent_saw(self):
        from copilot.server import _obs_detail
        d = _obs_detail("osp_inspect", {
            "objective": "Build PBPK", "parameters_to_determine": ["Lipophilicity"],
            "literature_physicochemical": [{"parameter": "fu", "value": 0.1, "unit": ""}]})
        self.assertIn("objective", d)
        self.assertIn("to determine (fit these)", d)
        o = _obs_detail("osp_options", {"editable_parameters": [
            {"name": "Lipophilicity", "tier": "estimate"}],
            "calculation_methods": {"partition": {"options": ["Schmitt"]}}})
        self.assertIn("editable parameters", o)
        r = _obs_detail("osp_read_report", {"report_markdown": "x" * 5000})
        self.assertTrue(list(r.values())[0].endswith("…"))     # trimmed
        self.assertIsNone(_obs_detail("osp_inspect", {}))      # nothing to show

    def test_no_token_is_open(self):
        # default (no token set) stays open — unchanged local behavior
        self.assertEqual(self.c.get("/api/models").status_code, 200)

    def test_run_without_key_streams_error_not_crash(self):
        import json
        prev = os.environ.pop("ANTHROPIC_API_KEY", None)
        try:
            rid = self.c.post("/api/run", json={"compound": "Triazolam"}).json()["run_id"]
            types = []
            with self.c.stream("GET", f"/api/stream/{rid}") as s:
                for line in s.iter_lines():
                    if line and line.startswith("data:"):
                        ev = json.loads(line[5:])
                        types.append(ev.get("type"))
                        if ev.get("type") == "end":
                            break
            self.assertIn("error", types)
            self.assertEqual(types[-1], "end")
        finally:
            if prev is not None:
                os.environ["ANTHROPIC_API_KEY"] = prev


if __name__ == "__main__":
    unittest.main()

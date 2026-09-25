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

    def test_catalog_lists_full_action_space(self):
        c = self.c.get("/api/catalog").json()
        self.assertEqual(len(c["partition_methods"]), 5)
        self.assertEqual(len(c["permeability_methods"]), 3)
        self.assertGreaterEqual(len(c["process_types"]), 15)
        self.assertTrue(c["parameters"]["estimate"])       # fittable set
        self.assertTrue(c["parameters"]["constant"])       # never-fit set
        self.assertTrue(all(m.get("description") for m in c["partition_methods"]))
        self.assertEqual(c["counts"]["fittable"], len(c["parameters"]["estimate"]))
        self.assertEqual(len(c["ddi_types"]), 6)                     # inhibition/induction kinetics
        areas = {a["area"]: a["status"] for a in c["pksim_areas"]}   # the wider PK-Sim/MoBi map
        self.assertIn("MoBi (open modeling)", areas)
        self.assertEqual(areas["MoBi (open modeling)"], "mobi")

    def test_built_view_gives_tunable_sliders(self):
        v = self.c.get("/api/built?compound=Tizanidine").json()
        self.assertTrue(v["methods"]["partition"]["current"])
        self.assertTrue(v["params"])
        for p in v["params"]:                                # each slider has value in range
            self.assertLessEqual(p["lo"], p["value"])
            self.assertLessEqual(p["value"], p["hi"])
        self.assertTrue(any("clearance" in p["name"].lower() for p in v["params"]))

    def test_pediatric_and_ddi_tasks_listed(self):
        ms = self.c.get("/api/models").json()
        kinds = {m["compound"]: m["kind"] for m in ms}
        peds = [m for m in ms if m["kind"] == "pediatric"]
        ddis = [m for m in ms if m["kind"] == "ddi"]
        self.assertTrue(peds and ddis)                       # existing tutorials wired in
        self.assertTrue(all(m.get("dir") for m in ms))       # each task resolves to a base dir
        # a DDI task (no clinical data) fails cleanly on the web runner, not a 500
        did = ddis[0]["compound"]
        j = self.c.post("/api/try", json={"compound": did, "edits": {}}).json()
        self.assertFalse(j["ok"])
        self.assertTrue(j.get("error"))

    def test_stdout_router_isolates_threads(self):
        import io as _io
        import threading
        from copilot.server import _install_router
        r = _install_router()
        outs = {}

        def worker(name):
            buf = _io.StringIO()
            r.register(buf)
            try:
                for i in range(4):
                    print(f"{name}{i}")
            finally:
                r.flush(); outs[name] = buf.getvalue(); r.unregister()
        ts = [threading.Thread(target=worker, args=(n,)) for n in "ABCD"]
        [t.start() for t in ts]; [t.join() for t in ts]
        for n, v in outs.items():                      # no cross-thread leakage
            self.assertTrue(all(line.startswith(n) for line in v.strip().splitlines()))

    def test_model_view_is_editable_action_space(self):
        r = self.c.get("/api/model?compound=Tizanidine")
        self.assertEqual(r.status_code, 200)
        v = r.json()
        self.assertEqual(len(v["methods"]["partition"]["options"]), 5)
        self.assertTrue(v["methods"]["partition"]["current"])
        tiers = {p["tier"] for p in v["parameters"]}
        self.assertTrue({"estimate", "measured_soft", "constant"} & tiers)
        self.assertEqual(self.c.get("/api/model?compound=NoSuch").status_code, 404)

    def test_try_model_degrades_without_cli(self):
        r = self.c.post("/api/try", json={"compound": "Tizanidine",
                                          "edits": {"calculation_methods": {"partition": "Schmitt"}}})
        self.assertEqual(r.status_code, 200)
        j = r.json()
        self.assertIn("ok", j)
        if not j["ok"]:
            self.assertTrue(j.get("error"))
        self.assertEqual(self.c.post("/api/try", json={"compound": "NoSuch"}).status_code, 404)

    def test_grade_endpoint_held_out(self):
        r = self.c.post("/api/grade", json={"compound": "Triazolam", "edits": {}})
        self.assertEqual(r.status_code, 200)
        self.assertIn("ok", r.json())                       # graceful (needs CLI to actually grade)
        self.assertEqual(self.c.post("/api/grade", json={"compound": "NoSuch"}).status_code, 404)

    def test_rw_topology_parses_structure(self):
        r = self.c.get("/api/rw/topology?project=Tizanidine")
        if r.status_code == 404:
            self.skipTest("snapshot not present")
        t = r.json()
        self.assertTrue(t.get("distribution"))          # partition method
        self.assertTrue(t.get("renal"))                 # has GFR
        self.assertTrue(t.get("oral"))                  # tablet formulation
        self.assertIn("CYP1A2", [m["enzyme"] for m in t["metabolism"]])
        # hover-info fields: physchem + per-process clearance parameters
        self.assertIsNotNone(t["physchem"]["lipophilicity"])
        self.assertIsNotNone(t["physchem"]["molecular_weight"])
        self.assertTrue(t["metabolism"][0]["params"])           # e.g. intrinsic clearance
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

    def test_report_compares_adopted_model_to_reference(self):
        # the post-build report grades the adopted model against the reference answer
        r = self.c.post("/api/report", json={"compound": "Alprazolam", "edits": {
            "calculation_methods": {"partition": "Rodgers and Rowland",
                                    "permeability": "PK-Sim Standard"},
            "add_processes": [{"type": "metabolization_mm", "molecule": "CYP3A4"}],
            "parameters": {"Lipophilicity": 2.5, "kcat@CYP3A4": 10.0}}})
        self.assertEqual(r.status_code, 200)
        j = r.json()
        comp = j.get("comparison")
        self.assertTrue(comp)
        by = {s["aspect"]: s for s in comp["structure"]}
        self.assertTrue(by["partition method"]["match"])          # matched the reference
        self.assertTrue(by["clearing molecules"]["match"])
        lip = next(p for p in comp["parameters"] if p["name"] == "Lipophilicity")
        self.assertIsNotNone(lip["fold"])                         # fold vs reference computed
        self.assertIn(lip["verdict"], ("recovered", "close", "off", "far"))
        # curves degrade gracefully without a CLI
        self.assertIn("in_sample", j)
        # the LLM interpretation field is present (None without an API key — best-effort)
        self.assertIn("interpretation", j)
        self.assertEqual(self.c.post("/api/report", json={"compound": "Nope"}).status_code, 404)

    def test_run_persists_and_reloads(self):
        import os as _os
        from copilot import server
        events = [{"type": "meta", "building": 5},
                  {"type": "decision", "text": "hi", "calls": []},
                  {"type": "done", "best_gmfe": 1.42, "best_edits": {"parameters": {"x": 1}},
                   "web_lookups": [], "blind": True},
                  {"type": "end"}]
        p = {"compound": "Triazolam", "model": "claude-sonnet-5"}
        path = server._run_record_path("Triazolam")
        try:
            server._persist_run("Triazolam", events, p)
            rec = server._load_run("Triazolam")
            self.assertEqual(rec["in_sample_gmfe"], 1.42)
            self.assertEqual(len(rec["events"]), 4)
            # /api/runs replays it
            r = self.c.get("/api/runs/Triazolam")
            self.assertEqual(r.status_code, 200)
            self.assertEqual(r.json()["best_edits"], {"parameters": {"x": 1}})
            # /api/runs/<c>/grade persists the held-out grade
            self.assertTrue(self.c.post("/api/runs/Triazolam/grade",
                                        json={"held_out": 1.77}).json()["ok"])
            self.assertEqual(server._load_run("Triazolam")["held_out_gmfe"], 1.77)
            # /api/models now reflects the persisted run (built + saved grade)
            m = next(x for x in self.c.get("/api/models").json() if x["compound"] == "Triazolam")
            self.assertEqual(m["status"], "done")
            self.assertEqual((m.get("saved") or {}).get("held_out"), 1.77)
            self.assertEqual(self.c.get("/api/runs/NoSuchCompound").status_code, 404)
        finally:
            if _os.path.isfile(path):
                _os.remove(path)

    def test_hard_mode_swaps_to_structure_blind_snapshot(self):
        # structure-blind run picks the hard_blanked snapshot (mechanism stripped)
        # and emits a `hard` meta before it needs a key.
        import json
        prev = os.environ.pop("ANTHROPIC_API_KEY", None)
        try:
            rid = self.c.post("/api/run", json={"compound": "Triazolam",
                                                "hard": True}).json()["run_id"]
            metas = []
            with self.c.stream("GET", f"/api/stream/{rid}") as s:
                for line in s.iter_lines():
                    if line and line.startswith("data:"):
                        ev = json.loads(line[5:])
                        if ev.get("type") == "meta":
                            metas.append(ev)
                        if ev.get("type") == "end":
                            break
            self.assertTrue(any(m.get("hard") for m in metas))     # structure-blind engaged
        finally:
            if prev is not None:
                os.environ["ANTHROPIC_API_KEY"] = prev

    def test_cancel_sets_flag_and_404s_unknown(self):
        from copilot import server
        self.assertEqual(self.c.post("/api/cancel/nosuchrun").status_code, 404)
        rid = "unit_cancel_run"
        server._RUNS[rid] = {"queue": None, "done": False}
        try:
            r = self.c.post(f"/api/cancel/{rid}")
            self.assertEqual(r.status_code, 200)
            self.assertTrue(r.json()["ok"])
            self.assertTrue(server._RUNS[rid]["cancel"])     # loop's should_stop reads this
        finally:
            server._RUNS.pop(rid, None)

    def test_clear_cache_removes_reffit_only(self):
        import os as _os
        from copilot import server
        server._REFFIT_CACHE["Zz"] = [1]
        # drop a fake .reffit.json into a compound's benchmark dir, if the library exists
        libs = [d for d in (_os.path.isdir(server._LIB) and _os.listdir(server._LIB) or [])
                if _os.path.isdir(_os.path.join(server._LIB, d, "benchmark"))]
        made = None
        if libs:
            made = _os.path.join(server._LIB, libs[0], "benchmark", ".reffit.json")
            existed = _os.path.exists(made)
            if not existed:
                open(made, "w").write("[]")
        n = server._clear_cache()
        self.assertEqual(server._REFFIT_CACHE, {})           # in-memory cache cleared
        if made:
            self.assertFalse(_os.path.exists(made))          # the disk cache file is gone
            self.assertGreaterEqual(n, 1)


if __name__ == "__main__":
    unittest.main()

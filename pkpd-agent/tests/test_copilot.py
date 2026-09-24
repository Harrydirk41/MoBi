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

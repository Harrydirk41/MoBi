"""Proof that the 2-7 tasks run off the CONFIG, not RA scaffolding.

Same machinery, different tasks.json: a psoriasis-flavoured config (PASI endpoints,
renamed columns, different drugs, no RA vocabulary) is driven through summarize_run,
the scorers, and every loop-tool family on SYNTHETIC run data - no API, no MATLAB. If
these pass, the tasks are generic: the vocabulary comes entirely from the config.

Two configs are exercised:
  * ALT  - a full psoriasis-shaped config (all five task families apply).
  * MIN  - a config with NO second-line / refractory block, to prove a member of the
           class that lacks that structure still loads and the other tasks still run.
"""

import unittest

from pkpd_agent.engines import qsp_config as C
from pkpd_agent.engines import qsp_tasks as T
from pkpd_agent.tools.registry import ToolRegistry


# --- a psoriasis-shaped config: DIFFERENT endpoints, columns, drugs, no RA words --- #
_ALT = {
    "name": "Psoriasis QSP",
    "project_aliases": ["pso"],
    "disease": "plaque psoriasis; severity PASI",
    "severity_readout": "PASI",
    "readout_desc": "PASI75 / PASI90 response",
    "readout_states": ["PASI75", "PASI90", "L2_IR", "L2_PASI75", "L2_PASI90",
                       "PASI_score", "PASI_base"],
    "run_columns": {
        "patient": "pid",
        "first_line": {"PASI75": "PASI75", "PASI90": "PASI90"},
        "subgroup_flag": "L2_IR",
        "second_line": {"PASI75": "L2_PASI75", "PASI90": "L2_PASI90"},
        "severity": {"baseline": "PASI_base", "readout": "PASI_read"},
    },
    "timeline": {"baseline_day": 0, "first_line_readout_day": 84,
                 "second_line_readout_day": 168},
    "drugs": {"SEC": {"drug": "secukinumab", "modality": "anti-IL-17A mAb",
                      "mechanism": "neutralizes IL-17A", "doses": ["SEC_300mg_Q4W"]}},
    "vpop_drivers": {"F_IL17": {"nominal": 10.0, "span": [0.1, 90],
                                "meaning": "IL-17 amplification"}},
    "vpop_target": {"mean": 12.0, "sd": 5.0, "band": [6, 40]},
    "fit_params": {"KD_SEC": {"unit": "M", "reference": 1e-10, "meaning": "IL-17 binding",
                              "search_range": [1e-12, 1e-8], "log_scale": True}},
    "design_targets": {"F_IL17": {"pathway": "IL-17A", "analogue": "secukinumab",
                                  "note": "central in psoriasis"}},
    "clinical_trials": {"SEC": {"trial": "UNCOVER-2", "population": "moderate-severe",
                                "weeks": {"12": {"drug": {"PASI75": 77.0, "PASI90": 54.0}}}}},
    "refractory_target": {"trial": "some refractory PsO trial", "PASI75": 60.0},
    "validate_arms": {"prior_therapies": ["SEC_300mg_Q4W"], "test_arm": "SEC_300mg_Q4W"},
    "flagship_protocol": {}, "calibrated_arms": [],
    "trial_objective": "predict second-line PASI response",
    "fit_default_arm": "SEC_300mg_Q4W", "fit_target": {"drug": "SEC", "week": 12},
    "design_background": "",
}

_MIN = {**_ALT, "name": "Single-readout PsO", "project_aliases": ["psomin"],
        "subgroup": None}
_MIN = dict(_ALT)
_MIN["run_columns"] = {**_ALT["run_columns"], "subgroup_flag": "", "second_line": {}}
_MIN["refractory_target"] = {}
_MIN["validate_arms"] = {}


def _cfg(d):
    return C.config_from_dict(d)


def _run(**cols):
    return {"columns": cols}


class _FakeSession:
    def __init__(self):
        self._d = {}

    def get(self, k, default=None):
        return self._d.get(k, default)

    def put(self, k, v):
        self._d[k] = v


class TestSummarizeOffAltColumns(unittest.TestCase):
    def test_reads_pasi_columns_not_ra(self):
        cfg = _cfg(_ALT)
        res = _run(pid=[1, 2, 3, 4],
                   PASI75=[1, 0, 1, 1], PASI90=[1, 0, 0, 1],
                   L2_IR=[0, 1, 1, 1],
                   L2_PASI75=[0, 1, 0, 1], L2_PASI90=[0, 1, 0, 0],
                   PASI_base=[12, 20, 15, 30], PASI_read=[4, 18, 14, 6])
        s = cfg.summarize_run(res)
        self.assertEqual(s["first_line"]["PASI75"], 75.0)      # 3/4
        self.assertEqual(s["first_line"]["PASI90"], 50.0)      # 2/4
        self.assertEqual(s["second_line"]["n_subgroup"], 3)    # L2_IR==1
        self.assertEqual(s["second_line"]["PASI75"], 66.7)     # 2/3 among IR
        self.assertEqual(s["severity"]["baseline_mean"], 19.25)
        # the summary carries NO RA keys
        self.assertNotIn("ACR20", s["first_line"])
        self.assertNotIn("das28", s)

    def test_trial_target_and_scoring_generic(self):
        cfg = _cfg(_ALT)
        tgt = cfg.trial_target("SEC", 12, "raw")
        self.assertEqual(tgt, {"PASI75": 77.0, "PASI90": 54.0})   # not ACR
        sc = T.score_flagship({"PASI75": 70.0, "PASI90": 50.0}, tgt)
        self.assertEqual(sc["n_endpoints"], 2)


class TestAllToolFamiliesRunOffAltConfig(unittest.TestCase):
    """Register every task family with the psoriasis config; the inspect output must
    speak PASI/secukinumab and contain NO RA vocabulary."""

    def _reg(self, register, extra=None):
        cfg = _cfg(_ALT)
        reg = ToolRegistry()
        ctx = {"cfg": cfg, "sb": None, "vpop": "V"}
        ctx.update(extra or {})
        register(reg, None, ctx)
        return reg

    def _no_ra(self, blob: str):
        low = blob.lower()
        for w in ("das28", "acr20", "acr50", "mtx", "tcz", "tocilizumab", "rheumatoid",
                  "radiate", "methotrexate"):
            self.assertNotIn(w, low, f"RA vocab '{w}' leaked into a tool output")

    def test_trial(self):
        from pkpd_agent.tools.qsp_trial_loop_tools import register_qsp_trial_loop_tools
        reg = self._reg(register_qsp_trial_loop_tools)
        for n in ("trial_inspect", "trial_run", "trial_finalize"):
            self.assertIn(n, reg)
        res = reg.dispatch("trial_inspect", {}, _FakeSession())
        self.assertTrue(res.ok)
        blob = str(res.data)
        self.assertIn("secukinumab", blob)
        self._no_ra(blob)

    def test_fit(self):
        from pkpd_agent.tools.qsp_fit_loop_tools import register_qsp_fit_loop_tools
        reg = self._reg(register_qsp_fit_loop_tools,
                        {"arm": "SEC_300mg_Q4W", "target": {"PASI75": 77.0}})
        res = reg.dispatch("fit_inspect", {}, _FakeSession())
        self.assertEqual(res.data["parameters_to_fit"][0]["name"], "KD_SEC")
        self._no_ra(str(res.data))

    def test_vpop(self):
        from pkpd_agent.tools.qsp_vpop_loop_tools import register_qsp_vpop_loop_tools
        reg = self._reg(register_qsp_vpop_loop_tools)
        res = reg.dispatch("vpop_inspect", {}, _FakeSession())
        self.assertIn("PASI", str(res.data))
        self.assertEqual(res.data["clinical_target"]["mean"], 12.0)
        self._no_ra(str(res.data))

    def test_design(self):
        from pkpd_agent.tools.qsp_design_loop_tools import register_qsp_design_loop_tools
        reg = self._reg(register_qsp_design_loop_tools,
                        {"sbproj": "p", "vpop": "v"})
        res = reg.dispatch("design_inspect", {}, _FakeSession())
        self.assertEqual(len(res.data["targetable_pathways"]), 1)   # F_IL17
        self._no_ra(str(res.data))

    def test_validate(self):
        from pkpd_agent.tools.qsp_validate_loop_tools import \
            register_qsp_validate_loop_tools
        reg = self._reg(register_qsp_validate_loop_tools)
        res = reg.dispatch("validate_inspect", {}, _FakeSession())
        self.assertIn("comparator", res.data)
        self._no_ra(str(res.data))


class TestMinimalConfigNoSecondLine(unittest.TestCase):
    """A class member WITHOUT the second-line / refractory structure still loads and the
    first-line tasks still run; the second-line arm is simply empty."""

    def test_loads_and_first_line_only(self):
        cfg = _cfg(_MIN)
        res = _run(pid=[1, 2], PASI75=[1, 0], PASI90=[1, 0],
                   PASI_base=[12, 20], PASI_read=[4, 18])
        s = cfg.summarize_run(res)
        self.assertEqual(s["first_line"]["PASI75"], 50.0)
        self.assertEqual(s["second_line"]["n_subgroup"], 0)     # no subgroup flag
        self.assertEqual(s["second_line"], {"n_subgroup": 0})   # nothing else to report

    def test_validate_tools_still_register(self):
        from pkpd_agent.tools.qsp_validate_loop_tools import \
            register_qsp_validate_loop_tools
        reg = ToolRegistry()
        register_qsp_validate_loop_tools(reg, None, {"cfg": _cfg(_MIN), "sb": None,
                                                     "vpop": "V"})
        for n in ("validate_inspect", "validate_run", "validate_finalize"):
            self.assertIn(n, reg)


class TestSameModelDifferentEndpoint(unittest.TestCase):
    """The shipped vantage_ra_rem project: SAME model/columns as vantage_ra, but the
    config tracks DAS28 remission instead of ACR. Proves endpoint is config-driven on
    the real (committed) config files, not just synthetic dicts."""

    def _run(self):
        return _run(patient=[1, 2, 3, 4],
                    ACR20=[1, 0, 0, 1], ACR50=[1, 0, 0, 0], ACR70=[0, 0, 0, 0],
                    Rem=[1, 0, 0, 1], MTX_NonResp=[0, 1, 1, 1],
                    TCZ_ACR20=[0, 1, 0, 1], TCZ_Rem=[0, 1, 1, 0],
                    DAS28_base=[5, 5, 5, 6], DAS28_read=[2.4, 4.9, 4.7, 2.5])

    def test_ra_reports_acr_ra_rem_reports_remission(self):
        ra = C.get("ra").summarize_run(self._run())
        rem = C.get("ra_rem").summarize_run(self._run())
        # ACR config exposes ACR endpoints; remission config exposes only remission
        self.assertIn("ACR20", ra["first_line"])
        self.assertNotIn("ACR20", rem["first_line"])
        self.assertEqual(list(rem["first_line"].keys()), ["n", "remission"])
        self.assertEqual(rem["second_line"], {"n_subgroup": 3, "remission": 66.7})

    def test_fit_target_follows_endpoint(self):
        # same trial machinery, endpoint from config: ACR keys vs a remission key
        ra_t = C.get("ra").trial_target(**C.get("ra").fit_target)
        rem_t = C.get("ra_rem").trial_target(**C.get("ra_rem").fit_target)
        self.assertIn("ACR20", ra_t)
        self.assertEqual(rem_t, {"remission": 38.0})


class _FakeCohortEngine:
    """Returns a synthetic cohort table (as sb_cohort.m would) so the multi-anchor
    tool path is testable without MATLAB."""
    def cohort_multi_arm(self, spec, arms_spec, bday, rday, n, seed, states=None):
        import random
        random.seed(1)
        sev, mtx = [], []
        for _ in range(int(n)):
            s = random.uniform(2.0, 8.0)
            sev.append(s)
            mtx.append(1 if s < 4.5 else 0)          # milder respond to MTX
        return {"columns": {"sample": list(range(int(n))), "sev_base": sev, "MTX": mtx},
                "matlab_log": ""}


class TestVpopMultiAnchorTool(unittest.TestCase):
    def test_multi_anchor_tool_matches_config_anchors(self):
        from pkpd_agent.tools.qsp_vpop_loop_tools import register_qsp_vpop_loop_tools
        reg = ToolRegistry()
        register_qsp_vpop_loop_tools(reg, None, {
            "cfg": C.get("ra"), "sb": _FakeCohortEngine(), "n_pool": 200})
        self.assertIn("vpop_select_multi", reg)        # registered because RA has anchors
        res = reg.dispatch("vpop_select_multi",
                           {"bounds": {"F_TNFa": [0.1, 50, "log"]}}, _FakeSession())
        self.assertTrue(res.ok, res.message)
        keys = {a["key"] for a in res.data["anchors"]}
        self.assertIn("severity", keys)               # baseline moment anchor
        self.assertIn("MTX", keys)                     # the MTX response-rate anchor
        self.assertGreater(res.data["effective_sample_size"], 5)


if __name__ == "__main__":
    unittest.main()

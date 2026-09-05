"""The clinical calibration ('training') step: the agent fits the transplanted model's free MTX
drug-effect to the FIRST-LINE MTX arm only. These tests are pure - the MATLAB engine (sb) is mocked -
and they assert two things: the bounded search finds the effect-scale that matches the training
target, and the held-out second-line arm is NEVER read while fitting."""

import json
import math
import unittest

from examples import run_qsp_full_pipeline as P


class TestBisectScale(unittest.TestCase):
    def test_finds_scale_matching_a_monotone_target(self):
        # a monotone-increasing readout in the scale; target is reachable inside the bracket
        def evaluate(s):
            return 50.0 * (1.0 - math.exp(-s))          # 0 -> 50, increasing
        scale, y, trace = P._bisect_scale(evaluate, target=33.7, lo=0.25, hi=6.0, budget=8)
        self.assertLess(abs(y - 33.7), 1.5)             # hit the target within tol
        self.assertGreater(len(trace), 2)

    def test_saturates_at_nearest_when_target_unreachable(self):
        # the one knob tops out at 40 < target 55 -> return the closest achievable (the high end),
        # honestly, instead of pretending it reached the target
        def evaluate(s):
            return min(40.0, 10.0 * s)
        scale, y, _ = P._bisect_scale(evaluate, target=55.0, lo=0.25, hi=6.0, budget=6)
        self.assertLessEqual(y, 40.0 + 1e-9)
        self.assertAlmostEqual(scale, 6.0)              # pinned at the max it could reach


class _FakeSB:
    """A mock SimBiology engine mirroring the real MTX wiring: the fittable knob is the Emax constant
    'Anti_CytSec_MaxbyMTX' (constant, feeds a rule); 'Anti_CytSec_MTX' is a rule OUTPUT (constant=
    False, must be ignored) and 'k12_MTX' is a PK constant (has no 'Max', must be ignored). A second-
    line dose would raise, so any read of the held-out arm during the fit fails the test loudly."""
    def __init__(self):
        self.meta = {"Anti_CytSec_MaxbyMTX": (0.5, True),       # Emax potency - the fit target
                     "Anti_CellInflux_MaxbyMTX": (0.5, True),   # Emax potency - the fit target
                     "Anti_CytSec_MTX": (0.0, False),           # rule output - NOT fittable
                     "k12_MTX": (0.18, True),                   # PK constant - NOT an efficacy knob
                     "kd_TNFa": (3.0, True)}                    # unrelated
        self.persisted = {}

    def list_parameters(self):
        return {"parameters": [{"name": n, "value": v, "constant": c}
                               for n, (v, c) in self.meta.items()]}

    def set_parameter(self, name, value):
        self.persisted[name] = value
        v, c = self.meta[name]
        self.meta[name] = (value, c)
        return value

    def run_vpop(self, xlsx, dose="", param_overrides="", **kw):
        assert "TCZ" not in dose, "held-out (second-line TCZ) arm read during the fit - LEAK"
        val = self.meta["Anti_CytSec_MaxbyMTX"][0]
        for kv in [x for x in param_overrides.split(";") if "=" in x]:
            n, v = kv.split("=", 1)
            if n == "Anti_CytSec_MaxbyMTX":
                val = float(v)
        k = int(round(val * 50))                        # ACR20% ~ Emax*50, as a flag column
        acr20 = [1.0] * k + [0.0] * (100 - k)
        return {"columns": {"ACR20": acr20}}


class TestFitClinical(unittest.TestCase):
    def _args_tasks(self):
        class A:
            limit, model = 100, "ra"
        tasks = {"run_columns": {"first_line": {"ACR20": "ACR20", "ACR50": "ACR50",
                                                "ACR70": "ACR70", "remission": "Rem"}},
                 "readout_states": None}
        return A(), tasks

    def test_selector_picks_emax_not_pk_or_rule_output(self):
        sb = _FakeSB()
        base = P._drug_potency_params(sb, "MTX")
        self.assertEqual(set(base), {"Anti_CytSec_MaxbyMTX", "Anti_CellInflux_MaxbyMTX"})
        self.assertNotIn("k12_MTX", base)          # PK constant excluded (no 'Max')
        self.assertNotIn("Anti_CytSec_MTX", base)  # rule output excluded (constant=False)

    def test_fits_mtx_to_first_line_target_and_persists(self):
        sb = _FakeSB()
        args, tasks = self._args_tasks()
        rc = tasks["run_columns"]
        r = P.fit_clinical(sb, "vpop.xlsx", args, tasks, b_day=200.0, r1=284.0,
                           mtx="MTX_15mg", rc=rc, target_acr20=33.7, budget=8, log=lambda *a: None)
        self.assertIsNotNone(r)
        self.assertLess(abs(r["acr20"] - 33.7), 2.0)                 # hit the training target
        self.assertIn("Anti_CytSec_MaxbyMTX", sb.persisted)         # Emax fit persisted into the model
        self.assertNotIn("k12_MTX", sb.persisted)                   # PK constant was NOT touched
        self.assertNotIn("TCZ", json.dumps(r))                      # no held-out leakage in the result

    def test_no_mtx_params_returns_none(self):
        sb = _FakeSB()
        sb.meta = {"kd_TNFa": (3.0, True)}                          # no MTX Emax potency constant
        args, tasks = self._args_tasks()
        r = P.fit_clinical(sb, "vpop.xlsx", args, tasks, 200.0, 284.0, "MTX", tasks["run_columns"],
                           33.7, log=lambda *a: None)
        self.assertIsNone(r)


if __name__ == "__main__":
    unittest.main()

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
    """A mock SimBiology engine: MTX PD params drive a monotone first-line ACR20; a second-line dose
    would raise, so any read of the held-out arm during the fit fails the test loudly."""
    def __init__(self):
        self.params = {"Anti_CytSec_MTX": 0.2, "Anti_CellInflux_MTX": 0.2, "kd_TNFa": 3.0}
        self.persisted = {}

    def list_parameters(self):
        return {"parameters": [{"name": n, "value": v, "constant": False}
                               for n, v in self.params.items()]}

    def set_parameter(self, name, value):
        self.persisted[name] = value
        self.params[name] = value
        return value

    def run_vpop(self, xlsx, dose="", param_overrides="", **kw):
        assert "TCZ" not in dose, "held-out (second-line TCZ) arm read during the fit - LEAK"
        val = self.params["Anti_CytSec_MTX"]
        for kv in [x for x in param_overrides.split(";") if "=" in x]:
            n, v = kv.split("=", 1)
            if n == "Anti_CytSec_MTX":
                val = float(v)
        k = int(round(val * 50))                        # ACR20% ~ val*50, as a flag column
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

    def test_fits_mtx_to_first_line_target_and_persists(self):
        sb = _FakeSB()
        args, tasks = self._args_tasks()
        rc = tasks["run_columns"]
        r = P.fit_clinical(sb, "vpop.xlsx", args, tasks, b_day=200.0, r1=284.0,
                           mtx="MTX_15mg", rc=rc, target_acr20=33.7, budget=8, log=lambda *a: None)
        self.assertIsNotNone(r)
        self.assertLess(abs(r["acr20"] - 33.7), 2.0)                 # hit the training target
        self.assertIn("Anti_CytSec_MTX", sb.persisted)              # fit was persisted into the model
        self.assertNotIn("TCZ", json.dumps(r))                      # no held-out leakage in the result

    def test_no_mtx_params_returns_none(self):
        sb = _FakeSB()
        sb.params = {"kd_TNFa": 3.0}                                # nothing ending in _MTX
        args, tasks = self._args_tasks()
        r = P.fit_clinical(sb, "vpop.xlsx", args, tasks, 200.0, 284.0, "MTX", tasks["run_columns"],
                           33.7, log=lambda *a: None)
        self.assertIsNone(r)


if __name__ == "__main__":
    unittest.main()

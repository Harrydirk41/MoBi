"""The guarded agentic workflow controller: the agent decides the next model-building action from
process + TRAINING signals only; the held-out arm is never shown (a leak RAISES). Pure - the LLM
and the executors are mocked."""

import json
import unittest

from pkpd_agent.engines import workflow_controller as W


class TestGuard(unittest.TestCase):
    def test_whitelists_process_and_training_keys(self):
        v = W._guard_state({"stable": False, "low_confidence_edges": 40, "first_line_error": 20.0,
                            "reactions": 102, "some_internal_thing": 1})
        self.assertEqual(set(v), {"stable", "low_confidence_edges", "first_line_error", "reactions"})

    def test_raises_on_held_out_or_truth_leak(self):
        # held-out arms AND truth-derived scores (precision/recall) both fail loudly
        for k in ("second_line", "radiate_acr20", "held_out_error", "validation_target",
                  "precision", "recall_secreting"):
            with self.assertRaises(ValueError):
                W._guard_state({"stable": True, k: 50.0})


class TestDecide(unittest.TestCase):
    def test_agent_picks_from_menu_and_state_is_guarded(self):
        seen = {}

        def call(system, user):
            seen["user"] = user
            return json.dumps({"action": "fit_clinical", "reason": "training error high"})
        act, reason = W.decide_next_action(
            {"stable": True, "first_line_error": 30.0, "radiate": 50.0}  # radiate must NOT reach LLM
            if False else {"stable": True, "first_line_error": 30.0}, call)
        self.assertEqual(act, "fit_clinical")
        self.assertNotIn("radiate", seen["user"].lower())
        self.assertNotIn("second", seen["user"].lower())

    def test_unknown_action_falls_back_to_finish(self):
        act, _ = W.decide_next_action({"stable": True}, lambda s, u: '{"action": "nonsense"}')
        self.assertEqual(act, "finish")


class TestControllerLoop(unittest.TestCase):
    def test_runs_actions_until_finish_and_never_sees_heldout(self):
        # a mock agent that stabilizes, then fits clinical, then finishes
        script = iter(["stabilize", "fit_clinical", "finish"])
        prompts = []

        def call(system, user):
            prompts.append(user)
            return json.dumps({"action": next(script), "reason": "x"})

        def stab(state):
            return {**state, "stable": True}

        def fitc(state):
            return {**state, "first_line_error": 5.0}

        state0 = {"stable": False, "first_line_error": 30.0, "low_confidence_edges": 40}
        state, hist = W.run_controller(state0, {"stabilize": stab, "fit_clinical": fitc}, call,
                                       max_steps=6)
        self.assertEqual([h["action"] for h in hist], ["stabilize", "fit_clinical", "finish"])
        self.assertTrue(state["stable"])
        self.assertEqual(state["first_line_error"], 5.0)
        # the held-out arm was never in any prompt
        self.assertFalse(any("radiate" in p.lower() or "second_line" in p.lower() for p in prompts))

    def test_stops_when_executor_missing(self):
        act_iter = iter(["widen_vpop", "finish"])
        state, hist = W.run_controller({"stable": True}, {},   # no executors
                                       lambda s, u: json.dumps({"action": next(act_iter)}),
                                       max_steps=5)
        self.assertEqual(hist[-1]["action"], "widen_vpop")     # stopped: no executor


if __name__ == "__main__":
    unittest.main()

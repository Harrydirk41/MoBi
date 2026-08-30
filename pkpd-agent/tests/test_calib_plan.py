"""Agent-driven L1: plan parsing + grading the agent's identifiability judgment. Stub LLM."""
import unittest
from pkpd_agent.engines import llm_calib_plan as CP


class TestCalibPlan(unittest.TestCase):
    def test_grade_catches_overclaim(self):
        call = lambda s, u: ('{"plan":[{"param":"kg","determinable":true,"reason":"scale"},'
                             '{"param":"K_IL6","determinable":true,"reason":"guessed"}]}')
        plan = CP.propose_plan(["Max"], ["kg", "K_IL6"], ["one steady state"], call)
        g = CP.grade_plan(plan, {"kg": True, "K_IL6": False})
        self.assertEqual(g["correct"], 1)                 # kg right, K_IL6 wrong
        self.assertEqual(g["overclaimed"], ["K_IL6"])     # claimed pinnable when it is not

    def test_perfect_judgment(self):
        call = lambda s, u: ('{"plan":[{"param":"kg","determinable":true},'
                             '{"param":"K_IL6","determinable":false}]}')
        plan = CP.propose_plan([], ["kg", "K_IL6"], [], call)
        g = CP.grade_plan(plan, {"kg": True, "K_IL6": False})
        self.assertEqual(g["accuracy"], 1.0)
        self.assertEqual(g["overclaimed"], [])


if __name__ == "__main__":
    unittest.main()


class TestOpenPlan(unittest.TestCase):
    def test_grade_open_rewards_naming_and_flagging_absence(self):
        call = lambda s, u: ('{"plan":[{"param":"kg","determinable":true,"needs":"steady state",'
                             '"needs_available":true},'
                             '{"param":"K6","determinable":false,"needs":"IL6 dose-response",'
                             '"needs_available":false}]}')
        plan = CP.propose_plan_open([{"name":"kg","has_literature_value":False},
                                     {"name":"K6","has_literature_value":False}], ["steady state"], call)
        g = CP.grade_open(plan, {"kg":True,"K6":False}, needs_kw=["dose"])
        self.assertEqual(g["id_accuracy"], 1.0)
        self.assertEqual(g["named_missing_experiment"], "1/1")   # named 'dose-response'
        self.assertEqual(g["flagged_data_absent"], "1/1")        # flagged it MISSING
        self.assertEqual(g["overclaimed"], [])

    def test_grade_open_catches_overclaim_on_shape(self):
        call = lambda s, u: ('{"plan":[{"param":"K6","determinable":true,'
                             '"needs":"none","needs_available":true}]}')
        plan = CP.propose_plan_open([{"name":"K6","has_literature_value":False}], ["x"], call)
        g = CP.grade_open(plan, {"K6":False}, needs_kw=["dose"])
        self.assertEqual(g["overclaimed"], ["K6"])

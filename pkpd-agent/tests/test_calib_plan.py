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

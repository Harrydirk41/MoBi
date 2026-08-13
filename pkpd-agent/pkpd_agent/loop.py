"""The decision loop: Observe -> Decide -> Act -> Evaluate, with a gate.

    session ─▶ policy.decide ─▶ tool calls ─▶ registry.dispatch
                    ▲                              │
                    │                              ▼
             observations  ◀── verification gates (act/evaluate results)

The loop is deliberately small and readable - it *is* the product. The
intelligence lives in the policy; the trust lives in the gates.
"""

from __future__ import annotations

from typing import Optional

from .config import AgentConfig
from .state import Decision, Finding, Finish, ModelingSession, Observation
from .system_prompt import SYSTEM_PROMPT
from .tools import ToolRegistry, build_default_registry
from .verification import run_gates


class DecisionLoop:
    def __init__(
        self,
        config: Optional[AgentConfig] = None,
        registry: Optional[ToolRegistry] = None,
        policy=None,
    ) -> None:
        self.config = config or AgentConfig()
        self.registry = registry or build_default_registry(self.config)
        self.policy = policy  # if None, an LLMPolicy is built per-run

    # ------------------------------------------------------------------ #
    def _build_policy(self):
        from .llm import LLMPolicy
        return LLMPolicy(self.config, self.registry, SYSTEM_PROMPT)

    def run(self, goal: str, session: Optional[ModelingSession] = None) -> ModelingSession:
        session = session or ModelingSession(goal=goal)
        policy = self.policy or self._build_policy()

        for _ in range(self.config.max_steps):
            if hasattr(policy, "observe"):
                policy.observe(session)

            step = policy.decide(session)

            # Finish?
            if step.__class__.__name__ == "FinishStep":
                session.record(Finish(step.text))
                break

            # Act: record the decision, then execute each call under the gate.
            calls = step.calls
            session.record(Decision(text=getattr(step, "text", ""), calls=calls))

            halted = False
            for call in calls:
                result = self.registry.dispatch(call.name, call.arguments, session)
                findings = self._verify(call.name, result)
                obs = Observation(
                    call_id=call.id,
                    tool=call.name,
                    ok=result.ok,
                    content=result.to_content(),
                    findings=findings,
                )
                session.record(obs)
                if obs.blocked and self.config.stop_on_block:
                    session.record(Finish(
                        "Halted: a verification BLOCK was raised and "
                        "stop_on_block is set.\n"
                        + "\n".join(f"[{f.level.upper()}] {f.gate}: {f.message}"
                                    for f in findings if f.level == "block")
                    ))
                    halted = True
                    break
            if halted:
                break
        else:
            session.record(Finish(f"Reached max_steps ({self.config.max_steps})."))

        return session

    # ------------------------------------------------------------------ #
    def _verify(self, tool_name: str, result) -> list[Finding]:
        # Only gate act/evaluate results (observations are inputs, not claims).
        if tool_name not in self.registry:
            return []
        phase = self.registry.get(tool_name).phase
        if phase == "observe":
            return []
        return run_gates(tool_name, result.to_content(), None)

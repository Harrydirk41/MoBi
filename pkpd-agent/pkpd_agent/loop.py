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

    def run(self, goal: str, session: Optional[ModelingSession] = None,
            on_event=None) -> ModelingSession:
        """Drive the loop. ``on_event(event)`` is called after each event is
        recorded, so callers can stream the trace live (each step includes a
        slow tool call, so batched output would look frozen)."""
        session = session or ModelingSession(goal=goal)
        policy = self.policy or self._build_policy()

        def emit(ev):
            session.record(ev)
            if on_event is not None:
                on_event(ev)

        for _ in range(self.config.max_steps):
            if hasattr(policy, "observe"):
                policy.observe(session)

            step = policy.decide(session)

            # Finish?
            if step.__class__.__name__ == "FinishStep":
                emit(Finish(step.text))
                break

            # Act: record the decision (emit before the slow dispatch), then
            # execute each call under the gate.
            calls = step.calls
            emit(Decision(text=getattr(step, "text", ""), calls=calls))

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
                emit(obs)
                if obs.blocked and self.config.stop_on_block:
                    emit(Finish(
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
            emit(Finish(f"Reached max_steps ({self.config.max_steps})."))

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

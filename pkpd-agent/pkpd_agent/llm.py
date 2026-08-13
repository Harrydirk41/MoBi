"""The decision-maker ("brain"), behind a swappable Policy interface.

The loop is engine- and model-agnostic: it asks a ``Policy`` what to do next
given the transcript, and the policy answers with either tool calls (Act) or a
final message (Finish).

  * ``LLMPolicy``      - Claude drives the decisions (real runs).
  * ``ScriptedPolicy`` - a fixed list of steps (tests / dry runs, no API key).

This split is the honest architecture: the LLM is the judgment layer and is
cleanly replaceable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Protocol

from .state import ToolCall


# --------------------------------------------------------------------------- #
# What a policy returns
# --------------------------------------------------------------------------- #

@dataclass
class ActStep:
    """The policy chose to act: reasoning text + one or more tool calls."""
    text: str
    calls: list[ToolCall]


@dataclass
class FinishStep:
    text: str


PolicyStep = ActStep | FinishStep


class Policy(Protocol):
    def decide(self, session) -> PolicyStep:  # noqa: D401
        """Given the session/transcript so far, choose the next step."""
        ...


# --------------------------------------------------------------------------- #
# Scripted policy (deterministic, no LLM)
# --------------------------------------------------------------------------- #

@dataclass
class ScriptedPolicy:
    """Replays a fixed list of steps. Each entry is one of:

        ("finish", "text")                          -> stop with a message
        ("call", "tool_name", {args})               -> one tool call
        ("calls", [("tool", {args}), ...])          -> parallel tool calls
    """

    steps: list[Any]
    _i: int = field(default=0, init=False)

    def decide(self, session) -> PolicyStep:
        if self._i >= len(self.steps):
            return FinishStep("scripted policy exhausted")
        step = self.steps[self._i]
        self._i += 1
        kind = step[0]
        if kind == "finish":
            return FinishStep(step[1])
        if kind == "call":
            raw_calls = [(step[1], step[2])]
        elif kind == "calls":
            raw_calls = list(step[1])
        else:
            raise ValueError(f"unknown scripted step kind: {kind!r}")
        calls = [
            ToolCall(id=f"call_{self._i}_{j}", name=name, arguments=args)
            for j, (name, args) in enumerate(raw_calls)
        ]
        return ActStep(text="", calls=calls)


# --------------------------------------------------------------------------- #
# LLM policy (Claude, manual tool-use loop)
# --------------------------------------------------------------------------- #

class LLMPolicy:
    """Wraps the Anthropic Messages API. Maintains its own message history so
    thinking blocks and tool_use/tool_result pairing are preserved across the
    manual loop."""

    def __init__(self, config, registry, system_prompt: str) -> None:
        self.config = config
        self.registry = registry
        self.system_prompt = system_prompt
        self._messages: list[dict[str, Any]] = []
        self._client = None  # lazy

    def _ensure_client(self):
        if self._client is None:
            import anthropic  # imported lazily so the package works without it
            self._client = anthropic.Anthropic()
        return self._client

    def observe(self, session) -> None:
        """Feed the newest observations (tool results) back as a user turn.

        The loop calls this after executing tool calls, before the next
        ``decide``. On the very first call it seeds the goal.
        """
        from .state import Observation

        if not self._messages:
            self._messages.append({
                "role": "user",
                "content": f"Modeling goal:\n{session.goal}\n\nBegin.",
            })
            return

        # Gather the observations produced since the last decision, as
        # tool_result blocks (with verification findings appended).
        pending: list[dict[str, Any]] = []
        for ev in reversed(session.transcript):
            if isinstance(ev, Observation):
                findings_txt = ""
                if ev.findings:
                    findings_txt = "\n\nVERIFICATION:\n" + "\n".join(
                        f"[{f.level.upper()}] {f.gate}: {f.message}" for f in ev.findings
                    )
                pending.append({
                    "type": "tool_result",
                    "tool_use_id": ev.call_id,
                    "content": _json(ev.content) + findings_txt,
                    "is_error": not ev.ok,
                })
            else:
                break
        pending.reverse()
        if pending:
            self._messages.append({"role": "user", "content": pending})

    def decide(self, session) -> PolicyStep:
        client = self._ensure_client()
        resp = client.messages.create(
            model=self.config.model,
            max_tokens=self.config.max_tokens,
            system=self.system_prompt,
            tools=self.registry.to_anthropic_schema(),
            thinking={"type": "adaptive"},
            output_config={"effort": self.config.effort},
            messages=self._messages,
        )
        # Preserve the full assistant turn (incl. thinking) for the next request.
        self._messages.append({"role": "assistant", "content": resp.content})

        if getattr(resp, "stop_reason", None) == "refusal":
            return FinishStep("The model declined this request (safety refusal).")

        text_parts, calls = [], []
        for block in resp.content:
            if block.type == "text":
                text_parts.append(block.text)
            elif block.type == "tool_use":
                calls.append(ToolCall(id=block.id, name=block.name,
                                      arguments=dict(block.input)))
        text = "\n".join(text_parts).strip()
        if calls:
            return ActStep(text=text, calls=calls)
        return FinishStep(text or "(no further action)")


def _json(obj: Any) -> str:
    import json
    return json.dumps(obj, default=str, indent=2)

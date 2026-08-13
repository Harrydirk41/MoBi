"""Session state and the event transcript that flows through the loop.

The transcript is the single source of truth. Both the LLM policy and the
scripted test policy consume it, and it doubles as the provenance record for a
modeling session (every decision, every result, every verification verdict).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Literal


# --------------------------------------------------------------------------- #
# Transcript events
# --------------------------------------------------------------------------- #

@dataclass
class Goal:
    """The scientific objective the loop is optimizing."""
    text: str


@dataclass
class ToolCall:
    id: str
    name: str
    arguments: dict[str, Any]


@dataclass
class Decision:
    """One Decide step: the model's reasoning plus the tool calls it chose."""
    text: str
    calls: list[ToolCall] = field(default_factory=list)


@dataclass
class Observation:
    """Result of executing one tool, plus any verification findings."""
    call_id: str
    tool: str
    ok: bool
    content: dict[str, Any]
    findings: list["Finding"] = field(default_factory=list)

    @property
    def blocked(self) -> bool:
        return any(f.level == "block" for f in self.findings)


@dataclass
class Finding:
    """A verification verdict on an action's result."""
    level: Literal["info", "warn", "block"]
    gate: str
    message: str


@dataclass
class Finish:
    text: str


Event = Goal | Decision | Observation | Finish


# --------------------------------------------------------------------------- #
# Session
# --------------------------------------------------------------------------- #

@dataclass
class ModelingSession:
    """Holds the transcript and a small key/value store of named artifacts
    (loaded models, fit results, simulation handles) that tools produce and
    later tools consume."""

    goal: str
    transcript: list[Event] = field(default_factory=list)
    artifacts: dict[str, Any] = field(default_factory=dict)
    started_at: float = field(default_factory=time.time)

    def __post_init__(self) -> None:
        if not any(isinstance(e, Goal) for e in self.transcript):
            self.transcript.insert(0, Goal(self.goal))

    # -- recording ------------------------------------------------------- #
    def record(self, event: Event) -> None:
        self.transcript.append(event)

    def put(self, name: str, value: Any) -> None:
        self.artifacts[name] = value

    def get(self, name: str, default: Any = None) -> Any:
        return self.artifacts.get(name, default)

    # -- queries --------------------------------------------------------- #
    @property
    def observations(self) -> list[Observation]:
        return [e for e in self.transcript if isinstance(e, Observation)]

    @property
    def finished(self) -> bool:
        return any(isinstance(e, Finish) for e in self.transcript)

    def summary(self) -> str:
        n_dec = sum(isinstance(e, Decision) for e in self.transcript)
        n_obs = len(self.observations)
        n_block = sum(o.blocked for o in self.observations)
        return (
            f"goal={self.goal!r} decisions={n_dec} actions={n_obs} "
            f"blocked={n_block} finished={self.finished}"
        )

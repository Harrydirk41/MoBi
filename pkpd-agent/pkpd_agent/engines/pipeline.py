"""Decomposable build pipeline: every modelling step is a pluggable Provider whose MODE says how
it is filled - human-given, a direct LLM read, an LLM+tool call, or a deterministic data/config
lookup. The pipeline runs the providers in dependency order and writes each layer into the spec,
tagging its source by the provider's mode.

This is the "each step, either LLM(+tool) or human/read" design: swap one layer's provider without
touching the others. No step is hardcoded - a baked-in default is just a `given(...)` provider, and
reading the model's own structure is just a `data(...)` provider; both are explicit choices, not
scaffolding.

    providers = {
        "frame":    given({...}) | from_llm(...) ,      # human tier / agent proposes
        "target":   given("IL6") | from_llm(...) ,      # which node to build
        "topology": from_llm(reg_prompt) | given([...]) | data(all_candidates),
        "form":     from_llm(motif_prompt) | given(default_motif),
    }
    spec = run(providers, ctx)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from pkpd_agent.engines import model_spec as MS

GIVEN, LLM, LLM_TOOL, DATA = "given", "llm", "llm+tool", "data"


@dataclass
class Provider:
    """One step. ``mode`` records HOW it is filled (for the source tag); ``fn(ctx)`` returns the
    layer's value. ``fn`` may use ctx['call'] (the LLM boundary), ctx tools, or config data."""
    mode: str
    fn: Callable[[dict], Any]


# ── constructors: the four ways to fill a step ──
def given(value: Any) -> Provider:
    """Human-done: the value is supplied directly."""
    return Provider(GIVEN, lambda ctx: value)


def from_llm(system: str, user: "str | Callable[[dict], str]", parse: Callable[[str], Any]) -> Provider:
    """Direct LLM read: call(system,user) -> parsed value. No tools."""
    def run(ctx):
        u = user(ctx) if callable(user) else user
        return parse(ctx["call"](system, u))
    return Provider(LLM, run)


def from_llm_tool(fn: Callable[[dict], Any]) -> Provider:
    """LLM + tool: fn drives an LLM that may call tools (web, data lookup, retrieval)."""
    return Provider(LLM_TOOL, fn)


def data(fn: Callable[[dict], Any]) -> Provider:
    """Deterministic config/model lookup (no LLM) - e.g. read the model's own structure."""
    return Provider(DATA, fn)


# ── the pipeline ──
_LAYERS = ("frame", "target", "topology", "form")


def run(providers: dict, ctx: dict) -> dict:
    """Fill the spec layer by layer. Each provider writes its layer; its mode is recorded so the
    spec says, per layer, whether a human, a direct LLM read, an LLM+tool, or data produced it."""
    missing = [l for l in _LAYERS if l not in providers]
    if missing:
        raise ValueError(f"pipeline missing providers for: {missing}")
    modes = {}
    frame = providers["frame"].fn(ctx); modes["frame"] = providers["frame"].mode; ctx["frame"] = frame
    target = providers["target"].fn(ctx); modes["target"] = providers["target"].mode
    ctx["target"] = target
    topology = providers["topology"].fn(ctx); modes["topology"] = providers["topology"].mode
    form = providers["form"].fn(ctx); modes["form"] = providers["form"].mode

    spec = MS.build_spec(frame, target, topology, form, ctx["prov"], ctx["levels"],
                         truth_regulators=ctx.get("truth"))
    spec["modes"] = modes                                  # per-layer: given|llm|llm+tool|data
    return spec

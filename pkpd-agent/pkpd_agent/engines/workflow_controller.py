r"""A GUARDED agentic workflow controller for QSP model-building.

The design flaw the fixed pipeline had: a human-written script decided the order of steps
(build -> calibrate -> vpop -> predict), so the agent never DECIDED to, say, calibrate to the
clinical training data - it just wasn't in the script. This controller closes that: at each step the
agent sees the current PROCESS signals (is the network stable? how much over-inclusion? calibration
drift?) and the TRAINING fit error (first-line, the calibration arm), and CHOOSES the next action to
improve the model - until it decides the model is good enough to predict the held-out test.

THE RED LINE (enforced, not just asked): the controller may show the agent process signals and the
TRAINING (first-line) error, but NEVER the held-out validation data (the second-line RADIATE trial).
Deciding or fitting from held-out would be cheating. ``_guard_state`` strips the state shown to the
agent to a whitelist and RAISES if any held-out-looking key is present, so a leak fails loudly
instead of silently training on the answer.

Actions are injected as ``executors`` (name -> fn(state) -> new_state), so the same controller runs
with real MATLAB executors in the pipeline and with mocks in tests.
"""

from __future__ import annotations

import json

from .llm_structure import _parse_json

# the bounded action menu the agent chooses from (descriptions are process/training framed only)
CONTROLLER_ACTIONS = {
    "stabilize": "remove the agent's own edges that make the coupled ODEs diverge (dynamics signal)",
    "prune": "drop your low-confidence, uncited edges to cut over-inclusion",
    "force_influx": "give a marginal cell a synthesized influx arm (fill a missing-data gap so an "
                    "influx-suppressing therapy can act)",
    "fit_clinical": "calibrate the free rates to the FIRST-LINE training response (the calibration "
                    "arm) - legitimate training, never the held-out arm",
    "widen_vpop": "widen the virtual-population sampling span so the baseline severity spreads",
    "finish": "stop refining and predict the held-out validation test",
}

# Only these keys may reach the agent. They are all PROCESS signals (from the dynamics) or the
# agent's OWN self-knowledge (its confidence) or a TRAINING error - none is derived from the model's
# truth. precision/recall are deliberately EXCLUDED: they are computed against the model's true edge
# set, so showing them would leak "how close your structure is to the answer" (a real, if mild, form
# of cheating). The agent decides over-inclusion from its own low-confidence-edge count instead.
_ALLOWED_KEYS = {"stable", "calibration_drift", "marginal_cells", "low_confidence_edges",
                 "first_line_error", "baseline_offset", "steps_taken", "last_action",
                 "last_action_effect", "diverged_species", "reactions", "actions_done"}
_FORBIDDEN_SUBSTR = ("second_line", "second-line", "radiate", "held_out", "heldout", "held-out",
                     "validation", "precision", "recall", "truth")

_CONTROLLER_SYS = (
    "You are the CONTROLLER of a QSP model-building workflow. At each step you are given PROCESS "
    "signals (network stability, which species diverge, calibration drift), your OWN over-inclusion "
    "signal (how many of YOUR edges you marked low-confidence), which cells are marginal (have no "
    "influx arm), and 'last_action_effect' telling you whether your previous action actually changed "
    "the model. Choose ONE next action. You are NEVER shown - and must never try to fit - any "
    "held-out validation data or any measure of how close your structure is to a reference (no "
    "precision/recall); decide only from the dynamics, your own confidence, and the training error. "
    "Address the worst signal first: fix an UNSTABLE network before anything; then reduce your "
    "low-confidence edges; give a marginal cell an influx arm if a therapy needs it. Do NOT repeat "
    "an action whose last_action_effect was 'no-op (saturated)' - move on. Choose 'finish' when the "
    "network is stable and few low-confidence edges remain. Return JSON "
    '{"action": one of the listed names, "reason": "one phrase"}.')


def _guard_state(state):
    """Whitelist the state shown to the agent; RAISE if any held-out-looking key is present so a
    leak fails loudly rather than silently letting the agent train on the answer."""
    bad = [k for k in state if any(s in k.lower() for s in _FORBIDDEN_SUBSTR)]
    if bad:
        raise ValueError(f"held-out data leaked into the controller state: {bad}")
    return {k: v for k, v in state.items() if k in _ALLOWED_KEYS}


def decide_next_action(state, call, actions=CONTROLLER_ACTIONS):
    """Ask the agent for the next workflow action given the guarded (process + training) state.
    Returns (action, reason); an unrecognized action falls back to 'finish'."""
    view = _guard_state(state)
    menu = "\n".join(f"  - {name}: {desc}" for name, desc in actions.items())
    user = (f"Current model signals (process + TRAINING only):\n{json.dumps(view, default=str)}\n\n"
            f"Available actions:\n{menu}\n\nChoose the next action.")
    d = _parse_json(call(_CONTROLLER_SYS, user))
    act = d.get("action")
    return (act if act in actions else "finish"), d.get("reason")


def run_controller(state, executors, call, max_steps=6, log=lambda *a: None,
                   actions=CONTROLLER_ACTIONS):
    """The guarded decision loop: the agent picks an action from process+training signals, the
    matching executor runs it and returns the updated (process+training) state, repeat until the
    agent chooses 'finish' or an action has no executor. ``actions`` restricts the menu the agent
    sees (e.g. only the clinical-stage actions). Held-out is NOT evaluated here - the caller runs it
    once, after, on the returned model. Returns (state, history)."""
    history, done = [], []
    for step in range(max_steps):
        act, reason = decide_next_action({**state, "steps_taken": step, "actions_done": list(done)},
                                         call, actions=actions)
        history.append({"step": step, "action": act, "reason": reason})
        log(f"  controller step {step}: {act}  ({reason})")
        if act == "finish":
            break
        ex = executors.get(act)
        if ex is None:
            log(f"  (no executor for '{act}'; stopping)")
            break
        state = ex(state)                                  # executor returns process+training state
        done.append(act)
    return state, history

"""Real run: Claude drives the loop over the (still mocked) engines.

Requires ANTHROPIC_API_KEY and the `anthropic` package. Engines stay in mock
mode by default so you can watch the *decision-making* without a pharmpy /
OSP / NONMEM install; set mock=False once the engines are wired to your setup.

    pip install anthropic
    export ANTHROPIC_API_KEY=...
    python -m examples.run_llm "Fit a one-compartment popPK model to warfarin.mod"
"""

import os
import sys

from pkpd_agent.config import AgentConfig
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, Observation


def main() -> None:
    goal = sys.argv[1] if len(sys.argv) > 1 else (
        "Load the builtin dataset. Fit it with the fast pkfit engine and with "
        "nlmixr2 (true NLME); compare the residual-error estimates and tell me "
        "which to trust and why. Qualify the chosen model with a VPC."
    )
    # REAL engines: pkfit always real; nlmixr2/OSP become real when an Rscript
    # with the backends is provided via the PKPD_RSCRIPT environment variable.
    cfg = AgentConfig(mock=False,
                      rscript_path=os.environ.get("PKPD_RSCRIPT", "Rscript"))
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set - see demo_dry_run.py for a no-key run.")
        return

    loop = DecisionLoop(config=cfg)  # policy=None -> LLMPolicy is built
    session = loop.run(goal)

    for ev in session.transcript:
        if isinstance(ev, Decision) and ev.text:
            print(f"\n[reason] {ev.text}")
        if isinstance(ev, Decision):
            for c in ev.calls:
                print(f"  -> {c.name} {c.arguments}")
        elif isinstance(ev, Observation):
            print(f"  <- {ev.tool}: {ev.content.get('message','')}")
            for f in ev.findings:
                print(f"     [{f.level.upper()}] {f.message}")
        elif isinstance(ev, Finish):
            print(f"\n=== SUMMARY ===\n{ev.text}")


if __name__ == "__main__":
    main()

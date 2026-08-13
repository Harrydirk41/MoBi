"""Runnable dry-run: the full Observe -> Decide -> Act -> Evaluate loop with a
scripted 'brain' and mocked engines. No API key, no pharmpy, no OSP required.

    python -m examples.demo_dry_run     (from the pkpd-agent/ directory)
"""

from pkpd_agent.config import AgentConfig
from pkpd_agent.llm import ScriptedPolicy
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.state import Decision, Finish, Observation


def main() -> None:
    # A scripted decision sequence standing in for the LLM:
    #   observe (NCA) -> observe (load model) -> act (fit) -> evaluate (VPC) -> finish
    policy = ScriptedPolicy(steps=[
        ("call", "nca_analyze",
         {"times": [0, 0.5, 1, 2, 4, 8, 12, 24],
          "concentrations": [0, 3.1, 5.0, 4.2, 2.8, 1.3, 0.6, 0.15]}),
        ("call", "pharmpy_load_model", {"path": "warfarin.mod"}),
        ("call", "pharmpy_fit", {"model_id": "model::warfarin.mod"}),
        ("call", "pharmpy_vpc", {"model_id": "model::warfarin.mod"}),
        ("finish",
         "NCA gave Cmax~5, AUC finite; fitted a 1-cpt oral model (minimization "
         "OK, RSEs reasonable); VPC coverage acceptable. A human should confirm "
         "the covariate model and check GOF plots before use."),
    ])

    loop = DecisionLoop(config=AgentConfig(mock=True), policy=policy)
    session = loop.run("Characterize warfarin PK and fit a population model.")

    print("=" * 70)
    print("GOAL:", session.goal)
    print("=" * 70)
    for ev in session.transcript:
        if isinstance(ev, Decision):
            for c in ev.calls:
                print(f"  DECIDE -> {c.name}({_short(c.arguments)})")
        elif isinstance(ev, Observation):
            tag = "OK " if ev.ok else "ERR"
            print(f"  ACT   [{tag}] {ev.tool}: {ev.content.get('message','')}")
            for f in ev.findings:
                print(f"          [{f.level.upper()}] {f.gate}: {f.message}")
        elif isinstance(ev, Finish):
            print("-" * 70)
            print("FINISH:", ev.text)
    print("=" * 70)
    print(session.summary())


def _short(d: dict) -> str:
    s = ", ".join(f"{k}={v}" for k, v in d.items())
    return s if len(s) < 60 else s[:57] + "..."


if __name__ == "__main__":
    main()

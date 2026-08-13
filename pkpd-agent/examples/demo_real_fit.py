"""Real end-to-end run: the decision policy drives the REAL pkfit engine.

Everything here actually computes - real maximum-likelihood fits, real AIC /
likelihood-ratio model selection, a real Monte-Carlo VPC. No API key, no
external backend. Requires numpy + scipy.

    python -m examples.demo_real_fit
"""

from pkpd_agent.config import AgentConfig
from pkpd_agent.loop import DecisionLoop
from pkpd_agent.policies import PharmacometricPolicy
from pkpd_agent.state import Decision, Finish, Observation


def main() -> None:
    policy = PharmacometricPolicy(covariate_param="CL", covariate="WT", covariate_ref=70)
    loop = DecisionLoop(config=AgentConfig(mock=False), policy=policy)
    session = loop.run("Build and qualify a population PK model for this dataset.")

    print("=" * 74)
    print("GOAL:", session.goal)
    print("=" * 74)
    for ev in session.transcript:
        if isinstance(ev, Decision):
            if ev.text:
                print("\n" + ev.text)
            for c in ev.calls:
                print(f"   -> {c.name}({_short(c.arguments)})")
        elif isinstance(ev, Observation):
            m = ev.content
            extra = ""
            if ev.tool == "pkfit_fit":
                extra = (f"  OFV={m.get('ofv')} AIC={m.get('aic')} "
                         f"ok={m.get('minimization_successful')} "
                         f"est={m.get('parameter_estimates')}")
            elif ev.tool == "pkfit_nca":
                extra = f"  Cmax={m.get('c_max')} Tmax={m.get('t_max')} AUC={m.get('auc')}"
            elif ev.tool == "pkfit_vpc":
                extra = f"  {m.get('pct_observations_within_90_pi')}% within 90% PI"
            print(f"   <- {ev.tool}: {m.get('message','')}{extra}")
            for f in ev.findings:
                print(f"      [{f.level.upper()}] {f.gate}: {f.message}")
        elif isinstance(ev, Finish):
            print("\n" + "-" * 74)
            print(ev.text)
    print("=" * 74)
    print(session.summary())


def _short(d: dict) -> str:
    return ", ".join(f"{k}={v}" for k, v in d.items())


if __name__ == "__main__":
    main()

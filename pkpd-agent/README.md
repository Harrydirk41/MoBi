# pkpd-agent

An **LLM decision loop** over pharmacometric engines. The model makes the
modeling *decisions* a human pharmacometrician would make; trusted engines do
the heavy lifting; a verification layer checks each decision for scientific
sanity before it is accepted.

This is the "executable loop of decision" idea, made concrete:

```
   ┌─────────────────────────────────────────────────────┐
   │ OBSERVE   load model / snapshot, run NCA             │
   ├─────────────────────────────────────────────────────┤
   │ DECIDE    the LLM picks ONE tool call + reasoning    │
   ├─────────────────────────────────────────────────────┤
   │ ACT       run it via pharmpy / OSP / NCA             │
   ├─────────────────────────────────────────────────────┤
   │ EVALUATE  verification gates check the result;       │
   │           a [BLOCK] is fed back for the LLM to fix   │
   └──────────────┬──────────────────────────────────────┘
                  │  not done → back to OBSERVE
                  ▼  done → plain-language summary
```

## Why this shape

The engines guarantee **mechanical** correctness (the fit ran, the ODE
integrated). They do **not** guarantee **scientific** correctness (is the fit a
junk local minimum? is the PBPK model mass-balanced?). So the design is:

> **LLM decides → the loop makes it *prove* the decision.**

The intelligence lives in the policy; the trust lives in the gates.

## Engines (the action space)

| Engine | Role | Status | Tools |
|---|---|---|---|
| **pkfit** | real MLE PK fitting + NCA + Monte-Carlo VPC | **real (runs here)** | `pkfit_load_data`, `pkfit_nca`, `pkfit_fit`, `pkfit_vpc` |
| **nlmixr2** | true NLME popPK (random effects, SAEM/FOCEi) | real via R backend | `nlmixr2_fit` |
| **pharmpy** | population PK/PD estimation, AMD, VPC | real via Python + backend | `pharmpy_load_model`, `pharmpy_fit`, `pharmpy_run_amd`, `pharmpy_vpc` |
| **OSP** (MoBi / PK-Sim) | mechanistic PBPK / QSP simulation | real via R (`ospsuite`) | `osp_load_snapshot`, `osp_set_parameter`, `osp_simulate` |
| **NCA** | model-free first pass (gap-filling binding) | real (builtin) | `nca_analyze` |

**Wiring the real backends:** see **[SETUP_WINDOWS.md](SETUP_WINDOWS.md)** for a
step-by-step (nlmixr2 + pharmpy + OSP from a conda env). Check what's live with:

```bash
python -m pkpd_agent.doctor                       # what's installed / what it unlocks
python -m pkpd_agent.doctor --rscript "C:\...\Rscript.exe"
```

The real R/.NET backends (nlmixr2, ospsuite) are reached through small worker
scripts in `engines/r_workers/`; each tool errors *clearly* (never crashes) when
its backend is absent, so the loop is always runnable.

### The real engine: `pkfit`

`pkfit` (`engines/pkfit.py`, numpy/scipy) is a genuine estimator that runs
in-process with **no NONMEM / nlmixr2 / R**:

- closed-form 1- and 2-compartment oral PK models,
- **maximum-likelihood** fitting (naive-pooled, proportional error),
- **standard errors** from the observed information matrix, condition number,
- real **model comparison** (OFV / AIC / BIC / likelihood-ratio test),
- optional covariate model (`param = param_pop * (cov/ref) ** coef`, `coef`
  estimated so a 1-df LRT is meaningful),
- Monte-Carlo **VPC**.

It ships a builtin dataset **simulated from a known truth** (1-cpt oral,
allometric WT on CL), so you can watch the fit *recover the truth*.

> Honest scope: `pkfit` is **naive-pooled**, not full nonlinear mixed-effects.
> Full NLME with random effects is exactly the job that needs an external
> backend (the pharmpy/NONMEM world). `pkfit` gives the loop a real estimator to
> drive end to end, with real convergence, precision, and model selection.

The other engines keep a **mock mode** (default) so the whole loop still runs
with no pharmpy / OSP install; their real call paths are written and guarded
behind `config.mock`.

## Run it

```bash
cd pkpd-agent
pip install numpy scipy                 # for the real engine
python -m unittest discover -s tests    # 20 tests (5 exercise the real fits)
python -m examples.demo_real_fit        # REAL fits + real model selection, no API key
python -m examples.demo_dry_run         # mock engines, scripted brain (stdlib only)
```

`demo_real_fit` runs the `PharmacometricPolicy` (a transparent, non-LLM expert
system) over the real engine. It fits a base 1-cpt model, tries 2-cpt, **rejects
it** (the fit is degenerate — the verification gate raises a `[BLOCK]`), **keeps
WT-on-CL** by a likelihood-ratio test (ΔOFV ≈ 32 ≫ 3.84), qualifies with a VPC,
and validates the estimates against the known truth. That is "good decision
making" you can run without an API key.

## Run it with Claude driving (real decisions)

```bash
pip install anthropic
export ANTHROPIC_API_KEY=...
python -m examples.run_llm "Load warfarin.mod, fit it, run a VPC, and tell me whether the fit is trustworthy."
```

Engines stay mocked by default so you can watch the *decision-making* before
wiring the real backends.

## Architecture

```
pkpd_agent/
  config.py            AgentConfig (mock switch, model, effort, gating policy)
  state.py             ModelingSession + the event transcript (also the provenance record)
  system_prompt.py     the 'modeler brain' instructions
  loop.py              DecisionLoop — the small, readable Observe→Decide→Act→Evaluate driver
  llm.py               Policy interface: LLMPolicy (Claude) | ScriptedPolicy (tests)
  policies.py          PharmacometricPolicy — a transparent, non-LLM decision policy
  engines/
    pkfit.py           the REAL estimation engine (numpy/scipy MLE, NCA, VPC)
  tools/
    registry.py        Tool schemas, dispatch, phase tagging (observe/act/evaluate)
    pkfit_tools.py     tools for the real engine
    pharmpy_tools.py   pharmpy adapter + tools (mock / backend)
    osp_tools.py       OSP (snapshot JSON / Rscript / MoBi.CLI) adapter + tools
    nca_tools.py       generic NCA binding
  verification/
    gates.py           scientific sanity checks (convergence, RSE, mass balance, VPC coverage)
```

The **brain is swappable** (`Policy`): `LLMPolicy` (Claude) for production,
`PharmacometricPolicy` (expert-system statistics) and `ScriptedPolicy` (tests)
for no-key runs — all interchangeable over the same loop, tools, and gates.
That is the honest statement that the LLM is the judgment layer and is cleanly
replaceable.

## Driving the real engine with Claude

```bash
pip install anthropic numpy scipy
export ANTHROPIC_API_KEY=...
python -m examples.run_llm "Load the builtin dataset, build a popPK model \
(compare 1- vs 2-compartment by AIC, test WT on CL by LRT), qualify it with a \
VPC, and tell me whether it is trustworthy."
```

Claude makes the same decisions the expert policy does, but reasons in natural
language and can depart from the fixed workflow — while every fit still passes
through the verification gates.

## Wiring the real engines

- **pharmpy** — `pip install pharmpy-core`, set `config.mock = False`. Needs an
  estimation backend (NONMEM / nlmixr2) for `pharmpy_fit` / `pharmpy_run_amd`.
- **OSP** — set `config.mock = False` and provide either `rscript_path` (with the
  `ospsuite` R package) or `mobi_cli_path` (MoBi.CLI). Runs on the Windows/.NET
  side; the LLM orchestrator can live anywhere and exchange **snapshot JSON**.

## Scope & honesty

The estimation engine (`pkfit`) is **real and runs here**. The population-NLME
(pharmpy) and mechanistic-ODE (OSP) engines are **mocked**, with real call paths
written and guarded — those are the worlds that need an external backend
(NONMEM/nlmixr2) or the Windows/.NET runtime. It does **not** yet cover MCP-Mod,
MBMA, optimal design, or PBBM — additional bindings, by design.

What it demonstrates is the load-bearing idea, end to end on real numbers: **a
decision-maker (LLM or expert policy) building a pharmacometric model, with
every fit verified before it counts** — recover the truth, reject the
over-parameterized model, keep the significant covariate, qualify with a VPC.

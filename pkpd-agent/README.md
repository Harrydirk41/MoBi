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

| Engine | Role | Tools |
|---|---|---|
| **pharmpy** | population PK/PD (NLME) estimation, AMD, VPC | `pharmpy_load_model`, `pharmpy_fit`, `pharmpy_run_amd`, `pharmpy_vpc` |
| **OSP** (MoBi / PK-Sim) | mechanistic PBPK / QSP simulation | `osp_load_snapshot`, `osp_set_parameter`, `osp_simulate` |
| **NCA** | model-free first pass (gap-filling binding) | `nca_analyze` |

Each engine has a **mock mode** (default) so the whole loop runs with no
pharmpy, no OSP/R install, and no API key. The real call paths are written
against the actual engine APIs and guarded behind `config.mock`.

## Run it (no dependencies, no API key)

```bash
cd pkpd-agent
python -m unittest discover -s tests   # 15 tests, stdlib only
python -m examples.demo_dry_run        # full loop, scripted 'brain', mock engines
```

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
  tools/
    registry.py        Tool schemas, dispatch, phase tagging (observe/act/evaluate)
    pharmpy_tools.py   pharmpy adapter + tools
    osp_tools.py       OSP (snapshot JSON / Rscript / MoBi.CLI) adapter + tools
    nca_tools.py       NCA binding (the gap pharmpy & OSP don't cover)
  verification/
    gates.py           scientific sanity checks (convergence, RSE, mass balance, VPC coverage)
```

The **brain is swappable** (`Policy`): `LLMPolicy` for real runs, `ScriptedPolicy`
for deterministic tests — the honest statement that the LLM is the judgment
layer and is cleanly replaceable.

## Wiring the real engines

- **pharmpy** — `pip install pharmpy-core`, set `config.mock = False`. Needs an
  estimation backend (NONMEM / nlmixr2) for `pharmpy_fit` / `pharmpy_run_amd`.
- **OSP** — set `config.mock = False` and provide either `rscript_path` (with the
  `ospsuite` R package) or `mobi_cli_path` (MoBi.CLI). Runs on the Windows/.NET
  side; the LLM orchestrator can live anywhere and exchange **snapshot JSON**.

## Scope & honesty

This is a **v0 skeleton**: real architecture, mocked engines. It covers the two
core engine worlds (NLME estimation + mechanistic ODE) plus one gap binding
(NCA). It does **not** yet cover MCP-Mod, MBMA, optimal design, or PBBM — those
are additional bindings, by design. What it does demonstrate is the load-bearing
idea: **an LLM making pharmacometric decisions, each one verified before it
counts.**

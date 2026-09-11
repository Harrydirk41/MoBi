# examples/ — index

Most of these are research/build scripts (the lab notebook). **A coding agent should not
wade through them** — the clean agent surface is the `pkpd-bench` CLI (see
`../.claude/skills/pkpd-modeling/SKILL.md` and `../SANDBOX_CC.md`). This index is for humans
navigating the pile.

> Run any script as `python -m examples.<name> ...`. Scripts marked ⚙️ need PK-Sim
> (`PKPD_PKSIM_CLI`); 🧪 need MATLAB + SimBiology; 🔑 need `ANTHROPIC_API_KEY`.

## ▶ Agent entry points (what you actually drive)

| Command | What it does |
| --- | --- |
| **`pkpd-bench …`** (not here — installed CLI) ⚙️ | The PBPK benchmark: seal → drive → judge. This replaces the old `seal_case`/`verify_no_cheat`/`judge_case`/`osp_agent_cli` scripts. |
| `run_qsp_full_pipeline` 🧪🔑 | QSP capstone: build immune network → calibrate → virtual pop → predict a real RA trial, graded vs the validated model. |
| `run_llm` 🔑 | The generic LLM decision-loop driver (popPK/PD and general). |
| `run_llm_build` ⚙️🔑 | **Legacy** scripted PBPK loop — superseded by the `pkpd-bench` sandbox flow; kept as a non-agent baseline. |

## Benchmark construction (PBPK) ⚙️

Turn published OSP models into sealed tasks and score them.

| Script | Role |
| --- | --- |
| `build_ddi_benchmark`, `build_metabolite_benchmark`, `build_biologic_benchmark` | Generate DDI / metabolite-cascade / biologic benchmark cases. |
| `blank_dissolution`, `reblank_methods` | Surgical blanking patches on existing snapshots. |
| `run_osp_scoreboard` | Run/aggregate the multi-model scoreboard. |
| `run_ddi`, `run_paper_compare`, `paper_compare` | DDI runs and paper-vs-model comparisons. |
| `check_harness_fidelity`, `published_gmfe`, `verify_benchmark`, `benchmark_coverage` | Harness sanity: reproduce published GMFEs, coverage checks. |
| `extract_snapshot`, `osp_run`, `fit_everything` | Snapshot extraction, a raw PK-Sim run, and the fit-everything anti-pattern baseline. |

## QSP pipeline 🧪🔑

The RA immune-network work. `run_qsp_full_pipeline` is the one-command capstone; the rest are
its stages and experiments.

| Group | Scripts |
| --- | --- |
| Build | `run_qsp_build_cell`, `run_qsp_build_network`, `run_qsp_build_general`, `run_qsp_build_honest`, `run_qsp_assemble`, `run_qsp_couple_general` |
| Calibrate / fit | `run_qsp_calibrate`, `run_qsp_calib_loop`, `run_qsp_fit`, `run_qsp_fair_fit`, `run_qsp_e2e_fit`, `run_qsp_l1_calibrate` |
| End-to-end / clinical | `run_qsp_end_to_end`, `run_qsp_e2e_clinical`, `run_qsp_e2e_coupled`, `run_qsp_agent_clinical`, `run_qsp_paper_pipeline`, `run_qsp_pipeline` |
| Virtual population | `run_qsp_vpop`, `run_qsp_agent_vpop`, `run_qsp_cohort`, `run_qsp_arm_targets` |
| Discover / analyze | `run_qsp_discover`, `run_qsp_spec`, `run_qsp_gsa`, `run_qsp_func_equiv`, `run_qsp_from_scratch_gap`, `run_qsp_probe_mtx`, `run_qsp_designed_drug`, `run_qsp_bench`, `run_qsp_l1_agent`, `run_qsp_l1_open` |
| LLM-driven QSP 🔑 | `run_llm_qsp_all`, `run_llm_qsp_design`, `run_llm_qsp_topology`, `run_llm_qsp_fit`, `run_llm_qsp_trial`, `run_llm_qsp_validate`, `run_llm_qsp_vpop_gen`, `run_llm_qsp_full` |

## LLM drivers (other) 🔑

| Script | Role |
| --- | --- |
| `run_llm_ddi` | DDI modeling loop. |
| `run_llm_extract`, `run_llm_extract_params`, `run_llm_draft_tasks`, `run_llm_task` | Data / parameter / task extraction from papers. |
| `run_llm_init`, `run_llm_onboard` | Project init and onboarding loops. |
| `run_llm_topology`, `run_llm_vpop_select` | Topology proposal, virtual-population selection. |

## Research / one-off

`demo_dry_run`, `demo_real_fit`, `dump_network`, `sbml_to_network_json`, `simbiology_smoke`,
`run_simbiology`, `extract_ra_provenance`, `extract_ra_targets`, `validate_new_tasks`,
`validate_processes` — demos, converters, and provenance/validation checks.

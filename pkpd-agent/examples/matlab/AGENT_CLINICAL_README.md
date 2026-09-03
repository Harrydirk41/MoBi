# Agent-built immune network → clinical sbproj → train / test / simulate

This is the bridge from the **agent-built immune network** (`mynet.xml`, emitted by
`run_qsp_build_network`) to the paper's **clinical train/test/simulate** workflow, run on
**your** machine (LLM agent + MATLAB/SimBiology).

## Why a transplant, not a from-scratch sbproj

The clinical shell — the DAS28-CRP composite (tender/swollen joint counts, CRP, patient global),
the CRP↔cytokine map, the PK compartments, and the drug dose objects — **is not in this repo**. It
lives only in the binary `.sbproj`. Reconstructing it would mean **inventing** the DAS28 formula and
the PK, which would be fabrication. So the honest bridge keeps the paper's given clinical shell and
only **swaps in the agent's immune mechanism**, connected by the fact that the agent network uses the
**same species names** as the paper (the nodes are the paper's own — IL6, TNFa, FLS, Th1, …). The
paper's `DAS28_CRP = f(immune state)` rule then re-reads the agent's dynamics automatically.

## Steps (in MATLAB, on your machine)

```matlab
% 0. emit the agent network first (Python side, with your API key):
%    python -m examples.run_qsp_build_network --model ra --live --prune --emit mynet.xml

% 1-4. build the agent-based clinical sbproj + baseline sanity sim:
addpath('examples/matlab');
sb_agent_clinical('Vantage RA QSP Model v1.0.sbproj', 'mynet.xml', 'agent_clinical.sbproj');

% then your usual train / test / simulate ON agent_clinical.sbproj:
sb_load('agent_clinical.sbproj');
% TRAIN    : sb_fit(paramSpec, 'trial_das28.csv', 'DAS28_CRP = das28', 'lsqnonlin', doses, 'fit.csv')
% TEST     : sb_run_vpop('Vpop1.xlsx', firstLineDoses, 285, 199, 284, 'first.csv', 300, '', '')
% SIMULATE : sb_run_vpop('Vpop1.xlsx', switchDoses,    601, 199, 600, 'second.csv', 300, '', '')
% COMPARE  : sb_paper_compare(...)   % agent-model response vs the RADIATE trial
```

`run_qsp_paper_pipeline.py --sbproj agent_clinical.sbproj --vpop Vpop1.xlsx --matlab` also drives the
whole PART II on the agent-based project.

## Seams to verify (this scaffold was written without a MATLAB instance)

Always run the **dry run first** — `sb_transplant_immune('mynet.xml', true)` — and read every field:

1. **Species-name match.** The transplant swaps only the species shared by name between `mynet.xml`
   and the sbproj. Watch aliased cells: the targets use `Macrophages`, but the sbproj species may be
   `Macrophage`. If `report.uncovered` lists a cell/cytokine you expected to swap, the names differ —
   rename in the emitter or the sbproj so they match, or that species keeps the paper's dynamics.
2. **`report.removeReactions`** — confirm these are all genuinely pure-immune (production/clearance of
   cytokines, cell life-cycle). Nothing clinical/PK should appear.
3. **`report.doubleDriven`** — reactions kept because they also touch a non-immune species but still
   produce/consume a transplanted one. Review by hand for double-counting.
4. **Baseline sim (step 3)** — `DAS28_CRP` must be finite and physiological. If it is non-finite, the
   transplanted network diverged; prune harder (`--prune`) and re-emit before fitting.

## Honest expectation

The agent network still over-includes edges (precision ≈ 0.3–0.5), so the clinical fit will be
**worse than the paper's own model**. That gap — quantified against the trial — is the result of the
experiment ("how close does a from-scratch agent build get to the validated model?"), not a bug.

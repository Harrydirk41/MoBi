# RA QSP virtual-trial agent task (Vantage RA model)

This is the SimBiology counterpart of the OSP PBPK agent benchmark. Where the
OSP tasks ask an agent to recover a *molecule's* parameters from PK data, this
task asks an agent to run an *in-silico trial* on a virtual population and
predict a clinical response rate — then checks it against a held-out arm.

The model is the **Vantage RA QSP Model v1.0** (`RA-QSP-Model/`): 59 species,
524 parameters, 101 reactions, and a virtual population of 300 patients whose
129 varied parameters live in `Vpop1.xlsx` / `Vpop2.xlsx`.

## The trial is encoded in the model as events (decoded from the .sbproj)

The `.sbproj` is a zip of HDF5 (`-v7.3`) `.mat` files. Decoding the event
triggers and rule expressions out of `simbiodata.mat` shows the entire clinical
trial is baked into the model as timed events — so the response is a **model
output we read, never something we recompute**:

| stage | event (from the model) | day |
|------|------------------------|-----|
| capture baseline | `DAS28_BL = DAS28_CRP` at `time>=199` | 199 |
| response math | `delta_DAS28_CRP = DAS28_BL - DAS28_CRP`; `ACR_Perc = 100*delta_DAS28_CRP/DAS28_BL` | — |
| first-line readout | `ACR20/50/70 = 1` when `ACR_Perc>=20/50/70` at `time>=284 & time<285` | 284 (wk 12) |
| remission | `Remission = 1` when `DAS28_CRP<=2.6` | 284 |
| flag MTX failure | `MTX_NonResp = 1` for inadequate responders | ~284 |
| second-line readout | `MTX_NonResp_TCZ_ACR20/50/70/Rem = 1` when `ACR_Perc>=… & MTX_NonResp==1 & time>=600` | 600 |

Key consequences:

* **Baseline is day 199, first-line readout is day 284 (week 12).** Treatment
  starts day 200 — every primary dose name carries a `_t200` suffix
  (`MTX_15mg_Q1W_SC_t200`, `TCZ8mgkg_Q4W_IV_t200`, `ADA40mg_Q2W_SC_t200`, …).
  A later switch dose exists too (`Ada40mg_t564`).
* **ACR is `100*(DAS28_BL - DAS28_CRP)/DAS28_BL`** — the Python side used to
  recompute exactly this; now it just reads the model's `ACR20/50/70` flags.
* **The flagship (TCZ in MTX-inadequate responders) is a built-in readout.**
  `MTX_NonResp_TCZ_ACR20/50/70/Rem` fire at day 600 for the `MTX_NonResp==1`
  subgroup — this is the paper's *held-out* validation, computed by the model.
  We do not need to reconstruct the two-step trial in Python; we run long enough
  (past day 600), then read the flags.

## Two arms, one runner

`examples/run_ra_vpop.py` runs the Vpop and reports the model's own flags. Both
readouts come out of a single run when the dose set includes MTX + TCZ and the
sim runs past day 600:

```
# first-line MTX only (week-12 ACR), quick 20-patient check
python -m examples.run_ra_vpop --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" \
    --vpop "..\RA-QSP-Model\Vpop1.xlsx" --dose MTX_15mg_Q1W_SC_t200 --limit 20

# full population, MTX first-line + TCZ — captures the second-line flagship too
python -m examples.run_ra_vpop --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" \
    --vpop "..\RA-QSP-Model\Vpop1.xlsx" \
    --dose "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200" --stop-time 700 --limit 300
```

The runner prints:
* DAS28 distribution at baseline (day 200) and first-line readout (day 284);
* **first-line** ACR20/ACR50/ACR70/remission rates (model flags, day 284);
* **second-line** TCZ ACR20/50/70/remission among the `MTX_NonResp==1` subgroup
  (model flags, day 600) — the flagship validation.

`sb_run_vpop.m` reads each patient's flags with a NaN-tolerant lookup
(`local_lastval`), so a run under a dose set that never triggers a given flag
just reports `n/a` for it rather than failing.

## Dose timing (decoded from the HDF5) and the sequential switch

Dumping the float64 scalars out of `simbiodata.mat` recovers the dose start
times even though SimBiology stores them in binary MCOS form: **13 doses start
at day 200** (every `_t200` name, including `TCZ8mgkg_Q4W_IV_t200`), exactly one
(`Ada40mg_t564`) starts at day 564, and `600` is the second-line readout gate.
So the shipped TCZ dose is **concurrent from day 200** — there is no late-start
TCZ dose in the project.

That concurrency is why applying `MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200`
gives a dead second-line arm (in a 300-patient run: first-line ACR20 43% /
remission 27% look right, but the flagship comes back ACR20 23%, ACR50/70 0%).
With both drugs from day 200, the `MTX_NonResp==1` patients are people who
failed the *combination* — the refractory tail — so there is no additional drug
to rescue them at day 600.

The fix is a **sequential switch**: give MTX from day 200 but start TCZ *after*
the day-284 first-line readout, so `MTX_NonResp` is a clean pure-MTX
classification and TCZ then acts on those non-responders through day 600. The
runner builds this with a `name@START` suffix that clones the dose with an
overridden `StartTime`:

```
python -m examples.run_ra_vpop --sbproj "...sbproj" --vpop "...Vpop1.xlsx" \
    --dose "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285" --stop-time 700 --limit 300
```

Sweep the TCZ switch day (285 = immediately after first-line; 564 = the switch
day the model itself uses for the ADA arm) and compare the flagship
ACR20/50/70 against the paper's reported held-out validation to pick the
protocol the authors used.

### Validated result (Vpop1, 300 patients, TCZ switched in at day 285)

`--dose "MTX_15mg_Q1W_SC_t200;TCZ8mgkg_Q4W_IV_t200@285" --stop-time 700 --limit 300`

| arm | ACR20 | ACR50 | ACR70 | Remission | n |
|-----|-------|-------|-------|-----------|---|
| MTX first-line (day 284, pure MTX) | 33.7% | 17.7% | 2.3% | 18.0% | 300 |
| TCZ in MTX-IR (day 600, flagship)  | 44.9% | 23.5% | 14.0% | 28.4% | 243 |

Both arms are clinically plausible (MTX-monotherapy week-12 ACR20 ~30-40%; real
TCZ-in-MTX/DMARD-IR trials report ACR20 ~48-59%, ACR50 ~32-44%, ACR70 ~12-22%).
The MTX-IR denominator (243/300) uses the stricter EULAR criterion
(deltaDAS28 < 1.2), not ACR20, so it is larger than the ACR20 non-responder count.

**Read the flags at the right time.** ACR20/50/70/Remission are CONTINUOUS flags
(trigger `ACR_Perc>=X`): they must be read AT the day-284 first-line readout, not
at sim end - reading at day 700 wanes for untreated patients and is contaminated
by any second-line drug, which silently made the first-line number track the
second-line drug's potency. `MTX_NonResp` (latched day 284) and the
`MTX_NonResp_TCZ_*` flagship (day 600) are read at sim end where they hold.

**Invariance check (the harness-correctness gate).** With the first-line flags
read at day 284, the first-line rates are IDENTICAL across second-line arms
(MTX-only = MTX+SEC@285 = MTX+TCZ@285 = ACR20 30% at n=50), proving the day-284
readout does not depend on a drug started at day 285. Only the second-line column
moves with the drug: SEC is correctly inert (the model reproduces secukinumab's
real-world RA failure), TCZ is alive (~42/26/16 at n=50, ~45/24/14 at n=300).
Two SimBiology gotchas were ruled out along the way: `OutputTimes` does not drive
event detection (the solver watches triggers at its own steps), and the
narrow-window events (`time>=284 & time<285`) fire fine without a `MaxStep` cap -
the only bug was the read time.

**One run reproduces both arms, correctly separated** - the MTX first-line and
the held-out TCZ-in-MTX-IR flagship - so the pipeline is validated and ready to
carry the agent loop.

## How this becomes the agent task

Fill-in-the-blank style (the style chosen for the OSP tasks):

1. Give the agent the model, the Vpop, and the **calibrated** arms (MTX, ADA,
   TCZ first-line) with their observed response rates.
2. Blank the **held-out** arm (TCZ in MTX-inadequate responders).
3. The agent proposes the regimen (dose set + schedule), runs the Vpop, reads
   the model flags, and predicts the second-line ACR20/50/70.
4. Score = agreement between the predicted rates and the model's held-out
   `MTX_NonResp_TCZ_*` rates (equivalently, the paper's reported validation).

The agent's decision loop is the same shape as the OSP loop — inspect the model
and available regimens, choose a protocol, run, evaluate the population readout,
and stop when the prediction is both close to the target and mechanistically
sensible (right drug, right line of therapy, plausible dose).

## Porting to a new QSP model

The engine (`SimBiologyEngine`), the loop, and the numerical routines
(`numeric_fit_1d`, `select_to_moments`, `ir_mask`, dose retime/scale, structural
drug-add) are **model-agnostic**. What is specific to the Vantage RA model is
gathered into one adapter — `engines/qsp_config.py` (`QSPModelConfig`, instance
`VANTAGE_RA`):

* `readout_states` — the 11 model state names, in role order, that `sb_run_vpop`
  reads. The CSV column *roles* are fixed (four first-line flags, five latched
  flags, a trajectory state, a baseline state); the config maps each role to *this
  model's* state name. `sb_run_vpop.m` takes these as its `stateSpec` argument and
  falls back to the RA defaults when none is given, so a new model changes the
  names without touching the `.m`.
* `timeline` — the days the model's readout events fire.
* `drugs` / `vpop_drivers` / `fit_params` / `design_targets` — the catalogs the
  task tools present to the agent.
* `clinical_reference` — the real trial data to score against.

So a new QSP model is **a new `QSPModelConfig` plus its reference data**, not a
code rewrite — the engine, loop, and numerical routines carry over unchanged. What
is *not* yet automated: deriving the catalogs from the model itself. `model_info()`
already lists any model's species / parameters / doses generically, so the next
step toward plug-and-play is having the agent propose the config from that
inventory rather than a human writing it.

## Two objectives (difficulty)

`run_llm_ra_trial.py --objective` selects the task:

* **`predict`** (default) — reproduce the held-out TCZ-in-MTX-IR rates. The
  decision space is small (~a dozen drug×timing combinations) and a "pick the
  strongest arm" heuristic happens to land on the right drug, so this is really a
  pipeline/​reasoning demonstration, not a hard search.
* **`min-dose`** — find the LOWEST second-line dose that still clears an ACR20 bar
  (`--min-dose-acr20`, default 35%). `dose_scale` makes the dose axis continuous,
  so "just max the dose" clears the bar but wastes drug (and scores worse on
  response-per-dose), while too little misses it. The agent must titrate to the
  efficient point — brute-force enumeration over a continuous dose no longer
  works, and the shortcut proxy fails.

The dose spec grew a `*scale` suffix (multiply the dose amount) beside `@day`
(retime the start): token form `NAME[*SCALE][@START]`, e.g.
`TCZ8mgkg_Q4W_IV_t200*0.5@285` is half-dose TCZ switched in at day 285.

**De-hinted tools.** The `ra_inspect`/`ra_run_trial` descriptions no longer name
the switch day (285), the sequential-vs-concurrent trade-off, or the expected
response magnitude. Those were giveaways; the agent now has to discover that a
concurrently-dosed second line contaminates the day-284 classification, and to
derive the drug/timing/dose itself from mechanism and experiment.

## Scoring against REAL clinical data (not just the model)

`engines/ra_clinical_reference.py` transcribes the paper's ESM1
(`Clinical_trials` sheet) - the actual trials the model was calibrated to:

| drug | trial | ACR20/50/70 (raw) | placebo-corrected |
|------|-------|-------------------|-------------------|
| MTX  | Strand 1999, wk12 | 46 / 23 / 9 | 20 / 15 / 5 |
| ADA  | OPTIMA, wk24 | 70 / 52 / 35 | 13 / 18 / 18 |
| TCZ  | ROSE, wk24 | 45 / 29 / 13.9 | 20 / 19 / 12 |

`run_llm_ra_trial.py --target-source clinical` (the default) scores the agent's
TCZ prediction against the **real ROSE trial**, not the model's own output. Two
things matter:

* The **agent's self-validation** still uses the *model's* MTX output (33.7%),
  because that check asks "does my pipeline reproduce the model" - swapping in the
  real 46% would make a correct harness look broken.
* The **final score** uses the *real* data, because that asks "does the
  model/agent match reality". The agent's committed prediction (42.1/26.3/15.8)
  lands **MAE 2.5 pp from real ROSE wk24 (45/29/13.9)**; the model's own output
  (44.9/23.5/14.0) is 1.9 pp. Placebo-corrected is ~11 pp off, confirming the QSP
  Vpop output aligns with *raw* clinical efficacy (it carries its own disease
  baseline). The run prints all three reference frames side by side.

## The agent loop (implemented)

Same three pieces as the OSP DDI loop:

* `pkpd_agent/engines/osp_ra_trial.py` — pure-Python: `summarize_run` (a Vpop CSV
  → first-line and second-line response rates), `build_dose_spec` (the
  `name@switch_day` sequential-switch string), `score_flagship` (predicted vs
  held-out, MAE in percentage points), and `DRUG_CATALOG` (the formulary with
  mechanisms).
* `pkpd_agent/tools/ra_trial_loop_tools.py` — `ra_inspect` (observe: disease,
  timeline, formulary, calibrated reference arms, the held-out objective) and
  `ra_run_trial` (act: apply a `{first_line, second_line, switch_day}` protocol,
  run the Vpop, return the model's response rates). The held-out target is **not**
  exposed to the agent.
* `examples/run_llm_ra_trial.py` — wires the SimBiology engine + tools + the LLM
  policy, then scores the agent's final protocol against the held-out truth.

```
set ANTHROPIC_API_KEY=...
python -m examples.run_llm_ra_trial ^
    --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" ^
    --vpop   "..\RA-QSP-Model\Vpop1.xlsx" --limit 50 --max-steps 8
```

The agent is told the disease, the timeline, the formulary (MTX + adalimumab /
tocilizumab / secukinumab / anakinra, each with its mechanism) and the calibrated
MTX first-line rates — but **not** which second-line drug to use, nor the held-out
answer. A correct run: validate the harness on first-line MTX, reason that
MTX-inadequate responders need a mechanistically distinct biologic switched in
after day 284, run `MTX_15mg_Q1W_SC_t200` + `TCZ8mgkg_Q4W_IV_t200@285`, and predict
ACR20 ~45% / ACR50 ~24% / ACR70 ~14%. Choosing concurrent dosing, or a
TNF/IL-17/IL-1 agent instead of the IL-6 blocker, moves the prediction away from
the held-out truth and raises the MAE. Keep `--limit` modest (50) during the loop —
each `ra_run_trial` simulates the whole subsampled population — then confirm the
winning protocol at `--limit 300` with `examples.run_ra_vpop`.

Scoring targets default to the validated model output (ACR20 44.9 / ACR50 23.5 /
ACR70 14.0); pass `--target-acr20/50/70` to score against the paper's reported
held-out numbers once you have them.

**Subsampling is representative, not first-N.** The Vpop rows are ordered by
disease severity — the first 50 patients have inflammatory drivers (IL-17,
RANTES, GM-CSF, VEGF secretion, FLS baseline) shifted +0.3 to +0.6 SD toward
*more severe* disease, and taking the first 50 floors the flagship at
ACR50/70 = 0 (even TCZ 8 mg/kg gives only ACR20 28.6% / ACR50 0% on that slice
vs 44.9% / 23.5% on the full 300). So `--limit N` runs an evenly-spaced sample
across the whole population, which tracks the full-population rates much more
closely and gives the agent a fair, scorable signal.

---

# Stage 1: can an LLM do the model *building*? (a benchmark, not a claim)

Everything above (the "2–7" tasks: trial design, calibration, Vpop generation,
drug design, validation) operates on the **finished** model. Honest accounting of
those tasks: the LLM is rarely load-bearing — a numerical routine (scipy) or the
hard-coded tool logic does the work, and the LLM mostly *selects among options we
laid out and narrates the result*. They show an agent can carry the workflow, not
that it does the science.

Stage 1 — building the model itself (its network + parameters) — is the opposite:
open-ended, no ground truth until the whole thing is done, and not brute-forceable.
It is the paper authors' real contribution and we did **none** of it. But because we
have the finished model as an **answer key**, we can benchmark whether an LLM could
*reconstruct* pieces of it — which is where LLM reasoning is finally load-bearing.

We measured several layers, top (qualitative) to bottom (quantitative). Every number
below is the model's own scope/wiring/values as the key; the LLM never sees it until
scoring. Headlines are over **5 runs** (`--repeat 5`) — variance is small enough that the
numbers are stable ceilings, not n=1 luck.

## Layer 0 — model scope (`run_llm_ra_scope`, `ra_scope.py`)

Given only the disease and the modeling goal (no cast), the agent proposes which cells and
mediators to include; scored precision/recall/F1 vs the model's real 26-node cast. Over-
inclusions that are real RA mediators the model deliberately omits (IL-2/8/15/18/32,
NK/dendritic/osteoclast, …) are flagged so the miss list is interpretable. This is the
Stage-2 scoping judgment the other tasks skipped by handing the agent the cast.

Result (5 runs): **F1 0.698 ± 0.028** (P 0.68 / R 0.72). The agent knows the cast; its
*genuine* errors are modeling-scope judgment, not biology recall: it under-includes the
recruitment/trafficking layer (Endo, CAM, chemokines MIP3/RANTES) and over-includes the
joint-erosion axis (chondrocyte, osteoclast, RANKL, MMP3) that the ACR/DAS28 endpoints
never read out — i.e. it scoped to "RA pathology" instead of "what drives *these*
endpoints." (First pass scored a spurious 0.44: the scorer, not the agent, was wrong —
free-text names like "Th1 cell", "Fibroblast-like synoviocyte", "CCL2 (MCP-1)" failed a
strict match and were double-penalised as miss+extra. `resolve_node` fixed it. This was
the third time the harness under-credited the agent; see the caveat below.)

## Layer 1 — network topology (`run_llm_ra_network`, `ra_network.py`)

Given only the cast (9 cells, 17 cytokines incl. TGFb/IL10), the agent proposes the
signed regulatory edges; scored precision/recall/F1 vs the model, both sign-aware and
topology-only. The answer key is parsed from the model's rule expressions — each
`MM(source, …)` term in a `(Pro|Anti)_<TargetProcess>_effect` rule is one edge — **not**
from the `_by` parameter names (that first key was 3× short: 30 vs 88 edges).

Result (full 88-edge key):
* biology-only prompt (5 runs): **topology F1 0.555 ± 0.017** (P 0.49 / R 0.63,
  ~50/80 recovered) — the sd is tiny: the LLM reconstructs almost the same network every
  run, so ~0.55 is a **stable ceiling**, not a lucky draw.
* convention-primed (`--conventions`: told the model is bipartite cell↔cytokine +
  negative self-feedback): **topology F1 ≈ 0.66** (R 0.75, 60/80)
* sign-aware F1 0.447 ± 0.013 — **signs are the weak point**; primed about signs,
  the agent hedges (proposes both signs) rather than localizing them.

Read: an LLM reconstructs the **skeleton** well (a genuinely useful first-pass draft),
and its self-diagnosis is correct and substantive (it identified that the model is a
bipartite graph with self-limiting feedback loops — a real modeling-convention insight).
But ~40% of its edges are spurious and half the signs are wrong, so it drafts, it
doesn't build.

## Layer 2 — parameter values (`run_llm_ra_params`, `ra_params.py`)

Given each parameter's name/units/cell-context, the agent predicts its value; scored
order-of-magnitude (log10 error) vs ESM2's 130 documented values. **The honesty check
is a naive unit-geomean baseline** (knows each unit's empirical scale, zero biology).

Result:
* overall "beats baseline: True" — but that is a **mirage**: it is carried by the 90
  **dimensionless** fold-effects, which cluster near 1 so the baseline already scores
  99% within 10×.
* on the 40 **dimensional** parameters (rates, secretion, concentrations — where
  physiology should help): LLM median log10 err **0.78 ± 0.05 (5 runs) vs baseline 0.62**
  — the LLM is **reliably worse** than the dumb guess, every run. This is the wall.
* the worst misses are not bad biology — they are model-internal normalization
  (molecule-scaled secretion ~1e-9, unknowable without the model) and unit conventions
  (koff given in 1/s vs the model's 1/day). The agent's own reasoning was often
  physiologically right and lost only on the model's bookkeeping.

## The ladder (the whole finding, all layers measured over 5 runs)

The distinguishing axis is **biology-determined vs model-committed**, not qualitative vs
quantitative. Where the answer is dictated by known RA biology, the LLM is strong; where it
is the modeler's own committed choice (a fitted value, a sign convention, an arbitrary
readout formula), it is weak — even when the task is otherwise a simple "selection".

| the answer is… | layer | LLM (5 runs) | baseline | verdict |
|---|---|---|---|---|
| **known RA biology** | sensitivity — *which* params matter | recall 0.87 ± 0.03 | random 0.44 | **strong** |
| | model scope (cast), primed | F1 0.82 ± 0.06 | — | strong |
| | model scope (cast), raw | F1 0.70 ± 0.03 | — | strong |
| | network topology, primed | F1 ~0.66 | — | useful |
| | network topology, raw | F1 0.555 ± 0.017 | — | useful **draft** |
| **this model's committed choice** | readout formula (DAS28 drivers) | F1 0.37 ± 0.06 | 9-cell key | weak |
| | edge signs (isolated) | acc 0.77 ± 0.005 | majority 0.74 | barely above guessing |
| | sensitivity — exact *rank* | Spearman 0.34 | — | weak |
| | parameter values (fair physiological subset) | 0.76 ± 0.04 | 0.68 | **below baseline (0/5)** |

**The LLM identifies what standard biology dictates (which cells belong, what wires to
what, which knobs dominate) and fails at what the model committed to on its own** — the
fitted parameter values, the exact feedback signs, the precise sensitivity order, and the
idiosyncratic readout formula (9 specific cells, T-cells split into 4 subsets, no cytokines,
Treg negative). The readout is the sharpest case: the agent gave the *biologically* sensible
answer (cell load + IL-6→CRP) and lost only because the model made *different* reductive
choices it had no way to guess. Priming lifts the biology-adjacent layers (scope 0.70→0.82,
topology 0.55→0.66) because it supplies the model's conventions; it cannot rescue the
committed-value layers, because those are not derivable.

Recurring methodological finding: the **harness under-credited the LLM four separate times**
(edge key 3× short; missing TGFb/IL10; scope free-text synonyms; readout free-text synonyms)
— a strict-match answer key systematically understates a free-text LLM, and each time the
agent's own self-diagnosis caught it. The tolerant `resolve_node` matcher and the naive/
majority/random baselines are most of what makes these numbers trustworthy.

Caveats kept explicit: (a) headlines are over 5 runs and the variance is small (F1 sd
0.02–0.03) — the numbers are stable, not n=1 noise; (b) this is **one** QSP model — the
ladder is a finding about this RA model, and would need a second model to generalize; (c)
the benchmark tests *reconstruction against an answer key*, which is easier than *de novo
construction* — recovering wiring you're scored against is not the same as building a model
that simulates and calibrates; (d) **the harness under-credited the agent three times**
(network key 3× short; missing TGFb/IL10 from the cast; scope free-text names failing a
strict match) — each time the raw number understated the LLM and the agent's own self-
diagnosis flagged the bug. A strict-string answer key systematically penalises a free-text
model; tolerant (or LLM-judge) matching is required for an honest score.

## More Stage-1 probes (built; run to measure)

Five further benchmarks push on the layers above and add two new skills. Each has a
`--repeat N` variance mode; the LLM runs need `ANTHROPIC_API_KEY` (no live MATLAB).

* **Fairer parameter test** (`run_llm_ra_params`, now splits the 40 dimensional params):
  a **physiological** subset (28: rates `1/day`/`sec-1`, concentrations `M` — groundable
  from biology; baseline median 0.68) vs a **model-scaling** subset (12: per-molecule /
  per-mL normalization — unknowable). `beats_physiological_baseline` is the *fair* verdict;
  the earlier all-dimensional 0.62 unfairly lumped in the unknowable units.
* **Isolated sign prediction** (`run_llm_ra_sign --network network.json`): hand the agent
  the true *unsigned* edges, ask only activate-vs-inhibit; scored vs the **majority-class
  baseline** (most edges activate, so "all +1" already scores high — the bar to beat).
* **Sensitivity ranking** (`run_llm_ra_sensitivity`): a *different skill* — which knobs
  matter, not biology recall. Rank the parameters driving DAS28-CRP from a pool of 50 (the
  paper's Fig-9 GSA top-20 hidden among 30 real distractors); scored overlap + rank
  correlation vs the GSA, against a **random blind-pick baseline (recall 0.40)**.
* **Readout mapping** (`run_llm_ra_readout --network network.json`): the mechanism→endpoint
  bridge — which nodes DAS28-CRP is computed from; scored vs the species the model's readout
  rule depends on. `--show-key` prints the extracted drivers + raw rule.
* **Scope priming** (`run_llm_ra_scope --conventions`): hands the agent the endpoint-focus +
  trafficking-layer conventions it diagnosed missing, to separate scope *judgment* from
  recall (analogous to `--conventions` on topology).

## Files

* `pkpd_agent/engines/ra_scope.py` / `tools/ra_scope_loop_tools.py` /
  `examples/run_llm_ra_scope.py` — Layer-0 scope/cast benchmark (`--repeat`,
  `--conventions`). `resolve_node` does tolerant free-text matching.
* `pkpd_agent/engines/ra_sensitivity.py` / `ra_readout.py` + their `tools/` and
  `examples/run_llm_ra_sensitivity.py` / `run_llm_ra_readout.py`; sign task in
  `tools/ra_sign_loop_tools.py` + `examples/run_llm_ra_sign.py` (`score_signs` in
  `ra_network.py`).
* `pkpd_agent/engines/ra_network.py` / `tools/ra_network_loop_tools.py` /
  `examples/run_llm_ra_network.py` — topology benchmark (`--conventions`, `--repeat`).
* `pkpd_agent/engines/ra_params.py` (+ `data/ra_params_esm2.json`) /
  `tools/ra_params_loop_tools.py` / `examples/run_llm_ra_params.py` — parameter
  benchmark with the naive-baseline honesty check (`--repeat`).
* `examples/matlab/sb_network_json.m` + `examples/dump_network.py` — dump the full
  wiring answer key from the model (run once).

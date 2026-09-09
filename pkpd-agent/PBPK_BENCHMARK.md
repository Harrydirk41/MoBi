# PBPK agent benchmark — complete guide

Tests whether an LLM agent can build whole-body PBPK models as well as the published
OSP expert models. It grades the agent the way OSP qualifies its own models: **fit on
part of the data, predict the held-out part, and be judged relative to the expert
model.** Everything below runs with PK-Sim (`--pksim <PKSim.CLI.exe>`).

---

## The one command

```powershell
cd pkpd-agent
python -m examples.run_osp_scoreboard --run --aggregate ^
    --pksim "C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe"
```

This runs the agent on every gradeable single-compound model (easy→hard), and writes
`osp_scoreboard.md`. Add `--models A,B,C` to run a subset, or `--hard` for the
mechanism-discovery variant (below). One model takes a few minutes to ~15 min.

---

## How grading works (held-out, relative to the reference)

1. **Split** each model's studies into **building** and **verification**, by whole
   study (never splitting a study). IV studies stay in building (they anchor
   disposition); ~1/3 of the rest are held out. OSP's own building/verification tag is
   used where present.
2. The **agent fits on the building studies only** — it never sees the verification
   concentrations, only every study's dosing design. Its model therefore *predicts* the
   held-out arms unseen.
3. **Grade** = the agent's GMFE on the verification studies vs the **reference model's**
   GMFE on the same verification studies. **PASS = agent ≤ 1.25× reference.** Any valid
   parameterization that predicts held-out data about as well as the expert passes — the
   grade does not require matching the reference's exact parameter values.

Scoreboard columns: `held-out agent | held-out ref | ratio | structure | all-data agent
| verdict`. The all-data (in-sample) GMFE is context only. Models with too few studies
to hold out (Sufentanil, Vancomycin, Montelukast) are graded on all data and flagged,
not counted as held-out passes.

---

## Two difficulty modes

- **Normal** (`run_osp_scoreboard --run`): the agent is given the mechanism (which
  enzymes clear the drug) and must choose distribution/permeability methods and fit the
  parameters.
- **Hard / discovery** (`run_osp_scoreboard --hard --run`): the metabolizing enzyme
  identity is **withheld**. The agent gets only a candidate pool of expressed enzymes and
  must reason out which one clears the drug from its chemistry/pharmacology. See
  `OSP-PBPK-Model-Library/HARD_DISCOVERY.md` — 10 of 19 models are genuine discovery
  tests (the rest have no distractor enzymes).

The agent tools **never** see the reference model or the answer — only the blanked
snapshot, the building data, and the biology facts. Verified by `verify_benchmark.py`.

---

## The judge is trustworthy (why the numbers mean something)

Before grading the agent, the harness must reproduce each **reference** model's own
published GMFE — otherwise a bad score could just be the judge miscounting. It does:

```powershell
python -m examples.check_harness_fidelity --pksim "...PKSim.CLI.exe"
```

14 of 15 adult references reproduce their published GMFE within ~1.1× (most exactly),
using the explicit observed→simulation linkage OSP records in each snapshot
(`OutputMappings`). The published yardstick is parsed by `examples.published_gmfe`.

The fix history that got here (all in `osp_score`): pairing each observation to its
simulation must respect study, route, dose (incl. unit-less and compact names),
formulation, food state, infusion duration, metabolizer phenotype, and dosing schedule
— and, above all, the **explicit `OutputMappings` linkage**, which makes pairing exact
and drops extra DDI-network arms the reference was never built on.

---

## Known limitations (honest)

- **Propofol** is excluded: its target-controlled-infusion arms are not reproduced by
  the harness, so its reference GMFE is not trustworthy. Would need infusion-protocol
  support.
- **Small-data models** (Sufentanil, Vancomycin, Montelukast, Raltegravir-pediatric)
  have too few studies to hold out and are graded on all data.
- **Absorption**: many models fit Weibull dissolution per formulation. The normal grade
  handles this through the held-out predictive fit; a dedicated absorption-fitting step
  in the agent loop is the remaining enhancement (phase 3).
- The held-out GMFE and the in-sample GMFE are usually close, because PBPK models have
  few (2–8) physiological parameters — so the held-out split confirms generalization
  rather than exposing large overfitting.

---

## Files

| file | role |
|---|---|
| `examples/run_osp_scoreboard.py` | run the agent across models, write the scoreboard |
| `examples/run_llm_build.py` | one model: build + fit (building) + grade (held-out) |
| `examples/check_harness_fidelity.py` | verify the judge reproduces published GMFEs |
| `examples/published_gmfe.py` | parse each report's published GMFE (the yardstick) |
| `examples/verify_benchmark.py` | confirm no answer leak; `--hard` checks discovery |
| `pkpd_agent/engines/osp_split.py` | building/verification study split |
| `pkpd_agent/engines/osp_score.py` | observation↔simulation linkage + GMFE scoring |
| `OSP-PBPK-Model-Library/HARD_DISCOVERY.md` | the mechanism-discovery variant |

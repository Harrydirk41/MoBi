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
| MTX first-line (day 284, pure MTX) | 43.3% | 23.7% | 13.7% | 26.7% | 300 |
| TCZ in MTX-IR (day 600, flagship)  | 44.9% | 23.5% | 14.0% | 28.4% | 243 |

Both arms are clinically plausible (real TCZ-in-MTX/DMARD-IR trials report ACR20
~48-59%, ACR50 ~32-44%, ACR70 ~12-22%). Two checks confirm the split is real,
not lucky: moving TCZ to day 285 cut the first-line DAS28 drop from 1.86 (combo)
to 0.93 (MTX-mono size), and the MTX-IR denominator rose from 173 to 243 -
because `MTX_NonResp` uses the stricter EULAR criterion (deltaDAS28 < 1.2), not
ACR20, so weaker pure MTX leaves more non-responders. The concurrent-dosing
variant gives a dead flagship (ACR50/70 = 0%) and is the wrong protocol.

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

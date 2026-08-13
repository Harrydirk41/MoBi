# Setup on Windows — wiring nlmixr2, pharmpy, and OSP into the agent

This gets the three real engines talking to the agent from a **conda
environment** on Windows. After each step, run the health check to watch
capabilities light up:

```bash
python -m pkpd_agent.doctor
python -m pkpd_agent.doctor --rscript "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
```

The agent is Python; the engines are Python (pharmpy) or R (nlmixr2, ospsuite)
or .NET (OSP native). The agent **bridges out** to each. You never open a GUI.

---

## Step 0 — the base (always works)

```bash
conda create -n pkpd python=3.11 -y
conda activate pkpd
cd pkpd-agent
pip install numpy scipy            # the real pkfit engine
python -m unittest discover -s tests
python -m examples.demo_real_fit   # real fits + model selection, no backend needed
python -m pkpd_agent.doctor        # expect numpy/scipy OK, rest pending
```

At this point you already have a **real** estimation engine (`pkfit`) and the
full decision loop. Everything below adds the heavier backends.

---

## Step 1 — Claude (the decision brain)

```bash
pip install anthropic
setx ANTHROPIC_API_KEY "sk-ant-..."   # new shell afterwards
```
Verify: `python -m pkpd_agent.doctor` → `python:anthropic OK`, `ANTHROPIC_API_KEY set`.

Run it: `python -m examples.run_llm "Load the builtin dataset and build a popPK model."`

---

## Step 2 — nlmixr2 (real NLME population PK)  ← https://github.com/nlmixr2/nlmixr2

nlmixr2 is R + a C/C++ compiler. No NONMEM needed.

1. **Install R** (https://cran.r-project.org) — note the path, e.g.
   `C:\Program Files\R\R-4.4.1\`.
2. **Install Rtools** (matching your R version, from CRAN) — this is the
   compiler nlmixr2/rxode2 need. Confirm:
   ```bash
   "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" -e "pkgbuild::has_build_tools(debug=TRUE)"
   ```
3. **Install nlmixr2** (CRAN):
   ```bash
   "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" -e "install.packages('nlmixr2')"
   ```
4. **Verify end-to-end:**
   ```bash
   python -m pkpd_agent.doctor --rscript "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
   ```
   Expect `R:Rscript OK` and `R:nlmixr2 OK`.

**Use it from the agent** (real NLME with random effects):
```python
from pkpd_agent import AgentConfig, DecisionLoop
from pkpd_agent.policies import PharmacometricPolicy

cfg = AgentConfig(
    mock=False,
    rscript_path=r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe",
    nlmixr2_est="focei",     # or "saem"
)
# an LLM run can now choose nlmixr2_fit; or call the tool directly:
from pkpd_agent.tools import build_default_registry
from pkpd_agent.state import ModelingSession
reg = build_default_registry(cfg); s = ModelingSession(goal="fit")
reg.dispatch("pkfit_load_data", {"source": "builtin"}, s)
print(reg.dispatch("nlmixr2_fit", {"model": "1cpt_oral"}, s).data)
```
> The R worker is `pkpd_agent/engines/r_workers/nlmixr2_fit.R`. It builds a
> 1-/2-compartment `linCmt()` model with IIV on CL and V. Tune the model block
> and result fields to your data if needed — it's plain R.

---

## Step 3 — pharmpy (model building on top of a backend)  ← https://github.com/pharmpy/pharmpy

pharmpy is pure Python (`pip`), but its **fitting** needs a backend — and
nlmixr2 (Step 2) is a valid one, so this composes with Step 2.

```bash
pip install pharmpy-core
python -c "import pharmpy; print(pharmpy.__version__)"
python -m pkpd_agent.doctor            # python:pharmpy OK
```

The agent's `pharmpy_*` tools call the real pharmpy API when `mock=False` and
pharmpy is importable (`pkpd_agent/tools/pharmpy_tools.py`). Point pharmpy's
model-fit / AMD at nlmixr2 as its estimation tool per pharmpy's docs. Use
pharmpy for **AMD / structural search / covariate search**; use `nlmixr2_fit`
for a direct single fit.

---

## Step 4 — OSP (mechanistic PBPK / QSP)  ← the `ospsuite` R package

OSP has no Python API; you drive it through R (`ospsuite`).

1. **Install the `ospsuite` R package + dependencies** (it is **not** on CRAN —
   the desktop installer does NOT install it):
   ```r
   install.packages("remotes")
   remotes::install_github("Open-Systems-Pharmacology/rSharp")
   remotes::install_github("Open-Systems-Pharmacology/OSPSuite.RUtils")
   remotes::install_github("Open-Systems-Pharmacology/OSPSuite-R")
   ```
   (`rSharp` bridges R to .NET — install the .NET runtime it asks for.)
2. **Verify:**
   ```bash
   python -m pkpd_agent.doctor --rscript "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
   ```
   Expect `R:ospsuite OK` and `.NET runtime OK`.
3. **Smoke-test a real simulation:**
   ```bash
   "C:\...\Rscript.exe" -e "library(ospsuite); s<-loadSimulation(system.file('extdata','simple.pkml',package='ospsuite')); cat(class(runSimulations(s)[[1]]))"
   ```

**Use it from the agent** (`mock=False`, same `rscript_path`):
```python
cfg = AgentConfig(mock=False, rscript_path=r"C:\...\Rscript.exe")
reg = build_default_registry(cfg); s = ModelingSession(goal="pbpk")
reg.dispatch("osp_load_snapshot", {"path": r"C:\models\my_sim.pkml"}, s)
print(reg.dispatch("osp_simulate", {"snapshot_id": r"snap::C:\models\my_sim.pkml"}, s).data)
```
> Worker: `pkpd_agent/engines/r_workers/osp_sim.R` (load → set params → run →
> summarize output curve). Adjust output/parameter paths to your model. For
> snapshot↔JSON conversion use `MoBi.CLI snap` (Windows exe) — set
> `config.mobi_cli_path`.

---

## The single switch

Everything is gated by two config fields:

```python
AgentConfig(
    mock=False,                       # use REAL engines instead of synthetic
    rscript_path=r"C:\...\Rscript.exe",  # the R that has nlmixr2 + ospsuite
)
```

- `mock=True` (default): synthetic results, no backend needed — good for
  developing the loop/logic.
- `mock=False`: pkfit stays real (numpy/scipy); nlmixr2 and OSP call R; pharmpy
  calls Python. Each errors clearly (not crashes) if its backend is missing.

## Which engine for which job

| Job | Tool | Backend |
|---|---|---|
| Model-free first pass | `pkfit_nca` / `nca_analyze` | none (real) |
| Fast structural fit, model comparison | `pkfit_fit` | numpy/scipy (real) |
| **True NLME popPK** (random effects) | `nlmixr2_fit` | R + nlmixr2 |
| Automatic model development, covariate search | `pharmpy_*` | pharmpy (+ nlmixr2) |
| Mechanistic **PBPK / QSP** simulation | `osp_*` | R + ospsuite (+ .NET) |

## Troubleshooting

| Symptom | Fix |
|---|---|
| `Rscript not found` | pass the full `rscript_path`; `where Rscript` to locate it |
| `nlmixr2 package not installed` | `Rscript -e "install.packages('nlmixr2')"` |
| nlmixr2 fit errors on compile | install/repair **Rtools** (Step 2.2) |
| `ospsuite ... not installed` | you installed the desktop app, not the R package (Step 4.1) |
| `rSharp`/.NET load error | install the .NET runtime rSharp requires; re-run `library(ospsuite)` |
| tool returns `ok=False` with a message | that's the guarded error path — read the message; the loop keeps going |

Run `python -m pkpd_agent.doctor --rscript <path>` after every step; when the
rows you need are `OK`, the agent uses those engines for real.

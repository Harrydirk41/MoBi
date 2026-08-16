# Validating on a machine with PK-Sim

Everything in this repo is unit-tested **without** PK-Sim (structure analyzers,
scorers, benchmark generators, leak checks - 173 tests). A few things depend on
what PK-Sim actually emits at runtime (Results.csv column names, exact process
InternalNames) and can only be confirmed on a box that has PK-Sim installed.
Run the checks below once there; each prints what matched so any mismatch is
visible, and how to read off the correct name if something is off.

```
set PKPD_PKSIM_CLI=C:\Program Files\Open Systems Pharmacology\PK-Sim 12.3\PKSim.CLI.exe
cd pkpd-agent
```

## 1. New task engines end-to-end (metabolite / biologic / DDI)

```
python -m examples.validate_new_tasks --lib ..\OSP-PBPK-Model-Library
```

Runs one real snapshot per task type and reports matched molecules / matrices /
arms + a GMFE. What each confirms:

| task | snapshot | confirms | main risk |
|---|---|---|---|
| metabolite | Itraconazole | parent + each daughter pulled from ONE run (per-molecule column match) | low - molecule name is in the column |
| biologic | BAY794620 | (organ, compartment) column TOKENS match the real biodistribution column names | **medium - tissue column names are inferred** |
| ddi | Erythromycin | victim plasma column resolves in a multi-compound run | low - already exercised |

If a **biologic matrix comes back UNMATCHED**, the tissue column tokens don't
match PK-Sim's real names. Read the real names off the CSV header:

```
python -m examples.validate_new_tasks --only biologic --keep
```

then adjust the tokens in `osp_biologic.matrix_specs` (or the organ/compartment
labels) to match, and re-run. The matcher is space-insensitive substring on the
column header, so usually only an organ spelling needs aligning.

## 2. Inferred process InternalNames (single-compound catalog)

Nine addable process types carry `internal_name_verified=False` (named by PK-Sim
convention, not yet read off a real snapshot). Confirm they build:

```
python -m examples.validate_processes --snapshot ..\OSP-PBPK-Model-Library\Alfentanil\benchmark\Alfentanil-Model.blanked.json --unvalidated
```

It adds each unvalidated mechanism to a suitable expressed molecule, builds, and
runs one simulation - printing `BUILDS` / `FAILED` and the exact flag flips to
apply in `osp_catalog`. Send the SUMMARY block back and the verified flags /
corrected names get set in one pass.

## 3. Full benchmark regeneration (optional)

Regenerate every task type's files from the library (dry run drops `--write`):

```
python -m examples.build_benchmark            ..\OSP-PBPK-Model-Library --write
python -m examples.build_benchmark --hard     ..\OSP-PBPK-Model-Library --write
python -m examples.build_metabolite_benchmark ..\OSP-PBPK-Model-Library --write
python -m examples.build_metabolite_benchmark --hard ..\OSP-PBPK-Model-Library --write
python -m examples.build_biologic_benchmark   ..\OSP-PBPK-Model-Library --write
python -m examples.build_ddi_benchmark        ..\OSP-PBPK-Model-Library --write
python -m examples.benchmark_coverage         ..\OSP-PBPK-Model-Library
```

## 4. One LLM loop per task type (needs ANTHROPIC_API_KEY)

```
python -m examples.run_llm_build --snapshot <...blanked.json> --input <...input.json> --report r.html
python -m examples.run_llm_ddi   --snapshot <...ddi_blanked.json> --input <...ddi_input.json> --victim Midazolam
```

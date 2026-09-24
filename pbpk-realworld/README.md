# PBPK real-world projects

What a modeler is actually handed — a report and data files — instead of the
benchmark's pre-digested `input.json`. Generated from `OSP-PBPK-Model-Library/`
by `pkpd-agent/examples/build_realworld_projects.py`; the original benchmark is
left untouched.

Every compound has **both versions**:

```
<Compound>/
  test/                 the compound AS the target — what the agent is given
    report/<C>_inputs.md  report with the ANSWER redacted: keeps background,
                          physicochemical data, clinical study descriptions, and
                          the literature in-vitro kinetics; drops the modeling
                          strategy, the chosen methods, and every fitted value.
    data/*.csv            the building studies as raw time,concentration files
    README.md
  context/              the compound AS a reference for OTHER compounds
    report/<C>.md         the complete evaluation report (full disclosure)
    data/*.csv            every clinical study
  answer/               sealed — for grading only
    held_out.txt          the verification studies used to grade prediction
    reference.json        the reference model's methods + fitted parameters
```

The agent extracts the inputs from `report/` + `data/` itself — nothing is
pre-parsed. Leave-one-out is enforced at run time: when a compound is the target,
the harness may read every OTHER compound's `context/`, never the target's own.

Regenerate with:

```
cd pkpd-agent
python -m examples.build_realworld_projects --write
```

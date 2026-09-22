# PBPK agent scoreboard (agent vs ground-truth model)

PRIMARY grade is HELD-OUT: the agent fits the building studies and is graded on how well it PREDICTS the held-out verification studies, RELATIVE to the reference model on the same held-out data (PASS = agent held-out GMFE <= 1.25x reference's). Structure = did the agent pick the reference's methods/processes. 'all-data GMFE' is the in-sample fit, shown for context.

| # | model | held-out agent | held-out ref | ratio | structure | all-data agent | verdict |
|---|---|---|---|---|---|---|---|
| 1 | Montelukast | - | - | - | - | - | (no report - run failed / not run) |
| 2 | Propofol | - | - | - | - | - | (no report - run failed / not run) |
| 3 | Raltegravir | - | - | - | - | - | (no report - run failed / not run) |

**0/3 models produced a scored run; of those, 0 had a held-out split and 0 PASSED (predicted held-out data within 1.25x of the reference); 0 had too few studies to hold out (graded on all data).** Missing rows are runs that did not finish - report them as gaps, not passes.
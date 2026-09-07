# HARD mode — mechanism discovery (build the model, don't fill in the blank)

The normal benchmark hands the agent the reference **structure** (which enzyme
metabolizes the drug, which distribution/permeability methods) and asks it to
recover the fitted numbers — a parameter-identification task. HARD mode removes
the structure so the agent must BUILD the model, and it withholds the
metabolizing-enzyme **identity** so the agent must DISCOVER the mechanism from the
drug's chemistry/pharmacology and a candidate pool.

Files (per model, alongside the normal ones — built by `build_benchmark.py --hard`):

- `benchmark/<Model>.hard_blanked.json` — value-blanked AND structure-stripped:
  no elimination processes, partition/permeability at neutral PK-Sim Standard.
- `json_input/<Model>.hard.input.json` — same physchem + observed data, but the
  biology facts do NOT name the enzyme; instead a `candidate_clearance_molecules`
  pool (the distinct molecules the model expresses) is given. The agent reasons out
  which candidate clears the drug, ADDs the process, and chooses the methods. The
  answer key is shared, so grading is `did the agent rediscover the reference
  mechanism` (the report's structure/processes comparison, now EARNED not given).

Run:  `python -m examples.run_osp_scoreboard --hard --run --aggregate --pksim <cli>`
Verify no leak:  `python -m examples.verify_benchmark --type single --hard <ref>.json`

## Discovery difficulty (honest)

Discovery is only a real test when the candidate pool has **distractors** beyond the
answer. The pool is the reference model's distinct expressed molecules (it must be —
`add_processes` requires an expressed enzyme). Where the reference only expressed the
enzyme it used, pool == answer and the discovery is trivial. This is a property of
each source model, reported plainly rather than hidden:

| model | answer (withheld) | pool | distractors | genuine test? |
|---|---|---|---|---|
| Vancomycin-Pediatrics | — (no metabolism; pure renal) | 17 | 17 | yes |
| Propofol-Pediatrics | CYP2B6, UGT1A9 | 19 | 17 | yes |
| Sildenafil-Model | CYP2C19, CYP2C9, CYP3A4 | 14 | 11 | yes |
| Dapagliflozin-Model | Hepatic-CYP, UGT1A9, UGT2B7 | 10 | 7 | yes |
| Triazolam-Model | CYP3A4 | 7 | 6 | yes |
| Tizanidine-Model | CYP1A2 | 7 | 6 | yes |
| Mexiletine-Model | CYP1A2, CYP2D6 | 8 | 6 | yes |
| Alprazolam-Model | CYP3A4 | 7 | 6 | yes |
| Alfentanil-Model | CYP3A4 | 7 | 6 | yes |
| Midazolam-Model | CYP3A4, GABRG2, UGT1A4 | 7 | 4 | yes |
| Vancomycin | — (no metabolism; pure renal) | 0 | 0 | no (pool==answer) |
| Sufentanil-Pediatrics | CYP3A4 | 1 | 0 | no (pool==answer) |
| Sufentanil | CYP3A4 | 1 | 0 | no (pool==answer) |
| Raltegravir-Pediatrics | UGT1A1, UGT1A9 | 2 | 0 | no (pool==answer) |
| Raltegravir-Model | UGT1A1, UGT1A9 | 2 | 0 | no (pool==answer) |
| Propofol | CYP2B6, UGT1A9 | 2 | 0 | no (pool==answer) |
| Montelukast-Pediatrics | CYP2C8, CYP2C9, CYP3A4, CYP3A5 | 4 | 0 | no (pool==answer) |
| Montelukast | CYP2C8, CYP2C9, CYP3A4, CYP3A5 | 4 | 0 | no (pool==answer) |
| Digoxin-Model | ATP1A2, P-gp | 2 | 0 | no (pool==answer) |

**10/19 models are genuine discovery tests** (≥1 distractor enzyme).
The rest still exercise model-BUILDING (add the process, pick methods, fit) but not
enzyme discrimination. Note Vancomycin-Pediatrics (17 candidates, answer = *no*
metabolism) is the hardest kind: the agent must resist adding any enzyme and model it
as purely renal. To make a trivial case genuine, its source expression system would
need distractor enzymes added — a source-model change, not a generator one.

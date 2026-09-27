"""Single source of truth for the order benchmark tasks run in.

There is no per-run "pick which compounds" UI yet, so the DEFAULT order is what
decides which tasks a partial run (limited time / budget / steps) actually covers.
It should therefore front-load the tasks that most inform the scoreboard, so an
early stop still yields a representative read of the agent vs the OSP experts.

Ordering principle (an INFORMED default, not a measured truth - easy to override
with an explicit selection):

  1. KIND: adult whole-body PBPK first (the core task the benchmark measures),
     then pediatric variants, then DDI. DDI is last because the web runner does
     not execute DDI tasks yet (they need the interaction pipeline), so running
     them first would just burn slots on a not-yet-supported path.

  2. MECHANISM DIVERSITY within adult: front-load a spread of clearance
     mechanisms - a CYP3A4 substrate, a pure renal (GFR) drug, a P-gp
     transporter substrate, a CYP2D6 (polymorphic) drug, a UGT glucuronidation
     drug, a CYP1A2 drug - so the first handful of runs already exercises the
     breadth of PBPK behaviour rather than five variations of one mechanism.

  3. REFERENCE TRUSTWORTHINESS: a compound whose published model the harness does
     NOT reproduce faithfully pollutes the scoreboard, so it is deprioritised.
     Propofol's target-controlled-infusion arms are not reproduced (no
     infusion-protocol support yet), so it goes last among adults.

`ADULT_IMPORTANCE` lists compound directory names most-informative first. Anything
not listed (e.g. a newly added library model) sorts after the listed ones,
alphabetically - so additions degrade gracefully instead of breaking.
"""

# adult compounds, most-informative first (see the principle above)
ADULT_IMPORTANCE = [
    "Midazolam",      # canonical CYP3A4 probe substrate - the reference DDI victim
    "Vancomycin",     # pure renal (GFR) clearance - a mechanism with no enzyme at all
    "Digoxin",        # P-gp transporter substrate - distribution/transport, not metabolism
    "Mexiletine",     # CYP2D6 (polymorphic metaboliser)
    "Raltegravir",    # UGT1A1 glucuronidation
    "Sildenafil",     # CYP3A4 / CYP2C9, data-rich
    "Alfentanil",     # CYP3A4, IV - a clean, well-behaved baseline
    "Montelukast",    # CYP2C8 / CYP3A4
    "Dapagliflozin",  # UGT-mediated
    "Fluvoxamine",    # CYP1A2 / CYP2D6
    "Tizanidine",     # CYP1A2
    "Alprazolam",     # CYP3A4
    "Triazolam",      # CYP3A4
    "Sufentanil",     # CYP3A4
    "Propofol",       # LAST: TCI infusion arms not faithfully reproduced (untrustworthy ref)
]

# adult first, then pediatric, then DDI (see principle #1)
KIND_ORDER = {"adult": 0, "pediatric": 1, "ddi": 2}
_BIG = 10_000


def task_sort_key(dir_name: str, kind: str) -> tuple:
    """Sort key for one task: (kind rank, importance index, dir name). Unlisted
    compounds fall after listed ones (importance index _BIG), then alphabetical."""
    try:
        imp = ADULT_IMPORTANCE.index(dir_name)
    except ValueError:
        imp = _BIG
    return (KIND_ORDER.get(kind, _BIG), imp, dir_name or "")


def order_tasks(tasks: list) -> list:
    """Return ``tasks`` (dicts with 'dir' and 'kind') ordered by importance,
    most-informative first. Stable and non-mutating."""
    return sorted(tasks, key=lambda t: task_sort_key(t.get("dir", ""), t.get("kind", "")))

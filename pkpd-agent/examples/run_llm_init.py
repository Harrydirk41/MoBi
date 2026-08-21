r"""One-command project init for a traditional modeler: describe the model, get a config.

Turns a model dump + a few plain-language sentences into a ready projects/<name>/ folder
(spec.json + tasks.json), so the modeler never edits JSON or memorizes a schema:

  1. reads network.json (the model structure dump);
  2. the LLM extracts the structure (any naming) and drafts the task-role candidates;
  3. the LLM fills a tasks.json from the modeler's DESCRIPTION (prunes candidates, places
     the drug/dose names, clinical numbers and timeline the description states);
  4. validates the result against the real model and prints plain-English problems;
  5. writes projects/<name>/{spec.json,tasks.json} and lists what still needs the modeler.

    set ANTHROPIC_API_KEY=...
    python -m examples.run_llm_init --network network.json --name my_model ^
        --describe describe.txt          (or --describe "a few sentences...")

`describe.txt` is free text - e.g. "Psoriasis model, severity readout PASI, active band
6-20. Drugs: secukinumab (SEC_300mg dose), the IL-17 blocker. Match the UNCOVER-2 trial
PASI75 at week 12 = 77%. Disease drivers are the F_* amplification factors; calibrate
KD_SEC (reference 1e-10 M)."  The more you state, the less is left as a stub.

Needs ANTHROPIC_API_KEY.
"""

from __future__ import annotations

import argparse
import json
import os

from pkpd_agent.config import AgentConfig
from pkpd_agent.engines import llm_structure as LS
from pkpd_agent.engines import llm_tasks as LT
from pkpd_agent.engines import llm_config_build as LC
from pkpd_agent.engines import project_validate as PV


def _read_describe(arg: str) -> str:
    if arg and os.path.isfile(arg):
        with open(arg, encoding="utf-8") as fh:
            return fh.read()
    return arg or ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--name", required=True, help="project folder name (projects/<name>/)")
    ap.add_argument("--describe", default="",
                    help="plain-language description, or a path to a .txt file")
    ap.add_argument("--llm-model", default=None)
    ap.add_argument("--out-dir", default=None, help="override the projects dir")
    args = ap.parse_args()

    with open(args.network, encoding="utf-8") as fh:
        network = json.load(fh)
    description = _read_describe(args.describe)
    if not description.strip():
        print("[warn] no --describe given; the config will be mostly stubs to fill by hand.")

    cfg = AgentConfig(mock=False)
    if args.llm_model:
        cfg.model = args.llm_model
    if not cfg.anthropic_key_present():
        print("ANTHROPIC_API_KEY not set.")
        return

    call = LT.default_call(cfg)          # larger max_tokens + empty-reply guard

    print("== 1/4 extracting structure (any naming) ==", flush=True)
    structure = LS.extract_structure(network, call)
    print(f"   {len(structure['nodes'])} nodes, {len(structure['edges'])} edges")

    print("== 2/4 drafting task-role candidates ==", flush=True)
    tasks_draft = LT.draft_tasks(network, call, name=network.get("name", args.name))
    print(f"   {len(tasks_draft['vpop_drivers'])} vpop / "
          f"{len(tasks_draft['design_targets'])} design / "
          f"{len(tasks_draft['fit_params'])} fit candidates")

    print("== 3/4 filling config from your description ==", flush=True)
    tasks = LC.build_tasks(network, description, tasks_draft, call)
    spec = LC.build_spec(network, description, structure)

    print("== 4/4 validating against the model ==", flush=True)
    report = PV.validate_project(tasks, spec, network)
    print(PV.format_report(report))

    base = args.out_dir or os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "projects"))
    folder = os.path.join(base, args.name)
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "tasks.json"), "w", encoding="utf-8") as fh:
        json.dump(tasks, fh, indent=2)
    with open(os.path.join(folder, "spec.json"), "w", encoding="utf-8") as fh:
        json.dump(spec, fh, indent=2)
    print(f"\nwrote {folder}/tasks.json and spec.json")
    print("Next: skim the WARNINGS above, fill any stubbed clinical numbers / dose "
          "names, then run e.g.  python -m examples.run_llm_qsp_full --model "
          f"{args.name} --sbproj ... --vpop ...")


if __name__ == "__main__":
    main()

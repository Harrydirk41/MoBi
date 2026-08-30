r"""End-to-end assembly (steps 3+4+7) WITHOUT reading any paper: take the topology + a library
rate-law motif + the KNOWN parameter values, assemble a runnable SBML subsystem, simulate it,
and compare its steady state to the real model. Measures the ASSEMBLY - how close a from-scratch
model built by the standard library convention (with parameters assumed known) gets to the
hand-built one - isolating what the modellers' fine rate-law tuning adds on top.

    python -m examples.run_qsp_assemble --model ra ^
        --sbproj "..\RA-QSP-Model\Vantage RA QSP Model v1.0.sbproj" --cell FLS

The cell's proliferation regulators (Maxby knobs), their strengths, the baseline growth/death
rates, and the regulator cytokines' steady-state levels are all read FROM the real model
(parameters known); only the STRUCTURE + motif are assembled. Needs the MATLAB engine.
"""

from __future__ import annotations

import argparse
import os
import re
import tempfile

from pkpd_agent.engines.simbiology import SimBiologyEngine
from pkpd_agent.engines.qsp_model import QSPModel, get_spec
from pkpd_agent.engines import model_assembly as MA, llm_discover as D


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="ra")
    ap.add_argument("--sbproj", required=True)
    ap.add_argument("--cell", default="FLS")
    ap.add_argument("--readout-day", type=float, default=199.0)
    ap.add_argument("--fit-k", action="store_true",
                    help="calibrate each half-effect K from a per-cytokine isolating experiment "
                         "on the real model (closes the residual gap the guessed K leaves)")
    ap.add_argument("--llm-motif", action="store_true",
                    help="let the LLM choose the rate-law form (order + how regulators combine) "
                         "from the biology + the reference rate law - the agent picks the motif")
    ap.add_argument("--blind", action="store_true",
                    help="with --llm-motif, do NOT show the LLM the reference rate law "
                         "(test whether it infers the form from biology alone)")
    ap.add_argument("--llm-model", default=None)
    args = ap.parse_args()

    sb = SimBiologyEngine()
    try:
        print("== starting MATLAB =="); sb.start()
        print(f"== loading {os.path.basename(args.sbproj)} =="); sb.load_project(args.sbproj)
        net = sb.network_json(os.path.join(os.getcwd(), "network.json"))
        nset = set(QSPModel(net, get_spec(args.model)).nodes)
        params = sb.list_parameters().get("parameters", [])
        pv = {p["name"]: float(p["value"]) for p in params}

        # discover the cell's baseline growth/death rates and its Maxby regulators
        base_p = next((n for n in pv if re.search(rf"(?i)k[gp]_{args.cell}.*base", n)), None)
        apop_p = next((n for n in pv if re.search(rf"(?i)kd_{args.cell}.*base", n)), None)
        medges = D.maxby_edges(list(pv), nset)
        regs = [(src, dst, knob) for (src, dst), knob in medges.items()
                if dst == args.cell and args.cell + "Prolif" in knob]
        if not (base_p and apop_p and regs):
            print(f"could not resolve FLS growth/death/regulators: base={base_p} apop={apop_p} "
                  f"regs={[r[2] for r in regs]}"); return
        print(f"assembling {args.cell}: growth {base_p}={pv[base_p]:g}, death {apop_p}="
              f"{pv[apop_p]:g}, {len(regs)} proliferation regulators")

        # show the REAL rate law + effect rule governing this cell's proliferation, so the
        # combination rule (how the regulators actually combine) is visible - our library motif
        # is only a guess for it.
        ref_rate = ""
        for rx in net.get("reactions", []):
            if args.cell in (rx.get("products") or []) and "prolif" in (rx.get("rate", "")
                                                                        + str(rx.get("reaction"))).lower():
                print(f"  REAL proliferation rate law: {rx.get('rate')}")
                ref_rate += (rx.get("rate") or "") + "\n"
        for ru in net.get("rules", []):
            expr = ru.get("rule", "") if isinstance(ru, dict) else str(ru)
            if re.search(rf"(?i){args.cell}.*prolif.*effect\s*=", expr):
                print(f"  REAL effect rule: {expr}")
                ref_rate += expr + "\n"

        # real model's steady-state cytokine levels (the clamp) and the real cell value (truth)
        prof = {k: v[-1] for k, v in sb.simulate(stop_time=args.readout_day + 1.0)
                .get("columns", {}).items() if v}
        real_cell = prof.get(args.cell)

        # build the regulator list + parameter values; K (half-effect) defaults to the cytokine's
        # own steady-state level (Hill at 0.5) - the library default, since a per-edge K is not
        # separately exposed. This is the standard-motif approximation the comparison measures.
        kg, kd = pv[base_p], pv[apop_p]

        # optional: calibrate each K from a per-cytokine ISOLATING experiment on the real model.
        # 5 K's cannot be pinned from one FLS target (under-determined), so isolate: set the
        # OTHER regulators' Maxby to 1.0 (no effect), read the real FLS driven by only this
        # cytokine -> its fold, then solve K analytically for the motif to reproduce it.
        fittedK = {}
        if args.fit_k:
            print("== calibrating each K from a per-cytokine isolating experiment ==", flush=True)
            for src, dst, knob in regs:
                saved = {k2: sb.set_parameter(k2, 1.0) for _, _, k2 in regs if k2 != knob}
                only = {k: v[-1] for k, v in sb.simulate(stop_time=args.readout_day + 1.0)
                        .get("columns", {}).items() if v}
                for k2, old in saved.items():
                    sb.set_parameter(k2, old)              # restore
                X, Max = max(prof.get(src, 0.0), 1e-12), pv[knob]
                fold = (only.get(args.cell, 0.0) * kd / kg) if kg else 1.0   # real fold from src
                frac = (fold - 1.0) / (Max - 1.0) if Max != 1 else 0.0        # = X/(K+X)
                K = X * (1.0 / frac - 1.0) if 0 < frac < 1 else (1e-9 if frac >= 1 else 1e12)
                fittedK[src] = max(K, 1e-9)
                print(f"    {src}: real fold {fold:.3f} -> K = {fittedK[src]:.3g}")

        regulators, values, clamp = [], {base_p: kg, apop_p: kd}, {}
        for src, dst, knob in regs:
            k_name = f"K_{dst}_{src}"
            regulators.append({"species": src, "max_param": knob, "k_param": k_name})
            values[knob] = pv[knob]
            values[k_name] = fittedK.get(src, max(prof.get(src, 1.0), 1e-9))  # fitted, or baseline
            clamp[src] = prof.get(src, 0.0)

        motif = None
        if args.llm_motif:
            from pkpd_agent.config import AgentConfig
            from pkpd_agent.engines import llm_tasks as LT
            cfg = AgentConfig(mock=False)
            if args.llm_model:
                cfg.model = args.llm_model
            if not cfg.anthropic_key_present():
                print("ANTHROPIC_API_KEY not set for --llm-motif."); return
            print("== LLM choosing the rate-law form (motif) "
                  + ("from biology alone ==" if args.blind else "from biology + reference ==") )
            motif = MA.propose_motif(args.cell, regulators,
                                     "" if args.blind else ref_rate, LT.default_call(cfg))
            print(f"  LLM motif: order={motif['proliferation_order']}, "
                  f"combine={motif['combination']}, cap={motif.get('cap')}"
                  + (f"  ({motif['reason']})" if motif.get('reason') else ""))

        spec = MA.build_subsystem(args.cell, base_p, apop_p, regulators, values, clamp,
                                  motif=motif)
        sbml = os.path.join(tempfile.gettempdir(), f"{args.cell}_assembled.sbml")
        with open(sbml, "w", encoding="utf-8") as fh:
            fh.write(MA.to_sbml(spec))
        print(f"  emitted SBML: {len(spec['species'])} species, {len(spec['parameters'])} "
              f"parameters, {len(spec['reactions'])} reactions -> {sbml}")

        asm_cell = sb.import_simulate(sbml, args.cell, stop_time=args.readout_day + 1.0)
        print(f"\n== assembled subsystem vs real model ==")
        print(f"  real model     {args.cell} steady state = {real_cell:g}")
        print(f"  assembled model {args.cell} steady state = {asm_cell:g}")
        if real_cell:
            r = asm_cell / real_cell
            print(f"  ratio assembled/real = {r:.2f}  "
                  + ("(close - the standard motif + known params reproduce it)" if 0.5 < r < 2
                     else "(off - the real rate-law form / half-effects are hand-tuned beyond "
                          "the library motif; this gap is what calibration adds)"))
        print("\n  NOTE: topology + motif + wiring were assembled; parameters were taken as "
              "known\n  (not read from papers). The residual gap isolates the rate-law-form / "
              "half-effect\n  tuning that a standard library motif does not capture.")
    finally:
        sb.stop()


if __name__ == "__main__":
    main()

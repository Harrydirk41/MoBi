"""PBPK Copilot backend — a small FastAPI server that drives the real agent.

Endpoints:
  GET  /                      -> the copilot web UI
  GET  /api/models            -> compounds with a benchmark snapshot + status
  GET  /api/library?compound  -> the leave-one-out reference library
  POST /api/run               -> start an agent build (returns {run_id})
  GET  /api/stream/{run_id}   -> Server-Sent Events: the live agent loop

The agent loop runs in a background thread; its Decision/Observation events and
the optimizer's per-eval stdout are pushed to a per-run queue and streamed to
the browser as SSE. One run at a time (stdout capture is process-global).
"""
from __future__ import annotations

import csv
import glob
import io
import json
import os
import queue
import threading
import uuid
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.dirname(_HERE)                       # pkpd-agent/
_LIB = os.path.abspath(os.path.join(_PKG, "..", "OSP-PBPK-Model-Library"))
_RW = os.path.abspath(os.path.join(_PKG, "..", "pbpk-realworld"))   # raw report+data tree


def _rw_projects() -> list[str]:
    if not os.path.isdir(_RW):
        return []
    return sorted(d for d in os.listdir(_RW)
                  if os.path.isdir(os.path.join(_RW, d, "test")))


def _rw_series(project: str, kind: str = "test", cap: int = 400) -> list[dict]:
    """Parse a project's data CSVs into plottable series (downsampled)."""
    ddir = os.path.join(_RW, project, kind, "data")
    man = os.path.join(ddir, "_manifest.json")
    entries = json.load(open(man, encoding="utf-8")) if os.path.isfile(man) else \
        [{"file": os.path.basename(f)} for f in glob.glob(os.path.join(ddir, "*.csv"))]
    out = []
    for e in entries[:24]:
        p = os.path.join(ddir, e["file"])
        if not os.path.isfile(p):
            continue
        pts = []
        for i, line in enumerate(open(p, encoding="utf-8")):
            if i == 0 or not line.strip():
                continue
            try:
                t, c = line.split(",")[:2]
                pts.append([float(t), float(c)])
            except ValueError:
                pass
        if len(pts) >= 2:
            out.append({"study": e.get("study") or e["file"], "route": e.get("route", ""),
                        "dose": e.get("dose", ""), "points": pts[:cap]})
    return out

# ---- runs registry -------------------------------------------------------- #
_RUNS: dict[str, dict] = {}


class _LoggingQueue(queue.Queue):
    """A recording BROADCAST queue for one run.

    It does two jobs so a run can be watched from several pages at once (e.g. a
    laptop and a second machine) and reloaded later:
      * records every dict event to ``.log`` (persisted to disk for replay), and
      * fans each event out to any number of live subscriber queues, so multiple
        SSE connections each get their OWN copy instead of stealing events from a
        single shared queue.
    A late joiner calls ``subscribe()`` to get the backlog-so-far plus a live
    queue, taken together under the lock so no event is dropped or duplicated
    across the hand-off."""
    def __init__(self):
        super().__init__()
        self.log: list = []
        self._subs: list[queue.Queue] = []
        self._lock = threading.Lock()

    def put(self, item, *a, **k):
        if isinstance(item, dict):
            with self._lock:
                self.log.append(item)
                subs = list(self._subs)
            for s in subs:
                s.put(item)
        return super().put(item, *a, **k)

    def subscribe(self):
        """Register a live subscriber. Returns (backlog_snapshot, sub_queue); the
        snapshot and the registration happen under one lock so the subscriber's
        queue receives exactly the events AFTER the snapshot."""
        s: queue.Queue = queue.Queue()
        with self._lock:
            backlog = list(self.log)
            self._subs.append(s)
        return backlog, s

    def unsubscribe(self, s) -> None:
        with self._lock:
            if s in self._subs:
                self._subs.remove(s)


def _run_record_path(compound: str) -> str:
    # resolve the base dir WITHOUT _find/_task_dir (which would recurse through
    # _compounds -> _load_run -> here); a variant cid is "<Dir> · pediatric/DDI".
    comp_dir = str(compound).split(" · ")[0]
    slug = "".join(ch if ch.isalnum() else "_" for ch in str(compound))
    return os.path.join(_LIB, comp_dir, "report", f"copilot_run__{slug}.json")


def _persist_run(compound: str, events: list, p: dict) -> None:
    """Write a run's outcome + full trace to <compound>/report/copilot_run.json so
    the page can reload/resume it later. Best-effort — never sinks a run."""
    try:
        done = next((e for e in reversed(events) if e.get("type") == "done"), {})
        # the pipeline report (streamed before 'end') carries the held-out grade
        rep = next((e for e in reversed(events) if e.get("type") == "report"), {})
        held = ((rep.get("held_out") or {}).get("gmfe")) if rep else None
        rec = {
            "compound": compound, "ts": __import__("time").time(),
            "model": p.get("model"), "hard": bool(p.get("hard")),
            "web": bool(p.get("web")), "self_extract": bool(p.get("self_extract")),
            "fit_vision": bool(p.get("fit_vision", True)),
            "minimal": bool(p.get("minimal")),
            "in_sample_gmfe": done.get("best_gmfe"),
            "held_out_gmfe": held,               # from the pipeline report; /api/runs/<c>/grade can still set it
            "best_edits": done.get("best_edits"),
            "givens": done.get("givens"),
            "web_lookups": done.get("web_lookups") or [],
            "blind": done.get("blind", True),
            "events": events,
        }
        path = _run_record_path(compound)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(rec, fh)
    except Exception:                            # noqa: BLE001
        pass


def _load_run(compound: str) -> "dict | None":
    path = _run_record_path(compound)
    if not os.path.isfile(path):
        return None
    try:
        return json.load(open(path, encoding="utf-8"))
    except (OSError, ValueError):
        return None
# Concurrent agent runs. Default 1 (serial, unchanged). Set PKPD_COPILOT_JOBS>1
# to run several compounds at once — each run's stdout is routed per-thread so
# their SSE streams don't interleave. Keep near the core count.
_RUN_SEM = threading.Semaphore(max(1, int(os.environ.get("PKPD_COPILOT_JOBS", "1"))))
_SIM_LOCK = threading.Lock()                        # serialize reference-fit sims


class _StdoutRouter:
    """A process-wide sys.stdout replacement that routes writes to the writer
    registered for the CURRENT thread (falling back to the real stdout). Lets
    parallel runs each capture their own output without a global redirect."""
    def __init__(self, real):
        self._real = real
        self._map: dict = {}

    def register(self, w):
        self._map[threading.get_ident()] = w

    def unregister(self):
        self._map.pop(threading.get_ident(), None)

    def _w(self):
        return self._map.get(threading.get_ident(), self._real)

    def write(self, s):
        return self._w().write(s)

    def flush(self):
        try:
            self._w().flush()
        except Exception:                               # noqa: BLE001
            pass

    def __getattr__(self, n):
        return getattr(self._real, n)


_ROUTER = None


def _install_router():
    global _ROUTER
    if _ROUTER is None:
        import sys
        _ROUTER = _StdoutRouter(sys.stdout)
        sys.stdout = _ROUTER
    return _ROUTER
_REFFIT_CACHE: dict[str, list] = {}                 # project -> per-dataset obs+sim


def _downsample(pts: list, n: int = 160) -> list:
    if len(pts) <= n:
        return pts
    step = len(pts) / n
    return [pts[int(i * step)] for i in range(n)]


_PART_DESC = {
    "PK-Sim Standard": "Lipophilicity-based empirical tissue:plasma partitioning (PK-Sim default).",
    "Rodgers and Rowland": "Mechanistic tissue:plasma for ionizable drugs — tissue lipids/water plus "
    "electrostatic binding to acidic phospholipids (bases) or to albumin (acids/neutrals).",
    "Schmitt": "Mechanistic partitioning from lipophilicity, pKa and tissue lipid/water/protein content.",
    "Poulin and Theil": "Lipid/water partition from logP and tissue composition (neutral-drug oriented).",
    "Berezhkovskiy": "Poulin & Theil corrected for the drug fraction in tissue vs plasma (volume correction).",
}
_PERM_DESC = {
    "PK-Sim Standard": "Cellular permeability from lipophilicity / molecular weight (PK-Sim default).",
    "Charge dependent Schmitt": "Permeability accounting for the charged fraction (pKa-dependent).",
    "Charge dependent Schmitt normalized to PK-Sim": "Charge-dependent Schmitt rescaled to the "
    "PK-Sim standard baseline.",
}


def _catalog() -> dict:
    """The whole (compound-independent) PK-Sim action space the agent chooses from:
    distribution/permeability methods, the mechanisms it may add, and the universe
    of parameters grouped by tier (never-fit / measured / fittable)."""
    from pkpd_agent.engines import osp_catalog as C
    from pkpd_agent.engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
    procs = []
    for key, spec in C.PROCESS_TYPES.items():
        procs.append({"key": key, "description": spec.get("description"),
                      "applies_to": spec.get("applies_to"),
                      "validated": spec.get("validated", False),
                      "parameters": [p.get("name") for p in spec.get("parameters") or []]})
    ddi = [{"type": d.get("type"), "description": d.get("description"),
            "parameters": [p.get("name") for p in d.get("parameters") or []]}
           for d in C.interaction_process_types()]
    params = {"constant": [], "measured_soft": [], "estimate": []}
    for name, v in C.PARAM_CATALOG.items():
        tier = v.get("tier", "estimate")
        params.setdefault(tier, []).append(
            {"name": name, "role": v.get("role"), "description": v.get("description")})
    # the FULL PK-Sim / MoBi capability map, each area tagged with how much of it the
    # agent actually wires up: active | partial | pksim (native, not wired) | mobi.
    pksim_areas = [
        {"area": "Administration routes", "status": "partial",
         "note": "IV & oral are active; the others are PK-Sim-native but not wired into the agent.",
         "items": ["IV bolus", "IV infusion", "Oral", "Dermal", "Inhalation",
                   "Intramuscular", "Subcutaneous", "Intra-arterial", "Intraperitoneal", "Ocular"]},
        {"area": "Formulations / dissolution", "status": "partial",
         "note": "Weibull dissolution parameters are fit; other dissolution/particle models are not exposed.",
         "items": ["Dissolved (solution)", "Weibull", "Lint80 particle dissolution",
                   "Tablet", "Capsule", "Suspension"]},
        {"area": "Populations", "status": "pksim",
         "note": "Pediatric snapshots exist; building virtual populations is not wired.",
         "items": ["Virtual population", "Age / weight / height / BMI", "Ethnicity database",
                   "Disease states", "Pregnancy", "Pediatric ontogeny", "Geriatric"]},
        {"area": "Species & physiology", "status": "pksim",
         "note": "Human physiology is fixed; editing it or switching species is not exposed.",
         "items": ["Human", "Rat", "Dog", "Monkey", "Mouse", "Minipig", "Rabbit",
                   "Organ volumes", "Blood flows", "Tissue composition",
                   "GFR / hematocrit / cardiac output"]},
        {"area": "Expression & ontogeny", "status": "pksim",
         "note": "Per-model enzyme/transporter expression is used; the full database & ontogeny are not visualized.",
         "items": ["Relative tissue expression", "Enzyme/transporter database",
                   "Age-dependent ontogeny", "Reference concentrations"]},
        {"area": "Large molecules / biologics", "status": "pksim",
         "note": "Small-molecule disposition only; biologic mechanisms are not modeled.",
         "items": ["FcRn recycling", "Target-mediated drug disposition (TMDD)", "Target binding",
                   "Lymph flow", "Size-based (2-pore) distribution", "Endosomal turnover"]},
        {"area": "MoBi (open modeling)", "status": "mobi",
         "note": "The open-ended engine behind PK-Sim — entirely outside the structured action space.",
         "items": ["Custom ODEs", "Arbitrary reaction networks", "Events",
                   "Custom compartments & molecules", "Observers", "Passive transports"]},
        {"area": "Simulation & analysis", "status": "pksim",
         "note": "Abstracted or automated by the harness.",
         "items": ["Solver settings", "Output intervals", "Parameter-identification algorithms",
                   "PK / NCA analysis", "Sensitivity analysis", "Population PK"]},
    ]
    return {
        "partition_methods": [{"name": m, "description": _PART_DESC.get(m)} for m in PARTITION_METHODS],
        "permeability_methods": [{"name": m, "description": _PERM_DESC.get(m)} for m in PERMEABILITY_METHODS],
        "process_types": procs,
        "ddi_types": ddi,
        "parameters": params,
        "pksim_areas": pksim_areas,
        "counts": {"partition": len(PARTITION_METHODS), "permeability": len(PERMEABILITY_METHODS),
                   "process_types": len(procs), "ddi": len(ddi),
                   "constant": len(params["constant"]), "measured": len(params["measured_soft"]),
                   "fittable": len(params["estimate"]),
                   "areas": len(pksim_areas)},
    }


def _model_view(f: dict) -> dict:
    """The selected compound's EDITABLE action space: current methods + all legal
    options, the parameters (with tier/range so the UI knows what may be fit), the
    active processes, and the molecules a mechanism can attach to."""
    from pkpd_agent.tools import osp_loop_tools as T
    from pkpd_agent.engines import osp_catalog as C
    from pkpd_agent.engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
    snap = f["snapshot"]
    model = T._current_model(snap)
    expressed = T._expressed_molecules(snap)                 # [{molecule,type}]
    inp = json.load(open(f["input"], encoding="utf-8"))
    lit = (inp.get("given_data", {}) or {}).get("literature_physicochemical", []) or []
    cms = model["calculation_methods"]

    def cur(word):
        for m in cms:
            if word in m.lower():
                return m.split(" - ")[-1].strip()
        return None

    editable = []
    for p in model["parameters"]:
        base = p["name"].split("@", 1)[0]
        cat = C.describe_parameter(base)
        tier = C.param_tier(base)
        e = {"name": p["name"], "value": p.get("value"), "unit": p.get("unit"),
             "value_status": p.get("value_status"), "tier": tier,
             "role": cat.get("role"), "range": cat.get("range"),
             "description": cat.get("description")}
        if tier == "measured_soft":
            mr = C.measured_range(p["name"], lit)
            if mr:
                e["measured_range"] = list(mr)
        editable.append(e)
    return {
        "compound": f["compound"],
        "methods": {
            "partition": {"current": cur("partition"), "options": list(PARTITION_METHODS)},
            "permeability": {"current": cur("permeability"), "options": list(PERMEABILITY_METHODS)}},
        "parameters": editable,
        "processes": model["processes"],
        "expressed": expressed,
        "addable": [{"type": p.get("type"), "can_attach_to": p.get("can_attach_to")}
                    for p in C.addable_process_types(expressed)],
    }


def _try_model(f: dict, edits: dict, subset: str = "building") -> dict:
    """Run the human-edited model. If `edits` has an `estimate` block, optimize it;
    otherwise just build+run with the given values. `subset` picks which studies to
    score on: 'building' (fit set) or 'held_out' (blind grade). Returns GMFE +
    per-study observed vs simulated series for the fit overlay. Needs PKSim.CLI."""
    from pkpd_agent.config import AgentConfig
    from pkpd_agent.engines.osp_cli import OSPCli
    from pkpd_agent.engines import osp_score, osp_split
    from pkpd_agent.engines import osp_optimize as OO
    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=cfg.pksim_cli_path, timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        return {"ok": False, "error": "PKSim.CLI not found — set PKPD_PKSIM_CLI to run a model."}
    inp = json.load(open(f["input"], encoding="utf-8"))
    observed = (inp.get("given_data", {}) or {}).get("clinical_observed_data")
    if not observed:
        return {"ok": False, "error": "This task has no plasma concentration data to fit "
                "(e.g. a DDI task) — it runs via the interaction pipeline, not the web runner yet."}
    snapd = json.load(open(f["snapshot"], encoding="utf-8"))
    split = osp_split.split_studies(snapd, observed)
    build = set(split.get("building") or [])
    if subset == "held_out":
        held = set(split.get("verification") or [])
        obs = [o for o in observed if o["dataset"] in held]
        if not obs:
            return {"ok": False, "error": "no held-out studies for this compound"}
    else:
        obs = [o for o in observed if o["dataset"] in build] or observed
    estimate = edits.get("estimate") or {}
    structure = {k: edits[k] for k in ("calculation_methods", "processes", "add_processes")
                 if k in edits}
    fixp = edits.get("parameters") or edits.get("fix") or {}
    optimized = None
    with _SIM_LOCK:
        if estimate:
            r = OO.run_optimization(cli, f["snapshot"], obs, estimate=estimate, fix=fixp,
                                    structure=structure)
            if not r.get("ok"):
                return {"ok": False, "error": r.get("message", "optimization failed")}
            optimized = r.get("optimized") or {}
            run_edits = {**structure, "parameters": {**fixp, **optimized}}
            gmfe = (r.get("fit") or {}).get("gmfe")
            by_route = r.get("by_route")
        else:
            run_edits = {**structure, "parameters": fixp}
            gmfe, by_route = None, None
        res = cli.build_and_run(f["snapshot"], edits=run_edits)
        if not res.get("ok"):
            return {"ok": False, "error": res.get("message", "PK-Sim run failed")}
        linkage = osp_score.linkage_from_snapshot(snapd)
        predicted, _ = osp_score.map_predictions(res["profiles"], obs, linkage)
        sc = osp_score.score_fit(obs, predicted)
        if gmfe is None:
            gmfe = sc["overall"].get("gmfe")
        by_route = by_route or sc.get("by_route")
        pk_gmfe = (sc.get("pk_parameters") or {}).get("gmfe_by_metric")
    pred_by = {p["dataset"]: p for p in predicted}
    series = []
    for o in obs:
        p = pred_by.get(o["dataset"])
        series.append({
            "study": o.get("study") or o["dataset"], "route": o.get("route", ""),
            "dose": o.get("dose", ""),
            "observed": [[t, c] for t, c in zip(o.get("time_h", []), o.get("conc_mg_L", []))
                         if c is not None],
            "simulated": _downsample([[t, c] for t, c in zip(p["time_h"], p["pred_conc_mg_L"])
                                      if c is not None]) if p else []})
    return {"ok": True, "gmfe": gmfe, "by_route": by_route, "pk_gmfe": pk_gmfe,
            "optimized": optimized, "series": series}


def _short_method(s: str) -> str:
    """'Cellular partition coefficient method - Rodgers and Rowland' -> 'Rodgers and Rowland'."""
    return str(s).split(" - ")[-1].strip() if s else s


def _reference_answer(f: dict) -> "dict | None":
    """The reference (ground-truth) structure + fitted values from the answer key.
    POST-HOC only — used to build the report's comparison, never handed to the agent."""
    import glob as _glob
    aks = [a for a in _glob.glob(os.path.join(_LIB, f["dir"], "answer_key", "*.answer_key.json"))
           if "ddi" not in os.path.basename(a)]
    if not aks:
        return None
    try:
        ak = json.load(open(aks[0], encoding="utf-8"))
    except (OSError, ValueError):
        return None
    sc = ak.get("structural_choices", {}) or {}
    methods = sc.get("calculation_methods") or []

    def m_of(kind):
        for s in methods:
            if kind in s.lower():
                return _short_method(s)
        return None
    procs = sc.get("metabolizing_processes") or sc.get("processes") or []
    return {
        "partition": m_of("partition"),
        "permeability": m_of("permeab"),
        "molecules": sorted({p.get("molecule") for p in procs if isinstance(p, dict) and p.get("molecule")}),
        "parameters": {e.get("parameter"): e.get("value") for e in ak.get("estimated_parameters") or []},
    }


def _meq(a, b) -> bool:
    return bool(a) and bool(b) and _short_method(a).lower() == _short_method(b).lower()


def _fold_verdict(fold):
    if fold is None:
        return "—"
    mag = fold if fold >= 1 else 1.0 / fold
    return ("recovered" if mag <= 1.25 else "close" if mag <= 2
            else "off" if mag <= 3 else "far")


def _report(f: dict, edits: dict, model: str | None = None) -> dict:
    """The post-build REPORT: in-sample + held-out fit curves and grades, the
    reference overlay, and a post-hoc comparison of the adopted model to the
    reference answer (structure + each fitted parameter). The comparison is a
    REVEAL for the human — it is never fed back to the agent."""
    edits = {k: v for k, v in (edits or {}).items() if k != "estimate"}   # forward run at the fitted values
    ins = _try_model(f, edits, subset="building")
    held = _try_model(f, edits, subset="held_out")
    try:
        ref = _reffit(f["dir"])
    except Exception:                                    # noqa: BLE001
        ref = {"ok": False}
    ans = _reference_answer(f)
    comparison = None
    if ans:
        cm = edits.get("calculation_methods") or {}
        amols = {p.get("molecule") for p in (edits.get("add_processes") or [])
                 if isinstance(p, dict) and p.get("molecule")}
        struct = [
            {"aspect": "partition method", "agent": cm.get("partition"),
             "reference": ans["partition"], "match": _meq(cm.get("partition"), ans["partition"])},
            {"aspect": "permeability method", "agent": cm.get("permeability"),
             "reference": ans["permeability"], "match": _meq(cm.get("permeability"), ans["permeability"])},
            {"aspect": "clearing molecules", "agent": sorted(amols),
             "reference": ans["molecules"], "match": (set(amols) == set(ans["molecules"]))},
        ]
        # The agent's model VALUE for a parameter is its fitted value OR the
        # literature value it FIXED — both are part of the adopted model. Comparing
        # only the optimized subset makes a correctly-fixed given (e.g. fu taken
        # from literature) show a misleading "—". And match by BASE name so a
        # molecule-qualified key ('Intrinsic clearance@CYP3A4') lines up with the
        # reference's plain 'Intrinsic clearance' instead of dropping to "—".
        def _basekey(n):
            return str(n).split("@", 1)[0].split("|", 1)[0].strip().lower()
        agent_vals: dict = {}
        for src in (edits.get("fix") or {}, edits.get("parameters") or {}):
            if isinstance(src, dict):
                for k, v in src.items():
                    agent_vals[_basekey(k)] = v      # a fitted value overrides a fixed one
        prows = []
        for name, rv in ans["parameters"].items():
            av = agent_vals.get(_basekey(name))
            fold = (av / rv) if (isinstance(av, (int, float)) and isinstance(rv, (int, float)) and rv) else None
            prows.append({"name": name, "agent": av, "reference": rv,
                          "fold": (round(fold, 3) if fold else None),
                          "verdict": _fold_verdict(fold)})
        comparison = {"structure": struct, "parameters": prows,
                      "reference_gmfe": None}
    def pack(r):
        return ({"gmfe": r.get("gmfe"), "pk_gmfe": r.get("pk_gmfe"), "series": r.get("series")}
                if r.get("ok") else {"error": r.get("error")})
    out = {"ok": True, "in_sample": pack(ins), "held_out": pack(held),
           "reference": ({"series": ref.get("series")} if ref.get("ok") else {"error": ref.get("error")}),
           "comparison": comparison}
    out["interpretation"] = _report_narrative(f["dir"], out, model)
    return out


def _report_narrative(compound: str, rep: dict, model: str | None) -> "str | None":
    """A short LLM interpretation of the finished model vs the reference — a
    post-hoc reviewer's read for the human (it sees the comparison, which is fine:
    it is never fed back into the modeling loop). Best-effort: no key -> None."""
    try:
        import anthropic
    except ImportError:
        return None
    if not os.environ.get("ANTHROPIC_API_KEY"):
        return None
    comp = rep.get("comparison") or {}
    facts = {
        "compound": compound,
        "in_sample_gmfe": (rep.get("in_sample") or {}).get("gmfe"),
        "held_out_gmfe": (rep.get("held_out") or {}).get("gmfe"),
        "pk_parameter_folds": (rep.get("held_out") or {}).get("pk_gmfe"),
        "structure_vs_reference": comp.get("structure"),
        "parameters_vs_reference": [p for p in (comp.get("parameters") or []) if p.get("reference") is not None],
    }
    prompt = (
        "You are a senior PBPK modeler reviewing a finished model AGAINST the "
        "reference answer (post-hoc review — the builder never saw this). Here are "
        "the facts as JSON:\n\n" + json.dumps(facts, default=str, indent=1) + "\n\n"
        "In 3–5 plain sentences, interpret honestly: did the STRUCTURE match the "
        "reference (methods, clearing molecules)? Is the HELD-OUT prediction good "
        "(GMFE ~<=1.5 is good, ~2 acceptable, >3 poor)? Which fitted PARAMETERS are "
        "far from the reference, and does it matter — a parameter far from the "
        "reference that the data cannot constrain (or that trades off with another) "
        "is acceptable, one that is far AND influential is a real miss. End with a "
        "one-line verdict. Do not restate the numbers as a table; interpret them.")
    try:
        client = anthropic.Anthropic()
        resp = client.messages.create(
            model=model or "claude-sonnet-5", max_tokens=600,
            messages=[{"role": "user", "content": prompt}])
        return "".join(b.text for b in resp.content if getattr(b, "type", None) == "text").strip() or None
    except Exception:                                # noqa: BLE001 - narrative must never sink the report
        return None


def _postmortem(compound: str, rep: dict, edits: dict, model: str | None) -> "str | None":
    """POST-STOP self-critique. After the blind build, the reference answer is
    REVEALED to the agent (the same model that built it) and it summarizes, in the
    first person, what it got wrong and why. Strictly post-hoc: this runs after the
    graded run has finished and never re-enters a scored build, so it cannot leak.
    Best-effort: no key -> None."""
    try:
        import anthropic
    except ImportError:
        return None
    if not os.environ.get("ANTHROPIC_API_KEY"):
        return None
    comp = rep.get("comparison") or {}
    edits = edits or {}
    facts = {
        "your_model": {
            "structure": {k: edits.get(k) for k in
                          ("calculation_methods", "add_processes", "processes") if edits.get(k)},
            "estimated_parameters": edits.get("parameters"),
            "fixed_parameters": edits.get("fix")},
        "your_in_sample_gmfe": (rep.get("in_sample") or {}).get("gmfe"),
        "your_held_out_gmfe": (rep.get("held_out") or {}).get("gmfe"),
        "structure_vs_reference": comp.get("structure"),
        "parameters_vs_reference": [p for p in (comp.get("parameters") or [])
                                    if p.get("reference") is not None],
    }
    prompt = (
        "You just built this PBPK model BLIND — you never saw the reference answer "
        "during the build. Now, post-hoc, the reference is REVEALED. Here is your "
        "model, your held-out grade, and the reference, as JSON:\n\n"
        + json.dumps(facts, default=str, indent=1) + "\n\n"
        "Write a short FIRST-PERSON post-mortem (3-6 sentences): where did your "
        "STRUCTURE or PARAMETERS differ from the reference; which differences "
        "actually mattered for the fit and which were harmless (a fitted value far "
        "from the reference that the data cannot constrain, or that trades off with "
        "another, is NOT a real miss — say so); and honestly, was the biological-"
        "reasonableness judgement you made during the build borne out or not? End "
        "with the ONE thing you would do differently next time. Do not restate the "
        "numbers as a table — interpret them.")
    try:
        client = anthropic.Anthropic()
        resp = client.messages.create(
            model=model or "claude-sonnet-5", max_tokens=700,
            messages=[{"role": "user", "content": prompt}])
        return "".join(b.text for b in resp.content if getattr(b, "type", None) == "text").strip() or None
    except Exception:                                # noqa: BLE001 - never sinks the run
        return None


def _built_view(project: str) -> dict:
    """A finished model's adopted methods + its key tunable parameters (with their
    fitted values and a plausible slider range), so the UI can perturb them and
    re-run to watch the curve move."""
    from pkpd_agent.engines.snapshot_edit import PARTITION_METHODS, PERMEABILITY_METHODS
    from pkpd_agent.engines import osp_score
    project = _task_dir(project)
    snap = os.path.join(_LIB, project, "json", f"{project}-Model.json")
    if not os.path.isfile(snap):
        return {"error": "no built model for this compound"}
    try:
        snapd = json.load(open(snap, encoding="utf-8"))
    except (OSError, ValueError):
        return {"error": "unreadable model"}
    comp = (snapd.get("Compounds") or [{}])[0]
    cms = comp.get("CalculationMethods") or []

    def cur(w):
        for m in cms:
            if isinstance(m, str) and w in m.lower():
                return m.split(" - ")[-1].strip()
        return None

    def first_val(block):
        for e in (block or []):
            for pp in (e.get("Parameters") or []):
                return pp.get("Value")
        return None

    params = []

    def add(name, val):
        if val is None:
            return
        try:
            v = float(val)
        except (TypeError, ValueError):
            return
        pb = None
        try:
            pb = osp_score.physical_bounds(name.split("@", 1)[0])
        except Exception:                               # noqa: BLE001
            pb = None
        if pb and pb[0] < pb[1]:
            lo, hi = pb
        elif v > 0:
            lo, hi = v / 4.0, v * 4.0
        else:
            lo, hi = v - 1, v + 1
        lo, hi = min(lo, v), max(hi, v)
        params.append({"name": name, "value": v, "lo": lo, "hi": hi})

    add("Lipophilicity", first_val(comp.get("Lipophilicity")))
    add("Fraction unbound (plasma, reference value)", first_val(comp.get("FractionUnbound")))
    add("Solubility at reference pH", first_val(comp.get("Solubility")))
    for p in comp.get("Processes") or []:
        for pp in (p.get("Parameters") or []):
            nm = pp.get("Name") or ""
            if nm in ("Intrinsic clearance", "GFR fraction", "kcat", "Km", "Vmax") \
                    or "clearance" in nm.lower():
                add(nm, pp.get("Value"))
    return {"compound": project,
            "methods": {"partition": {"current": cur("partition"), "options": list(PARTITION_METHODS)},
                        "permeability": {"current": cur("permeability"), "options": list(PERMEABILITY_METHODS)}},
            "params": params}



def _topology(project: str, snapd: dict) -> dict:
    """The model's drug-specific structure (the NON-fixed part), parsed from a
    finished snapshot: distribution/permeability methods, metabolizing enzymes,
    renal clearance, transporters, and whether an oral formulation dissolves.
    The whole-body organ network is fixed and drawn client-side."""
    comp = (snapd.get("Compounds") or [{}])[0]
    cms = comp.get("CalculationMethods") or []

    def pick(prefix):
        for m in cms:
            if isinstance(m, str) and m.lower().startswith(prefix):
                return m.split(" - ")[-1].strip()
        return None

    def proc_params(p):
        return {pp.get("Name"): pp.get("Value") for pp in (p.get("Parameters") or [])
                if pp.get("Name")}

    metab, transport = [], []
    renal, renal_params = False, {}
    for p in comp.get("Processes") or []:
        internal = (p.get("InternalName") or "").lower()
        mol = p.get("Molecule")
        if "metaboli" in internal:
            metab.append({"enzyme": mol, "source": p.get("DataSource"),
                          "params": proc_params(p)})
        elif "glomerular" in internal or "gfr" in internal or "renalclear" in internal:
            renal = True
            renal_params = proc_params(p)
        elif "transport" in internal or "active" in internal or "efflux" in internal \
                or "uptake" in internal:
            transport.append({"molecule": mol, "params": proc_params(p)})

    def first_val(block):
        for entry in (block or []):
            for pp in (entry.get("Parameters") or []):
                return pp.get("Value")
        return None
    mw = next((pp.get("Value") for pp in (comp.get("Parameters") or [])
               if pp.get("Name") == "Molecular weight"), None)
    physchem = {"lipophilicity": first_val(comp.get("Lipophilicity")),
                "fraction_unbound": first_val(comp.get("FractionUnbound")),
                "solubility": first_val(comp.get("Solubility")),
                "molecular_weight": mw}
    forms = snapd.get("Formulations") or []

    def ftype(f):
        return str(f.get("FormulationType") or f.get("Type") or f.get("Name") or "").lower()
    oral = any(any(k in ftype(f) for k in ("tablet", "capsule", "weibull", "particle"))
               for f in forms)
    return {"compound": project,
            "distribution": pick("cellular partition"),
            "permeability": pick("cellular permeability"),
            "metabolism": metab, "renal": renal, "renal_params": renal_params,
            "transport": transport, "physchem": physchem,
            "oral": oral, "formulations": [f.get("Name") for f in forms]}


def _reffit(project: str) -> dict:
    """Run a FINISHED reference model (forward simulation only — no fitting) and
    return per-dataset observed + simulated curves for an interactive fit overlay.
    Cached in-memory and on disk, so re-opening is instant. Needs PKSim.CLI."""
    project = _task_dir(project)                        # task id -> base dir
    if project in _REFFIT_CACHE:
        return {"ok": True, "series": _REFFIT_CACHE[project], "cached": True}
    snap = os.path.join(_LIB, project, "json", f"{project}-Model.json")
    inp = os.path.join(_LIB, project, "json_input", f"{project}-Model.input.json")
    if not (os.path.isfile(snap) and os.path.isfile(inp)):
        return {"ok": False, "error": f"no finished snapshot for {project}"}
    disk = os.path.join(_LIB, project, "benchmark", ".reffit.json")
    if os.path.isfile(disk):
        try:
            s = json.load(open(disk, encoding="utf-8"))
            _REFFIT_CACHE[project] = s
            return {"ok": True, "series": s, "cached": True}
        except (OSError, ValueError):
            pass
    from pkpd_agent.config import AgentConfig
    from pkpd_agent.engines.osp_cli import OSPCli
    from pkpd_agent.engines import osp_score
    cfg = AgentConfig(mock=False)
    cli = OSPCli(pksim_cli_path=cfg.pksim_cli_path, timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        return {"ok": False, "error": "PKSim.CLI not found — set PKPD_PKSIM_CLI to run the fit."}
    with _SIM_LOCK:
        if project in _REFFIT_CACHE:                # another request just built it
            return {"ok": True, "series": _REFFIT_CACHE[project], "cached": True}
        observed = json.load(open(inp, encoding="utf-8"))["given_data"]["clinical_observed_data"]
        snapd = json.load(open(snap, encoding="utf-8"))
        res = cli.build_and_run(snap, edits={})
        if not res["ok"]:
            return {"ok": False, "error": res["message"]}
        linkage = osp_score.linkage_from_snapshot(snapd)
        predicted, _ = osp_score.map_predictions(res["profiles"], observed, linkage)
        pred_by_ds = {p["dataset"]: p for p in predicted}
        series = []
        for o in observed:
            ds = o["dataset"]
            obs = [[t, c] for t, c in zip(o.get("time_h", []), o.get("conc_mg_L", []))
                   if c is not None]
            p = pred_by_ds.get(ds)
            sim = _downsample([[t, c] for t, c in zip(p["time_h"], p["pred_conc_mg_L"])
                               if c is not None]) if p else []
            series.append({"study": o.get("study") or ds, "route": o.get("route", ""),
                           "dose": o.get("dose", ""), "observed": obs, "simulated": sim})
        _REFFIT_CACHE[project] = series
        try:
            os.makedirs(os.path.dirname(disk), exist_ok=True)
            json.dump(series, open(disk, "w", encoding="utf-8"))
        except OSError:
            pass
        return {"ok": True, "series": series, "cached": False}


def _report_status(rep: str):
    status, gmfe = "new", None
    if os.path.isfile(rep):
        try:
            r = json.load(open(rep, encoding="utf-8"))
            g = r.get("held_out") or {}
            verdict = (r.get("verdict") or "").lower()
            status = "pass" if "pass" in verdict else ("fail" if "fail" in verdict else "done")
            gmfe = g.get("agent") or r.get("gmfe")
        except Exception:                               # noqa: BLE001
            pass
    return status, gmfe


def _compounds() -> list[dict]:
    """Every task with a blanked snapshot: adult PBPK, pediatric, and DDI variants.
    Each carries a unique `compound` id (dir + variant), the base `dir`, and `kind`."""
    out, seen = [], set()
    plain = sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*.blanked.json")))
    ddis = sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*.ddi_blanked.json")))
    for snap in plain + ddis:
        base = os.path.basename(snap)
        if ".hard_blanked" in base:
            continue
        comp = os.path.basename(os.path.dirname(os.path.dirname(snap)))
        if base.endswith(".ddi_blanked.json"):
            kind, cid = "ddi", f"{comp} · DDI"
            stem = base.replace(".blanked.json", "")
            inp = os.path.join(_LIB, comp, "json_input",
                               base.replace(".ddi_blanked.json", ".ddi_input.json"))
        elif "Pediatric" in base:
            kind, cid = "pediatric", f"{comp} · pediatric"
            stem = base.replace(".blanked.json", "")
            inp = os.path.join(_LIB, comp, "json_input", stem + ".input.json")
        elif base.endswith("-Model.blanked.json"):
            kind, cid = "adult", comp
            stem = base.replace(".blanked.json", "")
            inp = os.path.join(_LIB, comp, "json_input", stem + ".input.json")
        else:
            continue
        if cid in seen or not os.path.isfile(inp):
            continue
        seen.add(cid)
        status, gmfe = _report_status(os.path.join(_LIB, comp, "report", stem + ".json"))
        rec = _load_run(cid)                         # a persisted copilot run, if any
        saved = None
        if rec:
            saved = {"held_out": rec.get("held_out_gmfe"),
                     "in_sample": rec.get("in_sample_gmfe"), "ts": rec.get("ts")}
            if status == "new":                      # no CLI report — reflect the saved run
                status = "done"
                gmfe = rec.get("held_out_gmfe") or rec.get("in_sample_gmfe")
        out.append({"compound": cid, "dir": comp, "label": cid, "kind": kind,
                    "snapshot": snap, "input": inp, "status": status, "gmfe": gmfe,
                    "saved": saved})
    return out


def _find(compound: str) -> dict | None:
    cs = _compounds()
    return (next((c for c in cs if c["compound"] == compound), None)
            or next((c for c in cs if c["dir"] == compound), None))


def _task_dir(x: str) -> str:
    """Resolve a task id (or a bare directory) to its base compound directory."""
    d = _find(x)
    return d["dir"] if d else x


# ---- rerunnable deliverable ----------------------------------------------- #
def _deliverable_readme(compound, best_gmfe, building, held, givens, edits) -> str:
    n_given = sum(1 for g in (givens or []) if (g.get("provenance") or "given") == "given")
    n_fit = len((edits.get("parameters") or {})) if edits else 0
    methods = (edits or {}).get("calculation_methods") or {}
    return (
        f"# {compound} — PBPK model deliverable\n\n"
        "A rerunnable handoff produced by the PBPK Copilot agent. Every parameter "
        "is provenance-tagged so a reviewer can see what was **given** (read from "
        "the literature/report), what was **judged** (fixed by assumption), and "
        "what was **fitted** (identified from the building data).\n\n"
        "## Files\n"
        "- `adopted_model.json` — the structure the agent adopted (distribution / "
        "permeability methods, processes) and the edit spec that reproduces it.\n"
        "- `parameter_table.csv` — parameter, value, unit, source, provenance.\n"
        "- `observed_curves.csv` — the digitized clinical curves (time, conc).\n\n"
        "## Result\n"
        f"- Best in-sample GMFE: **{best_gmfe if best_gmfe is not None else '—'}**\n"
        f"- Distribution method: {methods.get('partition', '—')}; "
        f"permeability: {methods.get('permeability', '—')}\n"
        f"- Givens recorded: {len(givens or [])} ({n_given} read from source); "
        f"parameters fitted: {n_fit}\n\n"
        "## Held-out validation (leave-one-out)\n"
        f"- Building studies (fitted): {', '.join(building) if building else '—'}\n"
        f"- Held-out studies (graded, never seen in the fit): "
        f"{', '.join(held) if held else '—'}\n\n"
        "The agent never saw the reference model's chosen methods, fitted values, "
        "or held-out data. To produce the graded held-out report, run "
        "`examples/run_llm_build.py --report`.\n")


def _write_deliverable(compound, snapshot_path, edits, best_gmfe, givens, split, observed):
    """Assemble a rerunnable handoff folder for one finished run and return its
    path: the adopted model, a provenance-tagged parameter table, the observed
    curves, the held-out split, and a README. Packaging never raises into a run."""
    out = os.path.join(_LIB, compound, "deliverable")
    os.makedirs(out, exist_ok=True)
    edits = edits or {}
    with open(os.path.join(out, "adopted_model.json"), "w", encoding="utf-8") as fh:
        json.dump({"compound": compound, "best_gmfe": best_gmfe,
                   "calculation_methods": edits.get("calculation_methods"),
                   "structure": edits.get("structure"),
                   "estimated_parameters": edits.get("estimate"),
                   "adopted_edits": edits}, fh, indent=2)
    with open(os.path.join(out, "parameter_table.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["parameter", "value", "unit", "source", "provenance"])
        for g in givens or []:
            w.writerow([g.get("parameter"), g.get("value"), g.get("unit", ""),
                        g.get("source", ""), g.get("provenance", "given")])
        for k, v in (edits.get("parameters") or {}).items():
            w.writerow([k, v, "", "optimizer fit to building data", "fitted"])
    with open(os.path.join(out, "observed_curves.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["dataset", "route", "dose", "time_h", "conc_mg_L"])
        for d in observed or []:
            for t, c in zip(d.get("time_h", []), d.get("conc_mg_L", [])):
                w.writerow([d.get("dataset"), d.get("route"), d.get("dose"), t, c])
    with open(os.path.join(out, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(_deliverable_readme(compound, best_gmfe, split.get("building") or [],
                                     split.get("held_out_studies") or [], givens, edits))
    return out


# ---- the agent worker ----------------------------------------------------- #
def _push(q: queue.Queue, ev: dict) -> None:
    q.put(ev)


def _serialize_decision(ev) -> dict:
    return {"type": "decision", "text": getattr(ev, "text", "") or "",
            "calls": [{"name": c.name,
                       "args": {k: v for k, v in (c.arguments or {}).items()}}
                      for c in getattr(ev, "calls", []) or []]}


def _obs_detail(tool: str, c: dict):
    """A compact, size-capped view of what a tool actually returned, so the step
    card can be expanded to show what the agent read/saw (not just a one-liner)."""
    def trim(s, n=2000):
        s = str(s)
        return s if len(s) <= n else s[:n] + " …"

    def molname(m):
        return m.get("molecule") if isinstance(m, dict) else m

    if tool == "osp_read_report":
        md = c.get("report_markdown") or ""
        return {"report excerpt (agent extracts givens from this)": trim(md, 2400)} if md else None
    if tool == "osp_read_reference":
        d = {}
        if c.get("calculation_methods"):
            d["methods"] = c["calculation_methods"]
        procs = c.get("processes") or []
        if procs:
            d["processes"] = [f"{molname(p)}" + (f" ({p.get('type')})" if isinstance(p, dict) else "")
                              for p in procs]
        if c.get("estimated_parameters"):
            d["fitted parameters"] = [str(e if not isinstance(e, dict) else e.get("parameter"))
                                      for e in c["estimated_parameters"]]
        rep = c.get("report_markdown")
        if rep:
            d["analogue report (full write-up)"] = trim(rep, 2600)
        return d or None
    if tool == "osp_inspect":
        d = {}
        if c.get("objective"):
            d["objective"] = trim(c["objective"], 500)
        if c.get("parameters_to_determine"):
            d["to determine (fit these)"] = c["parameters_to_determine"]
        if c.get("given_input_parameters"):
            d["given inputs (trusted)"] = c["given_input_parameters"]
        lit = c.get("literature_physicochemical")
        if lit:
            d["literature values"] = [
                f"{x.get('parameter')} = {x.get('value')} {x.get('unit', '')}".strip()
                for x in lit][:40]
        cm = c.get("candidate_clearance_molecules")
        if cm:
            d["candidate enzymes/transporters"] = cm[:40]
        rl = c.get("reference_library")
        if isinstance(rl, dict):
            avail = rl.get("available") or rl.get("models") or []
            names = [a.get("compound") for a in avail if isinstance(a, dict)]
            d["reference library (open with osp_read_reference)"] = (
                f"{rl.get('mode')} · {len(avail)} analogues: " + ", ".join(n for n in names if n))
        return d or None
    if tool == "osp_options":
        d = {}
        ep = c.get("editable_parameters") or []
        if ep:
            d["editable parameters"] = [f"{p.get('name')} ({p.get('tier', '?')})" for p in ep][:60]
        methods = c.get("calculation_methods") or {}
        if (methods.get("partition") or {}).get("options"):
            d["partition methods"] = methods["partition"]["options"]
        if (methods.get("permeability") or {}).get("options"):
            d["permeability methods"] = methods["permeability"]["options"]
        if c.get("expressed_molecules"):
            d["expressed molecules"] = [molname(m) for m in c["expressed_molecules"]][:40]
        apt = c.get("addable_process_types")
        if apt:
            d["addable process types"] = [
                (p.get("type") if isinstance(p, dict) else str(p)) for p in apt][:40]
        return d or None
    if tool == "osp_list_studies":
        st = c.get("studies") or []
        return {"studies handed": [
            f"{s.get('study')} · {s.get('route', '')} {s.get('dose', '')} · {s.get('n_points')} pts"
            for s in st][:60]} if st else None
    if tool == "osp_read_study":
        pts = c.get("points") or []
        return {"study": f"{c.get('study')} · {c.get('route', '')} {c.get('dose', '')} · "
                f"{len(pts)} points"} if pts else None
    if tool in ("osp_optimize", "osp_sweep_methods"):
        d = {}
        rt = c.get("ranked_top")
        if isinstance(rt, list) and rt:
            d["methods ranked (by data)"] = [
                f"{x.get('partition')} / {x.get('permeability')} → GMFE {x.get('gmfe')}"
                for x in rt]
        opt = c.get("optimized") or (c.get("best") or {}).get("optimized")
        if isinstance(opt, dict) and opt:
            d["fitted parameters"] = [f"{k} = {v}" for k, v in list(opt.items())[:30]]
        if c.get("advice"):
            d["verdict"] = trim(c["advice"], 400)
        if c.get("recommendations"):
            d["recommendations"] = [trim(r, 200) for r in c["recommendations"]][:12]
        if c.get("params_at_bound"):
            d["params at bound"] = c["params_at_bound"]
        br = c.get("by_route")
        if isinstance(br, dict):
            d["by route"] = [f"{k}: GMFE {v.get('gmfe') if isinstance(v, dict) else v}"
                             for k, v in br.items()]
        return d or None
    return None


def _serialize_observation(ev) -> dict:
    c = getattr(ev, "content", {}) or {}
    tool = getattr(ev, "tool", "")
    return {"type": "observation", "tool": tool,
            "message": c.get("message", ""),
            "gmfe": (c.get("best") or {}).get("gmfe") if isinstance(c.get("best"), dict) else c.get("gmfe"),
            "optimized": c.get("optimized"),
            "by_route": c.get("by_route"),
            "ranked_top": c.get("ranked_top"),
            "recorded": c.get("recorded"),        # agent's self-extracted givens table
            "saw_fit_curve": c.get("saw_fit_curve"),   # model was shown a fit overlay
            "defaulted_params": c.get("defaulted_process_params"),  # params silently at default
            "detail": _obs_detail(tool, c)}       # what the agent actually read/saw


class _QueueWriter(io.TextIOBase):
    """Forwards optimizer stdout lines to the SSE queue as {type:'log'}."""
    def __init__(self, q: queue.Queue):
        self.q = q
        self._buf = ""

    def write(self, s):
        self._buf += s
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            line = line.strip()
            if line:
                self.q.put({"type": "log", "line": line})
        return len(s)


def _run_agent(run_id: str, p: dict) -> None:
    q: queue.Queue = _RUNS[run_id]["queue"]
    try:
        with _RUN_SEM:                                  # up to PKPD_COPILOT_JOBS at once
            _run_agent_locked(run_id, p, q)
    except Exception as e:                              # noqa: BLE001
        q.put({"type": "error", "message": f"{type(e).__name__}: {e}"})
    finally:
        q.put({"type": "end"})
        _RUNS[run_id]["done"] = True
        # persist the full trace + outcome so the page can reload/resume it later
        try:
            _persist_run(p.get("compound"), getattr(q, "log", []), p)
        except Exception:                          # noqa: BLE001
            pass


def _run_agent_locked(run_id, p, q) -> None:
    compound = p["compound"]
    model = p["model"]
    max_steps = p["max_steps"]
    library = p.get("library")
    only = p.get("only")
    self_extract = p.get("self_extract")
    note = (p.get("note") or "").strip()
    from pkpd_agent.config import AgentConfig
    from pkpd_agent.engines.osp_cli import OSPCli
    from pkpd_agent.engines import osp_split
    from pkpd_agent.llm import LLMPolicy
    from pkpd_agent.loop import DecisionLoop
    from pkpd_agent.state import Decision, Observation
    from pkpd_agent.tools.registry import ToolRegistry
    from pkpd_agent.tools.osp_loop_tools import register_osp_loop_tools
    from pkpd_agent.tools.context_tools import register_context_tools
    import examples.run_llm_build as R                 # reuse the system prompt

    f = _find(compound)
    if not f:
        q.put({"type": "error", "message": f"no benchmark snapshot for {compound}"})
        return

    # STRUCTURE-BLIND (hard) mode: swap in the hard_blanked snapshot (processes /
    # enzyme identity stripped) + the hard input (a candidate_clearance_molecules
    # pool instead of the given enzyme), so the agent must DISCOVER the mechanism,
    # not just fit values into a given skeleton.
    if p.get("hard"):
        hs = f["snapshot"].replace(".blanked.json", ".hard_blanked.json")
        hi = f["input"].replace(".input.json", ".hard.input.json")
        if os.path.isfile(hs) and os.path.isfile(hi):
            f = {**f, "snapshot": hs, "input": hi}
            q.put({"type": "meta", "hard": True})
        else:
            q.put({"type": "meta", "hard_unavailable": True})

    cfg = AgentConfig(mock=False, max_steps=max_steps)
    cfg.model = model
    if not cfg.anthropic_key_present():
        q.put({"type": "error", "message": "ANTHROPIC_API_KEY not set on the server."})
        return
    cli = OSPCli(pksim_cli_path=cfg.pksim_cli_path, timeout_s=cfg.pksim_timeout_s)
    if not cli.pksim_cli_path or not os.path.exists(cli.pksim_cli_path):
        q.put({"type": "error", "message": f"PKSim.CLI not found ({cli.pksim_cli_path!r}); "
                                           "set PKPD_PKSIM_CLI."})
        return

    inp = json.load(open(f["input"], encoding="utf-8"))
    observed = (inp.get("given_data", {}) or {}).get("clinical_observed_data")
    if not observed:
        q.put({"type": "error", "message": "This task has no plasma concentration data to fit "
               "(e.g. a DDI task) — it runs via the interaction pipeline, not the web runner yet."})
        return
    split = osp_split.split_studies(json.load(open(f["snapshot"], encoding="utf-8")), observed)
    build_set = set(split.get("building") or [])
    if split.get("verification"):
        build_obs = [o for o in observed if o["dataset"] in build_set]
        inp_agent = dict(inp)
        inp_agent["given_data"] = {**inp["given_data"], "clinical_observed_data": build_obs}
    else:
        build_obs, inp_agent = observed, inp

    context_reports: dict = {}                   # compound -> its full context report .md
    if library:
        from pkpd_agent.engines import reference_library as RL
        lib = RL.library_for_snapshot(f["snapshot"],
                                      same_type=(library == "same-type"), full=False,
                                      only=only or None)
        inp_agent = dict(inp_agent)
        inp_agent["reference_library"] = lib
        for mdl in lib.get("models") or []:      # let the agent open each analogue's own report
            comp = mdl.get("compound")
            reps = glob.glob(os.path.join(_RW, comp, "context", "report", "*.md"))
            if reps:
                context_reports[comp] = reps[0]
        q.put({"type": "meta", "library": lib["mode"], "n_ref": len(lib["models"])})

    q.put({"type": "meta", "building": len(build_obs),
           "held_out": len(split.get("verification") or []),
           "studies": len(split.get("held_out_studies") or [])})

    # SELF-EXTRACT: hand the agent the raw report+data so it builds its own
    # context lake, and withhold the pre-digested givens. Confined to the handed
    # materials (no web) so leave-one-out holds.
    ctx_report = ctx_data = None
    if self_extract:
        tdir = os.path.join(_RW, f["dir"], "test")
        reps = glob.glob(os.path.join(tdir, "report", "*.md"))
        ddir = os.path.join(tdir, "data")
        if reps and os.path.isdir(ddir):
            ctx_report, ctx_data = reps[0], ddir
            inp_agent = dict(inp_agent)
            inp_agent["given_data"] = {k: v for k, v in inp_agent["given_data"].items()
                                       if k != "literature_physicochemical"}
            q.put({"type": "meta", "self_extract": True})

    registry = ToolRegistry()
    register_osp_loop_tools(registry, cfg, {
        "cli": cli, "snapshot_path": f["snapshot"],
        "observed": build_obs, "input": inp_agent,
        "context_reports": context_reports,
        "self_extract": bool(ctx_report),
        "fit_vision": p.get("fit_vision", True),     # let the model SEE its fit curves
        "minimal": p.get("minimal", False)})         # purist mode: drop modeling-scaffold
    if ctx_report:
        register_context_tools(registry, cfg, {
            "report_path": ctx_report, "data_dir": ctx_data, "input": inp_agent})
    web = bool(p.get("web"))                 # Anthropic's native, server-side web tools
    if web:                                  # real open-web lookup (NOT leave-one-out)
        q.put({"type": "meta", "web": True})
    if p.get("fit_vision", True):            # the model sees its own fit overlays
        q.put({"type": "meta", "fit_vision": True})
    minimal = bool(p.get("minimal"))
    if minimal:                              # purist mode: tools + invariants + stop rule only
        q.put({"type": "meta", "minimal": True})

    goal = f"{inp.get('objective', 'Build the PBPK model.')}\n\n"
    if note:                                 # modeler steering from the composer
        goal += f"Modeler note: {note}\n\n"
    if web:
        goal += ("You have web search and web fetch: look up this compound's DMPK "
                 "(clearing enzyme, physchem, transporters) and any already-published "
                 "PBPK model, and build on what has been done - cite what you use.\n\n")
    goal += "Start with osp_inspect, then determine the model and call osp_optimize."
    policy = LLMPolicy(cfg, registry,
                       R._system_prompt(1.6, self_extract=bool(ctx_report), web=web,
                                        minimal=minimal),
                       web=web,
                       on_web=lambda e: q.put({"type": "web", **e}))   # stream web activity live
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    def on_event(ev):
        if isinstance(ev, Decision):
            q.put(_serialize_decision(ev))
        elif isinstance(ev, Observation):
            q.put(_serialize_observation(ev))

    q.put({"type": "status", "phase": "running"})
    writer = _QueueWriter(q)
    router = _install_router()                          # route THIS thread's prints to q
    router.register(writer)
    # cooperative cancel: /api/cancel sets this flag; the loop checks it at each
    # step/tool boundary and stops cleanly (frees the worker slot).
    should_stop = lambda: bool(_RUNS.get(run_id, {}).get("cancel"))
    try:
        session = loop.run(goal, on_event=on_event, should_stop=should_stop)
    finally:
        router.flush()
        router.unregister()

    best = session.get("osp_best_gmfe")
    edits = session.get("osp_best_edits")
    givens = session.get("extracted_givens") or \
        (inp_agent.get("given_data", {}) or {}).get("literature_physicochemical")
    deliverable = None
    try:
        deliverable = _write_deliverable(f["dir"], f["snapshot"], edits, best,
                                         givens, split, build_obs)
    except Exception:                        # noqa: BLE001 - packaging must never sink a run
        pass
    web_lookups = session.get("web_lookups") or []
    q.put({"type": "done", "best_gmfe": best, "best_edits": edits,
           "givens": givens, "deliverable": bool(deliverable),
           # transparency: what the run consulted on the open web (non-blind if any)
           "web_lookups": web_lookups, "blind": (not web_lookups)})

    # The REPORT is part of the pipeline, not a separate click/HTTP call. Compute it
    # HERE, on the same worker, right after the build: its PK-Sim runs serialize with
    # the build via _SIM_LOCK (no pile-up), it streams back over the SAME SSE as a
    # 'report' event before 'end' (so it can't hit an HTTP timeout), and because the
    # batch worker waits for 'end', the next compound only starts once this report is
    # done. Best-effort: a report failure never sinks the run.
    if edits:
        try:
            rep = _report(f, edits, model)
            q.put({"type": "report", "compound": compound, **rep})
            # POST-STOP reveal: the agent sees the reference and critiques its own
            # build. Strictly post-hoc (the graded run is already done) - never leaks.
            pm = _postmortem(compound, rep, edits, model)
            if pm:
                q.put({"type": "postmortem", "compound": compound, "text": pm})
        except Exception as e:                       # noqa: BLE001
            q.put({"type": "report", "compound": compound,
                   "error": f"{type(e).__name__}: {e}"})


# ---- FastAPI app ---------------------------------------------------------- #
def create_app():
    import hmac
    from fastapi import FastAPI
    from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse

    app = FastAPI(title="PBPK Copilot")

    @app.on_event("startup")
    async def _bump_threadpool():
        # each open SSE stream and each sync PK-Sim endpoint holds a threadpool
        # worker; the default (~40) can starve so a new /api/model call hangs on
        # "loading…". Raise the ceiling so the UI stays responsive during runs.
        try:
            import anyio
            anyio.to_thread.current_default_thread_limiter().total_tokens = 256
        except Exception:                               # noqa: BLE001
            pass

    # Optional access token: set PBPK_COPILOT_TOKEN to require it (needed before
    # exposing the server beyond localhost / a trusted LAN). Pass it once as
    # ?token=... — the server sets a cookie, so links, fetch and the SSE stream
    # all carry it afterwards. Unset -> no auth (unchanged local behavior).
    TOKEN = os.environ.get("PBPK_COPILOT_TOKEN")

    @app.middleware("http")
    async def _auth(request, call_next):
        if not TOKEN:
            return await call_next(request)
        supplied = request.cookies.get("pbpk_token") or request.query_params.get("token")
        if supplied and hmac.compare_digest(supplied, TOKEN):
            resp = await call_next(request)
            if request.query_params.get("token"):        # first visit -> remember it
                resp.set_cookie("pbpk_token", TOKEN, httponly=True, samesite="lax")
            return resp
        return HTMLResponse(
            "<body style='font:15px system-ui;max-width:34em;margin:12vh auto;padding:0 6vw'>"
            "<h2>PBPK Copilot</h2><p>This server requires an access token. Open it as "
            "<code>?token=YOUR_TOKEN</code> (the token set in <code>PBPK_COPILOT_TOKEN</code> "
            "on the host).</p></body>", status_code=401)

    @app.get("/", response_class=HTMLResponse)
    def index():
        with open(os.path.join(_HERE, "static", "index.html"), encoding="utf-8") as fh:
            return fh.read()

    @app.get("/api/models")
    def models():
        return JSONResponse([{k: c.get(k) for k in
                              ("compound", "status", "gmfe", "label", "kind", "dir", "saved")}
                             for c in _compounds()])

    @app.get("/api/library")
    def library(compound: str, same_type: bool = True):
        from pkpd_agent.engines import reference_library as RL
        f = _find(compound)
        if not f:
            return JSONResponse({"models": []})
        return JSONResponse(RL.library_for_snapshot(f["snapshot"],
                                                    same_type=same_type, full=False))

    @app.get("/api/rw/projects")
    def rw_projects():
        return JSONResponse(_rw_projects())

    @app.get("/api/rw/report")
    def rw_report(project: str, kind: str = "test"):
        d = os.path.join(_RW, _task_dir(project), kind, "report")
        files = glob.glob(os.path.join(d, "*.md"))
        if not files:
            return JSONResponse({"markdown": ""})
        return JSONResponse({"markdown": open(files[0], encoding="utf-8").read()})

    @app.get("/api/rw/series")
    def rw_series(project: str, kind: str = "test"):
        return JSONResponse({"series": _rw_series(_task_dir(project), kind)})

    @app.get("/api/catalog")
    def catalog():
        return JSONResponse(_catalog())

    @app.get("/api/model")
    def model_view(compound: str):
        f = _find(compound)
        if not f:
            return JSONResponse({"error": "unknown compound"}, status_code=404)
        return JSONResponse(_model_view(f))

    @app.post("/api/try")
    def try_model(payload: dict):
        compound = payload.get("compound")
        edits = payload.get("edits") or {}
        f = _find(compound)
        if not f:
            return JSONResponse({"error": "unknown compound"}, status_code=404)
        return JSONResponse(_try_model(f, edits))

    @app.post("/api/grade")
    def grade(payload: dict):
        # blind held-out grade of an adopted model (forward run on the verification
        # studies the agent never fit). edits should carry fitted `parameters`.
        compound = payload.get("compound")
        edits = payload.get("edits") or {}
        f = _find(compound)
        if not f:
            return JSONResponse({"error": "unknown compound"}, status_code=404)
        return JSONResponse(_try_model(f, edits, subset="held_out"))

    @app.post("/api/report")
    def report(payload: dict):
        # post-build report: in-sample + held-out curves/grades + comparison vs the
        # reference answer (revealed to the human only, never fed to the agent).
        f = _find(payload.get("compound"))
        if not f:
            return JSONResponse({"error": "unknown compound"}, status_code=404)
        return JSONResponse(_report(f, payload.get("edits") or {}, payload.get("model")))

    @app.get("/api/runs/{compound}")
    def run_record(compound: str):
        # the persisted trace + outcome of a finished run, so the page can replay it
        # after a reload without re-running the agent.
        rec = _load_run(compound)
        if not rec:
            return JSONResponse({"error": "no saved run"}, status_code=404)
        return JSONResponse(rec)

    @app.delete("/api/runs/{compound}")
    def run_delete(compound: str):
        # discard a saved run so the compound is treated as un-built (a fresh
        # re-run starts from scratch and the batch no longer skips it).
        path = _run_record_path(compound)
        try:
            if os.path.isfile(path):
                os.remove(path)
        except OSError:
            return JSONResponse({"ok": False}, status_code=500)
        return JSONResponse({"ok": True, "deleted": compound})

    @app.post("/api/runs/{compound}/grade")
    def run_grade(compound: str, payload: dict):
        # persist the held-out grade the scoreboard computed, so a reload shows it.
        rec = _load_run(compound)
        if not rec:
            return JSONResponse({"ok": False, "error": "no saved run"}, status_code=404)
        rec["held_out_gmfe"] = payload.get("held_out")
        try:
            with open(_run_record_path(compound), "w", encoding="utf-8") as fh:
                json.dump(rec, fh)
        except OSError:
            return JSONResponse({"ok": False}, status_code=500)
        return JSONResponse({"ok": True})

    @app.get("/api/built")
    def built(compound: str):
        return JSONResponse(_built_view(compound))

    @app.get("/api/rw/topology")
    def rw_topology(project: str):
        project = _task_dir(project)
        snap = os.path.join(_LIB, project, "json", f"{project}-Model.json")
        if not os.path.isfile(snap):
            return JSONResponse({"error": "no snapshot"}, status_code=404)
        try:
            snapd = json.load(open(snap, encoding="utf-8"))
        except (OSError, ValueError):
            return JSONResponse({"error": "unreadable snapshot"}, status_code=404)
        return JSONResponse(_topology(project, snapd))

    @app.get("/api/rw/simulate")
    def rw_simulate(project: str):
        # run the finished reference model to get its simulated fit (cached)
        return JSONResponse(_reffit(project))

    @app.get("/api/rw/figure")
    def rw_figure(project: str, path: str):
        # the context report's simulated-vs-observed fit figures live in the OSP
        # library (images/...); serve them so the reference fit renders in-viewer.
        from fastapi.responses import FileResponse, Response
        base = os.path.normpath(os.path.join(_LIB, _task_dir(project)))
        fp = os.path.normpath(os.path.join(base, path))
        if (not fp.startswith(base + os.sep) or not fp.lower().endswith(".png")
                or not os.path.isfile(fp)):
            return Response(status_code=404)
        return FileResponse(fp, media_type="image/png")

    @app.post("/api/run")
    def run(payload: dict):
        compound = payload.get("compound")
        model = payload.get("model") or "claude-sonnet-5"
        max_steps = int(payload.get("max_steps") or 10)
        # default: the agent sees ALL other projects' context; the selector narrows it.
        only = payload.get("context_projects")          # None (= all) | [compounds]
        lib = payload.get("library") or ("all" if only is None or only else None)
        # default: agent builds its own context lake when a realworld tree exists.
        se = payload.get("self_extract")
        self_extract = os.path.isdir(_RW) if se is None else bool(se)
        run_id = uuid.uuid4().hex[:12]
        _RUNS[run_id] = {"queue": _LoggingQueue(), "done": False, "compound": compound}
        params = {"compound": compound, "model": model, "max_steps": max_steps,
                  "library": lib, "only": only, "self_extract": self_extract,
                  "hard": bool(payload.get("hard")), "web": bool(payload.get("web")),
                  # fit-vision defaults ON: the modeler chose "best modeling"
                  "fit_vision": bool(payload.get("fit_vision", True)),
                  "minimal": bool(payload.get("minimal")),   # purist mode (default off)
                  "note": payload.get("note")}
        threading.Thread(target=_run_agent, args=(run_id, params), daemon=True).start()
        return JSONResponse({"run_id": run_id})

    @app.post("/api/cancel/{run_id}")
    def cancel(run_id: str):
        # cooperatively stop the server-side agent run: the loop checks this
        # flag at each step/tool boundary and finishes cleanly, freeing the slot.
        run = _RUNS.get(run_id)
        if not run:
            return JSONResponse({"ok": False, "error": "unknown run"}, status_code=404)
        run["cancel"] = True
        return JSONResponse({"ok": True, "cancelled": run_id, "done": run.get("done", False)})

    @app.get("/api/deliverable")
    def deliverable(compound: str):
        from fastapi.responses import Response
        d = os.path.join(_LIB, _task_dir(compound), "deliverable")
        if not os.path.isdir(d):
            return JSONResponse({"error": "no deliverable yet — run the agent first"},
                                status_code=404)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for name in sorted(os.listdir(d)):
                fp = os.path.join(d, name)
                if os.path.isfile(fp):
                    z.write(fp, arcname=f"{compound}_deliverable/{name}")
        return Response(buf.getvalue(), media_type="application/zip",
                        headers={"Content-Disposition":
                                 f'attachment; filename="{compound}_deliverable.zip"'})

    @app.get("/api/stream/{run_id}")
    def stream(run_id: str):
        run = _RUNS.get(run_id)
        if not run:
            return JSONResponse({"error": "unknown run"}, status_code=404)
        q = run["queue"]

        def gen():
            # multi-watcher: each connection gets its own subscriber queue, and a
            # late joiner (a second page, or one opened after the run started)
            # first replays the backlog, then streams live - so two pages on the
            # same run stay in sync instead of stealing each other's events.
            if hasattr(q, "subscribe"):
                backlog, sub = q.subscribe()
            else:                                     # pragma: no cover - plain Queue fallback
                backlog, sub = [], q
            try:
                for ev in backlog:
                    yield f"data: {json.dumps(ev)}\n\n"
                    if ev.get("type") == "end":       # run already finished -> replay is complete
                        return
                while True:
                    try:
                        ev = sub.get(timeout=30)
                    except queue.Empty:
                        yield ": keep-alive\n\n"
                        continue
                    yield f"data: {json.dumps(ev)}\n\n"
                    if ev.get("type") == "end":
                        break
            finally:
                if hasattr(q, "unsubscribe"):
                    q.unsubscribe(sub)
        return StreamingResponse(gen(), media_type="text/event-stream")

    @app.get("/api/running")
    def running():
        """Live runs, so a page (even on another machine) can DISCOVER and attach
        to a build already in flight - it does not need to have started it."""
        return JSONResponse([
            {"run_id": rid, "compound": r.get("compound"), "done": bool(r.get("done"))}
            for rid, r in _RUNS.items() if not r.get("done")])

    return app


def _lan_ips() -> list[str]:
    """Best-effort local IPv4 addresses, for the 'open from another computer' case."""
    import socket
    ips = []
    try:                                            # the IP that routes outbound
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ips.append(s.getsockname()[0])
        s.close()
    except OSError:
        pass
    try:
        for info in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET):
            ip = info[4][0]
            if ip not in ips and not ip.startswith("127."):
                ips.append(ip)
    except OSError:
        pass
    return ips


def _clear_cache() -> int:
    """Delete the on-disk reference-fit cache (the *.reffit.json overlay files
    the '▶ Run the model' button writes). The agent build itself is never
    cached; /api/try and Explore run live in throwaway temp dirs. The in-memory
    cache clears whenever the server restarts."""
    import glob as _glob
    n = 0
    for fp in _glob.glob(os.path.join(_LIB, "*", "benchmark", ".reffit.json")):
        try:
            os.remove(fp); n += 1
        except OSError:
            pass
    _REFFIT_CACHE.clear()
    return n


def main():
    import sys
    if "--clear-cache" in sys.argv[1:]:
        n = _clear_cache()
        print(f"cleared {n} cached reference-fit file(s) under {_LIB}")
        print("(the agent build is never cached; the in-memory cache clears on restart)")
        return
    import uvicorn
    host = os.environ.get("PBPK_COPILOT_HOST", "127.0.0.1")
    port = int(os.environ.get("PBPK_COPILOT_PORT", "8765"))
    token = os.environ.get("PBPK_COPILOT_TOKEN")
    suffix = f"/?token={token}" if token else ""
    print(f"PBPK Copilot -> http://{host}:{port}{suffix}")
    if host in ("0.0.0.0", "::"):                   # bound to all interfaces
        for ip in _lan_ips():
            print(f"   from another computer on this network -> http://{ip}:{port}{suffix}")
    else:
        print("   local only. Set PBPK_COPILOT_HOST=0.0.0.0 to reach it from other machines.")
    if token:
        print("   access token is ON — only URLs with ?token=... get in.")
    else:
        print("   NO access token. Set PBPK_COPILOT_TOKEN before exposing it beyond a "
              "trusted LAN (there is no other login; the API key/compute run on THIS host).")
    print("   To reach it from ANY network: run a tunnel to this port, e.g. "
          "`cloudflared tunnel --url http://localhost:%d` (public https URL) or Tailscale "
          "(private). Keep the token on when you do." % port)
    uvicorn.run(create_app(), host=host, port=port, log_level="info")


if __name__ == "__main__":
    main()

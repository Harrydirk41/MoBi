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

import glob
import io
import json
import os
import queue
import threading
import uuid
from contextlib import redirect_stdout

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
_RUN_LOCK = threading.Lock()                        # serialize runs (stdout capture)


def _compounds() -> list[dict]:
    """Compounds that have a normal blanked snapshot, with a report status if any."""
    out = []
    for snap in sorted(glob.glob(os.path.join(_LIB, "*", "benchmark", "*-Model.blanked.json"))
                       + glob.glob(os.path.join(_LIB, "*", "benchmark", "*.blanked.json"))):
        base = os.path.basename(snap)
        if ".hard_blanked" in base or "Pediatric" in base:
            continue
        comp = os.path.basename(os.path.dirname(os.path.dirname(snap)))
        if any(o["compound"] == comp for o in out):
            continue
        stem = base.replace(".blanked.json", "")
        inp = os.path.join(_LIB, comp, "json_input", stem + ".input.json")
        if not os.path.isfile(inp):
            continue
        rep = os.path.join(_LIB, comp, "report", stem + ".json")
        status, gmfe = "new", None
        if os.path.isfile(rep):
            try:
                r = json.load(open(rep, encoding="utf-8"))
                g = r.get("held_out") or {}
                verdict = (r.get("verdict") or "").lower()
                status = "pass" if "pass" in verdict else ("fail" if "fail" in verdict else "done")
                gmfe = g.get("agent") or r.get("gmfe")
            except Exception:
                pass
        out.append({"compound": comp, "snapshot": snap, "input": inp,
                    "status": status, "gmfe": gmfe})
    return out


def _find(compound: str) -> dict | None:
    return next((c for c in _compounds() if c["compound"] == compound), None)


# ---- the agent worker ----------------------------------------------------- #
def _push(q: queue.Queue, ev: dict) -> None:
    q.put(ev)


def _serialize_decision(ev) -> dict:
    return {"type": "decision", "text": getattr(ev, "text", "") or "",
            "calls": [{"name": c.name,
                       "args": {k: v for k, v in (c.arguments or {}).items()}}
                      for c in getattr(ev, "calls", []) or []]}


def _serialize_observation(ev) -> dict:
    c = getattr(ev, "content", {}) or {}
    return {"type": "observation", "tool": getattr(ev, "tool", ""),
            "message": c.get("message", ""),
            "gmfe": (c.get("best") or {}).get("gmfe") if isinstance(c.get("best"), dict) else c.get("gmfe"),
            "optimized": c.get("optimized"),
            "by_route": c.get("by_route"),
            "ranked_top": c.get("ranked_top")}


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


def _run_agent(run_id: str, compound: str, model: str, max_steps: int,
               library: str | None, only: list | None, self_extract: bool) -> None:
    q: queue.Queue = _RUNS[run_id]["queue"]
    try:
        with _RUN_LOCK:
            _run_agent_locked(run_id, compound, model, max_steps, library, only,
                              self_extract, q)
    except Exception as e:                              # noqa: BLE001
        q.put({"type": "error", "message": f"{type(e).__name__}: {e}"})
    finally:
        q.put({"type": "end"})
        _RUNS[run_id]["done"] = True


def _run_agent_locked(run_id, compound, model, max_steps, library, only,
                      self_extract, q) -> None:
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
    observed = inp["given_data"]["clinical_observed_data"]
    split = osp_split.split_studies(json.load(open(f["snapshot"], encoding="utf-8")), observed)
    build_set = set(split.get("building") or [])
    if split.get("verification"):
        build_obs = [o for o in observed if o["dataset"] in build_set]
        inp_agent = dict(inp)
        inp_agent["given_data"] = {**inp["given_data"], "clinical_observed_data": build_obs}
    else:
        build_obs, inp_agent = observed, inp

    if library:
        from pkpd_agent.engines import reference_library as RL
        lib = RL.library_for_snapshot(f["snapshot"],
                                      same_type=(library == "same-type"), full=False,
                                      only=only or None)
        inp_agent = dict(inp_agent)
        inp_agent["reference_library"] = lib
        q.put({"type": "meta", "library": lib["mode"], "n_ref": len(lib["models"])})

    q.put({"type": "meta", "building": len(build_obs),
           "held_out": len(split.get("verification") or []),
           "studies": len(split.get("held_out_studies") or [])})

    # SELF-EXTRACT: hand the agent the raw report+data so it builds its own
    # context lake, and withhold the pre-digested givens. Confined to the handed
    # materials (no web) so leave-one-out holds.
    ctx_report = ctx_data = None
    if self_extract:
        tdir = os.path.join(_RW, compound, "test")
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
        "observed": build_obs, "input": inp_agent})
    if ctx_report:
        register_context_tools(registry, cfg, {
            "report_path": ctx_report, "data_dir": ctx_data, "input": inp_agent})

    goal = (f"{inp.get('objective', 'Build the PBPK model.')}\n\nStart with osp_inspect, "
            "then determine the model and call osp_optimize.")
    policy = LLMPolicy(cfg, registry, R._system_prompt(1.6, self_extract=bool(ctx_report)))
    loop = DecisionLoop(config=cfg, registry=registry, policy=policy)

    def on_event(ev):
        if isinstance(ev, Decision):
            q.put(_serialize_decision(ev))
        elif isinstance(ev, Observation):
            q.put(_serialize_observation(ev))

    q.put({"type": "status", "phase": "running"})
    writer = _QueueWriter(q)
    with redirect_stdout(writer):
        session = loop.run(goal, on_event=on_event)

    best = session.get("osp_best_gmfe")
    edits = session.get("osp_best_edits")
    q.put({"type": "done", "best_gmfe": best, "best_edits": edits})


# ---- FastAPI app ---------------------------------------------------------- #
def create_app():
    from fastapi import FastAPI
    from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse

    app = FastAPI(title="PBPK Copilot")

    @app.get("/", response_class=HTMLResponse)
    def index():
        with open(os.path.join(_HERE, "static", "index.html"), encoding="utf-8") as fh:
            return fh.read()

    @app.get("/api/models")
    def models():
        return JSONResponse([{k: c[k] for k in ("compound", "status", "gmfe")}
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
        d = os.path.join(_RW, project, kind, "report")
        files = glob.glob(os.path.join(d, "*.md"))
        if not files:
            return JSONResponse({"markdown": ""})
        return JSONResponse({"markdown": open(files[0], encoding="utf-8").read()})

    @app.get("/api/rw/series")
    def rw_series(project: str, kind: str = "test"):
        return JSONResponse({"series": _rw_series(project, kind)})

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
        _RUNS[run_id] = {"queue": queue.Queue(), "done": False}
        threading.Thread(target=_run_agent,
                         args=(run_id, compound, model, max_steps, lib, only, self_extract),
                         daemon=True).start()
        return JSONResponse({"run_id": run_id})

    @app.get("/api/stream/{run_id}")
    def stream(run_id: str):
        run = _RUNS.get(run_id)
        if not run:
            return JSONResponse({"error": "unknown run"}, status_code=404)
        q: queue.Queue = run["queue"]

        def gen():
            while True:
                try:
                    ev = q.get(timeout=30)
                except queue.Empty:
                    yield ": keep-alive\n\n"
                    continue
                yield f"data: {json.dumps(ev)}\n\n"
                if ev.get("type") == "end":
                    break
        return StreamingResponse(gen(), media_type="text/event-stream")

    return app


def main():
    import uvicorn
    host = os.environ.get("PBPK_COPILOT_HOST", "127.0.0.1")
    port = int(os.environ.get("PBPK_COPILOT_PORT", "8765"))
    print(f"PBPK Copilot -> http://{host}:{port}")
    uvicorn.run(create_app(), host=host, port=port, log_level="info")


if __name__ == "__main__":
    main()

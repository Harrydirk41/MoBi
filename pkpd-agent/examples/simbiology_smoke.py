r"""Smoke-test the Python -> MATLAB Engine -> SimBiology toolchain.

Proves, step by step and entirely headless (no GUI):
  1. the MATLAB Engine for Python is installed and starts,
  2. SimBiology is licensed/available,
  3. a SimBiology model can be BUILT + SIMULATED from code (sb_smoke.m),
  4. (optional) an existing .sbproj QSP project loads and simulates
     (sb_project_info.m) - e.g. the Vantage RA model.

Every step is flushed to the console AND mirrored to simbiology_smoke.log, and each
engine call is bracketed by ->/<- markers, so a silent stop (buffered output or a
native engine crash) still pinpoints exactly which call died.

    python -m examples.simbiology_smoke
    python -m examples.simbiology_smoke --sbproj "C:\path\to\Vantage RA QSP Model v1.0.sbproj"
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import traceback

_MATLAB_ROOT = r"C:\Program Files\MATLAB\R2026a"
_LOG = None


def log(msg: str) -> None:
    print(msg, flush=True)
    if _LOG is not None:
        try:
            _LOG.write(msg + "\n")
            _LOG.flush()
        except Exception:                              # noqa: BLE001
            pass


def _install_hint():
    log("[FAIL] `matlab.engine` is not installed for this Python.")
    log("       Install the MATLAB Engine for Python (in-tree, from an ADMIN shell):")
    log(f'         cd "{os.path.join(_MATLAB_ROOT, "extern", "engines", "python")}"')
    log("         python -m pip install .")


def _flat(md) -> list[float]:
    try:
        return [float(v) for v in md._data]
    except Exception:                                  # noqa: BLE001
        out = []
        try:
            for row in md:
                if hasattr(row, "__iter__"):
                    out.extend(float(x) for x in row)
                else:
                    out.append(float(row))
        except TypeError:
            out.append(float(md))
        return out


def run(args) -> None:
    try:
        import matlab.engine
    except Exception:                                  # noqa: BLE001
        _install_hint()
        return
    log("[ok] matlab.engine imported")

    log("[..] starting MATLAB engine (first start can take ~15-30s)...")
    t0 = time.time()
    eng = matlab.engine.start_matlab()
    log(f"[ok] engine started in {time.time() - t0:.1f}s")

    try:
        log("-> eng.version")
        try:
            log(f"[ok] MATLAB {eng.version(nargout=1)}")
        except Exception as exc:                       # noqa: BLE001
            log(f"[warn] version: {exc}")

        log("-> eng.license test SimBiology")
        try:
            has = eng.license("test", "SimBiology", nargout=1)
            log(f"[{'ok' if bool(has) else 'FAIL'}] SimBiology license test -> {has}")
            if not bool(has):
                log("       SimBiology not licensed; stopping.")
                return
        except Exception as exc:                       # noqa: BLE001
            log(f"[warn] license test errored (continuing): {exc}")

        here = os.path.dirname(os.path.abspath(__file__))
        mldir = os.path.join(here, "matlab")
        log(f"-> eng.addpath {mldir}")
        eng.addpath(mldir, nargout=0)
        log("<- addpath done")

        # engine sanity: a trivial non-SimBiology call after addpath
        log("-> eng.eval '1+1' (engine sanity)")
        log(f"<- 1+1 = {eng.eval('1+1', nargout=1)}")

        # build the 1-cmt model INCREMENTALLY from Python so a hard engine crash
        # pinpoints the exact SimBiology call (the last '->' with no '<-').
        steps = [
            ("create model",    "m = sbiomodel('smoke');"),
            ("add compartment", "cc = addcompartment(m,'central',1.0);"),
            ("add species",     "ss = addspecies(cc,'drug',10.0);"),
            ("add parameter",   "pp = addparameter(m,'ke',0.5);"),
            ("add reaction",    "rr = addreaction(m,'drug -> null');"),
            ("add kinetic law", "kl = addkineticlaw(rr,'MassAction'); kl.ParameterVariableNames = {'ke'};"),
            ("configset",       "cs = getconfigset(m); cs.StopTime = 24;"),
            ("simulate",        "sd = sbiosimulate(m);"),
            ("select drug",     "dd = selectbyname(sd,'drug');"),
        ]
        crashed = False
        for label, code in steps:
            log(f"-> {label}: {code}")
            try:
                eng.eval(code, nargout=0)
                log(f"<- {label} ok")
            except Exception as exc:                   # noqa: BLE001
                log(f"[FAIL] {label}: {exc}")
                log(traceback.format_exc())
                crashed = True
                break
        if crashed:
            return
        # Returning ARRAYS to Python (bare-expr eval OR eng.workspace) crashes this
        # engine build; scalars/strings marshal fine. So reduce to SCALARS in
        # MATLAB and read those. (The real engine will move arrays via files.)
        log("-> reduce results to scalars in MATLAB")
        try:
            eng.eval("smoke_n = numel(smoke_t); smoke_c0 = smoke_c(1); "
                     "smoke_cend = smoke_c(end); smoke_tend = smoke_t(end);",
                     nargout=0)
            log("<- reduced")
        except Exception as exc:                       # noqa: BLE001
            log(f"[FAIL] reduce: {exc}")
            log(traceback.format_exc())
            return
        log("-> read smoke_n via workspace (scalar)")
        n = eng.workspace["smoke_n"]
        log(f"<- smoke_n = {n}")
        c0 = eng.workspace["smoke_c0"]
        cend = eng.workspace["smoke_cend"]
        tend = eng.workspace["smoke_tend"]
        log(f"[ok] built + simulated a 1-cmt model: {int(n)} timepoints, "
            f"C(0)={float(c0):.3f} -> C(end)={float(cend):.3f} over t=0..{float(tend):g}h")
        log("<- model build/simulate done")

        if args.sbproj:
            if not os.path.exists(args.sbproj):
                log(f"[FAIL] --sbproj not found: {args.sbproj}")
            else:
                log(f"-> eng.sb_project_info {args.sbproj}")
                # let the .m do the reporting via fprintf and CAPTURE MATLAB's
                # stdout in Python (no struct marshaled back - another spot the
                # engine could choke on).
                import io
                out, err = io.StringIO(), io.StringIO()
                try:
                    eng.sb_project_info(args.sbproj, nargout=0,
                                        stdout=out, stderr=err)
                    txt = out.getvalue().strip()
                    log(txt if txt else "(no output captured from MATLAB)")
                    etxt = err.getvalue().strip()
                    if etxt:
                        log("   MATLAB stderr: " + etxt)
                    log("[ok] project loaded + baseline simulated")
                except Exception as exc:               # noqa: BLE001
                    log(f"[FAIL] sb_project_info: {exc}")
                    so = out.getvalue().strip()
                    se = err.getvalue().strip()
                    if so:
                        log("   stdout so far: " + so)
                    if se:
                        log("   stderr: " + se)
                    log(traceback.format_exc())
                log("<- sb_project_info done")

        log("=== SMOKE END (reached the end cleanly) ===")
    finally:
        if not args.keep_open:
            try:
                eng.quit()
            except Exception:                          # noqa: BLE001
                pass


def main() -> None:
    global _LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", default=None)
    ap.add_argument("--keep-open", action="store_true")
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(line_buffering=True)     # no lost buffered output
    except Exception:                                   # noqa: BLE001
        pass
    logpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "simbiology_smoke.log")
    try:
        _LOG = open(logpath, "w", encoding="utf-8")
    except Exception:                                   # noqa: BLE001
        _LOG = None
    log(f"(logging to {logpath})")
    try:
        run(args)
    except Exception as exc:                             # noqa: BLE001
        log(f"[FATAL] {type(exc).__name__}: {exc}")
        log(traceback.format_exc())
    finally:
        if _LOG is not None:
            _LOG.close()


if __name__ == "__main__":
    main()

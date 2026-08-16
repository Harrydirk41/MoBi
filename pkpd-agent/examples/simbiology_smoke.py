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

        log("-> eng.sb_smoke (build + simulate a 1-cmt model)")
        try:
            t, c = eng.sb_smoke(nargout=2)
            tv, cv = _flat(t), _flat(c)
            log(f"[ok] built + simulated a 1-cmt model: {len(tv)} timepoints, "
                f"C(0)={cv[0]:.3f} -> C(end)={cv[-1]:.3f} over t=0..{tv[-1]:g}h")
        except Exception as exc:                       # noqa: BLE001
            log(f"[FAIL] sb_smoke build/simulate failed: {exc}")
            log(traceback.format_exc())
            return
        log("<- sb_smoke done")

        if args.sbproj:
            if not os.path.exists(args.sbproj):
                log(f"[FAIL] --sbproj not found: {args.sbproj}")
            else:
                log(f"-> eng.sb_project_info {args.sbproj}")
                try:
                    info = eng.sb_project_info(args.sbproj, nargout=1)
                    log(f"[ok] loaded project '{info.get('name')}': "
                        f"{int(info.get('nSpecies', 0))} species, "
                        f"{int(info.get('nParameters', 0))} params, "
                        f"{int(info.get('nReactions', 0))} reactions, "
                        f"{int(info.get('nVariants', 0))} variants, "
                        f"{int(info.get('nDoses', 0))} doses; "
                        f"baseline sim ok={info.get('baselineOk')}")
                    vn = info.get("variantNames") or []
                    dn = info.get("doseNames") or []
                    if vn:
                        log(f"     first variants: {list(vn)[:5]}")
                    if dn:
                        log(f"     doses: {list(dn)}")
                except Exception as exc:               # noqa: BLE001
                    log(f"[FAIL] sb_project_info failed: {exc}")
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

r"""Smoke-test the Python -> MATLAB Engine -> SimBiology toolchain.

Proves, step by step and entirely headless (no GUI):
  1. the MATLAB Engine for Python is installed and starts,
  2. SimBiology is licensed/available,
  3. a SimBiology model can be BUILT + SIMULATED from code (sb_smoke.m),
  4. (optional) an existing .sbproj QSP project loads and simulates
     (sb_project_info.m) - e.g. the Vantage RA model.

Each step prints [ok]/[FAIL] so a failure pinpoints exactly where the toolchain
breaks. Run it on the machine that has MATLAB + SimBiology:

    python -m examples.simbiology_smoke
    python -m examples.simbiology_smoke --sbproj "C:\path\to\Vantage RA QSP Model v1.0.sbproj"

If the engine import fails, the script prints the exact install commands for your
MATLAB release.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

_MATLAB_ROOT = r"C:\Program Files\MATLAB\R2026a"


def _install_hint():
    print("[FAIL] `matlab.engine` is not installed for this Python.")
    print("       Install the MATLAB Engine for Python (version-proof, in-tree):")
    print(f'         cd "{os.path.join(_MATLAB_ROOT, "extern", "engines", "python")}"')
    print("         python -m pip install .")
    print("       (or from PyPI, matching the release: pip install matlabengine==26.1.*)")
    print("       Note: the Python running THIS script must be a version your")
    print("       MATLAB release supports, and be the same interpreter you pip-install into.")


def _flat(md) -> list[float]:
    """Flatten a matlab.double (Nx1 / 1xN / NxM) to a Python float list."""
    try:
        return [float(v) for v in md._data]          # fast path: internal flat buffer
    except Exception:
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sbproj", default=None, help="optional .sbproj to load + probe")
    ap.add_argument("--keep-open", action="store_true",
                    help="leave the MATLAB engine running (skip quit)")
    args = ap.parse_args()

    # 1. import
    try:
        import matlab.engine
    except Exception:                                  # noqa: BLE001
        _install_hint()
        sys.exit(1)
    print("[ok] matlab.engine imported")

    # 2. start engine
    print("[..] starting MATLAB engine (first start can take ~15-30s)...", flush=True)
    t0 = time.time()
    try:
        eng = matlab.engine.start_matlab()
    except Exception as exc:                           # noqa: BLE001
        print(f"[FAIL] could not start MATLAB engine: {exc}")
        sys.exit(1)
    print(f"[ok] engine started in {time.time() - t0:.1f}s")

    try:
        try:
            print(f"[ok] MATLAB {eng.version(nargout=1)}")
        except Exception as exc:                       # noqa: BLE001
            print(f"[warn] could not read MATLAB version: {exc}")

        # 3. SimBiology license/availability
        try:
            has = eng.license("test", "SimBiology", nargout=1)
            ok = bool(has)
            print(f"[{'ok' if ok else 'FAIL'}] SimBiology license test -> {has}")
            if not ok:
                print("       SimBiology is not licensed for this MATLAB; stopping.")
                return
        except Exception as exc:                       # noqa: BLE001
            print(f"[warn] SimBiology license test errored (continuing): {exc}")

        # put our .m helpers on the MATLAB path
        here = os.path.dirname(os.path.abspath(__file__))
        eng.addpath(os.path.join(here, "matlab"), nargout=0)

        # 4. build + simulate a trivial model from code
        try:
            t, c = eng.sb_smoke(nargout=2)
            tv, cv = _flat(t), _flat(c)
            print(f"[ok] built + simulated a 1-cmt model: {len(tv)} timepoints, "
                  f"C(0)={cv[0]:.3f} -> C(end)={cv[-1]:.3f} over t=0..{tv[-1]:g}h")
        except Exception as exc:                       # noqa: BLE001
            print(f"[FAIL] sb_smoke build/simulate failed: {exc}")
            return

        # 5. optional: load + probe a real .sbproj (e.g. Vantage RA)
        if args.sbproj:
            if not os.path.exists(args.sbproj):
                print(f"[FAIL] --sbproj not found: {args.sbproj}")
            else:
                try:
                    info = eng.sb_project_info(args.sbproj, nargout=1)
                    print(f"[ok] loaded project '{info.get('name')}': "
                          f"{int(info.get('nSpecies', 0))} species, "
                          f"{int(info.get('nParameters', 0))} params, "
                          f"{int(info.get('nReactions', 0))} reactions, "
                          f"{int(info.get('nVariants', 0))} variants, "
                          f"{int(info.get('nDoses', 0))} doses; "
                          f"baseline sim ok={info.get('baselineOk')}")
                except Exception as exc:               # noqa: BLE001
                    print(f"[FAIL] sb_project_info failed: {exc}")

        print("\nToolchain smoke complete.")
    finally:
        if not args.keep_open:
            try:
                eng.quit()
            except Exception:                          # noqa: BLE001
                pass


if __name__ == "__main__":
    main()

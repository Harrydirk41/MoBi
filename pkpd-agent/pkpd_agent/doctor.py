"""Engine health-check: what's installed, what it unlocks, what's missing.

Run after each setup step to see the agent's real capabilities light up:

    python -m pkpd_agent.doctor
    python -m pkpd_agent.doctor --rscript "C:\\Program Files\\R\\R-4.4.1\\bin\\Rscript.exe"

It never raises - every probe is guarded - so it is safe to run on a fresh box.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from dataclasses import dataclass


@dataclass
class Probe:
    name: str
    ok: bool
    detail: str
    unlocks: str


def _py_pkg(name: str) -> tuple[bool, str]:
    try:
        mod = __import__(name)
        return True, getattr(mod, "__version__", "installed")
    except Exception:  # noqa: BLE001
        return False, "not importable"


def _run(cmd: list[str], timeout: int = 30) -> tuple[bool, str]:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=timeout)
        out = (p.stdout or p.stderr).strip().splitlines()
        return p.returncode == 0, (out[0] if out else "")
    except FileNotFoundError:
        return False, "executable not found"
    except Exception as exc:  # noqa: BLE001
        return False, f"{type(exc).__name__}: {exc}"


def _rscript_path(explicit: str | None) -> str | None:
    if explicit and os.path.exists(explicit):
        return explicit
    if explicit:
        return explicit  # let the probe report the failure
    return shutil.which("Rscript")


def run_probes(rscript: str | None = None) -> list[Probe]:
    probes: list[Probe] = []

    # --- Python-side ---
    for pkg, unlocks in [
        ("numpy", "real pkfit engine"),
        ("scipy", "real pkfit engine"),
        ("anthropic", "Claude-driven decisions (LLMPolicy)"),
        ("pharmpy", "population NLME / AMD (pharmpy tools)"),
    ]:
        ok, detail = _py_pkg(pkg)
        probes.append(Probe(f"python:{pkg}", ok, detail, unlocks))

    probes.append(Probe(
        "env:ANTHROPIC_API_KEY",
        bool(os.environ.get("ANTHROPIC_API_KEY")),
        "set" if os.environ.get("ANTHROPIC_API_KEY") else "unset",
        "real LLM runs",
    ))

    # --- R-side ---
    rs = _rscript_path(rscript)
    if not rs:
        probes.append(Probe("R:Rscript", False, "not on PATH (pass --rscript)",
                            "nlmixr2 + ospsuite backends"))
    else:
        ok, detail = _run([rs, "-e", "cat(R.version.string)"])
        probes.append(Probe("R:Rscript", ok, detail if ok else "found but failed to run",
                            "nlmixr2 + ospsuite backends"))
        if ok:
            for pkg, unlocks in [("nlmixr2", "real NLME popPK (nlmixr2 tools)"),
                                 ("ospsuite", "mechanistic PBPK/QSP (OSP tools)")]:
                pok, pdetail = _run([rs, "-e", f"suppressMessages(library({pkg}));"
                                              f"cat(as.character(packageVersion('{pkg}')))"])
                probes.append(Probe(f"R:{pkg}", pok,
                                    pdetail if pok else "not installed", unlocks))

    # --- .NET / OSP CLI ---
    ok, detail = _run(["dotnet", "--version"])
    probes.append(Probe(".NET runtime", ok, detail if ok else "not found",
                        "OSP native engine / MoBi.CLI"))

    return probes


def summary(probes: list[Probe]) -> str:
    width = max(len(p.name) for p in probes)
    lines = ["", "pkpd-agent engine health check", "=" * 60]
    for p in probes:
        mark = "OK " if p.ok else "-- "
        lines.append(f"[{mark}] {p.name:<{width}}  {p.detail}")
        lines.append(f"        └─ unlocks: {p.unlocks}")
    ok_n = sum(p.ok for p in probes)
    lines += ["=" * 60, f"{ok_n}/{len(probes)} checks passing.",
              "The mock-engine loop + real pkfit engine run with just numpy/scipy.",
              "See SETUP_WINDOWS.md to install the remaining backends.", ""]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description="pkpd-agent engine health check")
    ap.add_argument("--rscript", default=None,
                    help="Full path to the Rscript.exe that has nlmixr2/ospsuite")
    args = ap.parse_args()
    print(summary(run_probes(args.rscript)))


if __name__ == "__main__":
    main()

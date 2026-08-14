"""Configuration for the pkpd-agent decision loop."""

from __future__ import annotations

import os
from dataclasses import dataclass, field


@dataclass
class AgentConfig:
    """Runtime configuration.

    ``mock`` is the important switch: when True (the default), the engine
    adapters return synthetic-but-plausible results so the whole loop runs
    without pharmpy, an OSP/R install, or an Anthropic API key. Flip it off to
    drive the real engines.
    """

    # --- the LLM "brain" ---
    model: str = "claude-opus-5"
    effort: str = "medium"          # low | medium | high | xhigh | max
    max_tokens: int = 16000
    adaptive_thinking: bool | None = None  # None = auto (on for Claude-5, off for
                                           # Haiku/older which reject it)

    # --- loop control ---
    max_steps: int = 24             # hard stop on Decide->Act->Evaluate iterations
    stop_on_block: bool = False     # if True, a BLOCK verdict halts the loop instead
                                    # of being fed back to the model to fix

    # --- engines ---
    mock: bool = True               # synthetic engine results; no real deps required
    rscript_path: str = "Rscript"   # Rscript.exe that has nlmixr2 / ospsuite
    nlmixr2_est: str = "focei"      # nlmixr2 estimation method: focei | saem
    mobi_cli_path: str | None = None  # path to MoBi.CLI executable (non-mock, Windows)
    pksim_cli_path: str | None = field(
        default_factory=lambda: os.environ.get("PKPD_PKSIM_CLI"))
    # ^ PKSim.CLI.exe: builds a .pksim5 from a snapshot and runs it (headless PBPK)
    pksim_timeout_s: int = 900      # PK-Sim snap+export can take minutes
    stream_optimizer: bool = True   # print each optimizer evaluation (live progress)
    nonmem_available: bool = False  # whether pharmpy can reach a NONMEM install

    # --- provenance ---
    workdir: str = field(default_factory=lambda: os.getcwd())

    def anthropic_key_present(self) -> bool:
        return bool(os.environ.get("ANTHROPIC_API_KEY"))

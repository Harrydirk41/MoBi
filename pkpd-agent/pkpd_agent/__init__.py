"""pkpd-agent: an LLM decision loop over pharmacometric engines (pharmpy + OSP).

The package wires a Large Language Model into an explicit
Observe -> Decide -> Act -> Evaluate loop, where the model makes the modeling
decisions a human pharmacometrician would make and the heavy lifting is
delegated to trusted engines:

  * pharmpy  -> population PK/PD (NLME) estimation, AMD, VPC, bootstrap
  * OSP      -> mechanistic PBPK / QSP simulation (MoBi / PK-Sim)
  * NCA      -> non-compartmental analysis (gap-filling binding)

The engines guarantee *mechanical* correctness; a verification layer checks
each decision for *scientific* sanity before it is accepted.
"""

__version__ = "0.1.0"

from .config import AgentConfig
from .loop import DecisionLoop
from .state import ModelingSession

__all__ = ["AgentConfig", "DecisionLoop", "ModelingSession", "__version__"]

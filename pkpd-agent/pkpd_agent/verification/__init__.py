"""Verification layer: scientific sanity checks on each action's result.

The engines guarantee *mechanical* correctness (the fit ran, the ODE
integrated). These gates check *scientific* correctness before a result is
accepted - the difference between "plausible modeling" and "trustworthy
science". A gate returns zero or more Findings; a ``block`` finding is fed back
to the model as something it must address.
"""

from .gates import DEFAULT_GATES, run_gates

__all__ = ["DEFAULT_GATES", "run_gates"]

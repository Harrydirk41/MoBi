"""A transparent, non-LLM decision policy for population PK model building.

This is 'good decision making' you can run without an API key: a small expert
system that drives the *real* pkfit engine through a defensible workflow and
makes each choice from real statistics -

  * structural model selection by AIC (with a parsimony tie-break),
  * covariate inclusion by a likelihood-ratio test (drop OFV by >3.84 on 1 df),
  * model adequacy by a visual predictive check.

It reads prior results out of the session transcript, so each decision is
grounded in numbers actually computed the step before. It is interchangeable
with ``LLMPolicy`` - same loop, same tools, same verification gates.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .llm import ActStep, FinishStep, PolicyStep
from .state import ToolCall

LRT_CHI2_1DF_95 = 3.84


@dataclass
class PharmacometricPolicy:
    covariate_param: str = "CL"
    covariate: str = "WT"
    covariate_ref: float = 70.0

    _stage: int = field(default=0, init=False)
    _chosen_model: str | None = field(default=None, init=False)
    _final_id: str | None = field(default=None, init=False)
    _notes: dict[str, Any] = field(default_factory=dict, init=False)

    # ------------------------------------------------------------------ #
    def _fits(self, session) -> dict[str, dict[str, Any]]:
        return {
            o.content.get("fit_id"): o.content
            for o in session.observations
            if o.tool == "pkfit_fit" and o.ok
        }

    def _call(self, name: str, **args) -> ToolCall:
        return ToolCall(id=f"pm_{self._stage}", name=name, arguments=args)

    # ------------------------------------------------------------------ #
    def decide(self, session) -> PolicyStep:
        stage = self._stage
        self._stage += 1
        fits = self._fits(session)

        if stage == 0:
            return ActStep("OBSERVE: load the concentration-time data.",
                           [self._call("pkfit_load_data", source="builtin")])

        if stage == 1:
            return ActStep("OBSERVE: model-free NCA first pass for orientation.",
                           [self._call("pkfit_nca")])

        if stage == 2:
            return ActStep("ACT: fit the base one-compartment oral model.",
                           [self._call("pkfit_fit", model="1cpt_oral")])

        if stage == 3:
            return ActStep(
                "ACT: fit a two-compartment alternative to test whether the "
                "extra disposition phase is justified.",
                [self._call("pkfit_fit", model="2cpt_oral")])

        if stage == 4:
            base = fits.get("1cpt_oral")
            alt = fits.get("2cpt_oral")
            chosen, reason = self._choose_structural(base, alt)
            self._chosen_model = chosen
            return ActStep(
                f"DECIDE (structure): {reason} ACT: add {self.covariate}-on-"
                f"{self.covariate_param} to the {chosen} model and re-fit.",
                [self._call("pkfit_fit", model=chosen,
                            covariate_param=self.covariate_param,
                            covariate=self.covariate,
                            covariate_ref=self.covariate_ref)])

        if stage == 5:
            base = fits.get(self._chosen_model)
            cov_id = f"{self._chosen_model}+{self.covariate}_on_{self.covariate_param}"
            cov = fits.get(cov_id)
            final_id, reason = self._choose_covariate(base, cov, cov_id)
            self._final_id = final_id
            return ActStep(
                f"DECIDE (covariate): {reason} EVALUATE: run a VPC on the "
                f"final model ({final_id}).",
                [self._call("pkfit_vpc", fit_id=final_id)])

        if stage == 6:
            return FinishStep(self._summary(session, fits))

        return FinishStep("done")

    # ------------------------------------------------------------------ #
    def _choose_structural(self, base, alt) -> tuple[str, str]:
        if base is None:
            return "2cpt_oral", "base 1-cpt fit unavailable; falling back to 2-cpt."
        if alt is None or not alt.get("minimization_successful"):
            self._notes["structural"] = "2cpt unstable/failed -> keep 1cpt"
            return "1cpt_oral", (
                f"the 2-compartment fit did not minimize reliably; the "
                f"1-compartment model (AIC {base['aic']}) is retained on "
                "parsimony grounds.")
        d_aic = alt["aic"] - base["aic"]
        if d_aic < -2.0:   # meaningfully better
            self._notes["structural"] = f"2cpt better by dAIC={d_aic:.1f}"
            return "2cpt_oral", (
                f"the 2-compartment model lowers AIC by {-d_aic:.1f} "
                f"({alt['aic']} vs {base['aic']}) - the extra phase is justified.")
        self._notes["structural"] = f"1cpt kept (dAIC={d_aic:.1f})"
        return "1cpt_oral", (
            f"two compartments do not improve AIC ({alt['aic']} vs "
            f"{base['aic']}, dAIC={d_aic:+.1f}); the simpler 1-compartment "
            "model is preferred.")

    def _choose_covariate(self, base, cov, cov_id) -> tuple[str, str]:
        if base is None or cov is None:
            return (cov_id if cov else self._chosen_model), "covariate fit unavailable."
        d_ofv = base["ofv"] - cov["ofv"]
        self._notes["cov_dofv"] = round(d_ofv, 2)
        if cov.get("minimization_successful") and d_ofv > LRT_CHI2_1DF_95:
            coef = cov["parameter_estimates"].get(f"{self.covariate}_on_{self.covariate_param}")
            self._notes["covariate"] = f"kept ({self.covariate}_on_{self.covariate_param}={coef})"
            return cov_id, (
                f"the covariate drops OFV by {d_ofv:.1f} (> {LRT_CHI2_1DF_95}, "
                f"1 df, p<0.05); keep {self.covariate}-on-{self.covariate_param} "
                f"(coefficient {coef}).")
        self._notes["covariate"] = "dropped (not significant)"
        return self._chosen_model, (
            f"the covariate drops OFV by only {d_ofv:.1f} (<= {LRT_CHI2_1DF_95}); "
            "not significant, so it is dropped.")

    def _summary(self, session, fits) -> str:
        truth = session.get("truth") or {}
        final = fits.get(self._final_id, {})
        vpc = next((o.content for o in session.observations if o.tool == "pkfit_vpc"), {})
        est = final.get("parameter_estimates", {})
        lines = [
            "FINAL MODEL: " + str(self._final_id),
            f"  structural decision : {self._notes.get('structural','-')}",
            f"  covariate decision  : {self._notes.get('covariate','-')} "
            f"(LRT dOFV={self._notes.get('cov_dofv','-')})",
            f"  estimates           : {est}",
            f"  OFV/AIC             : {final.get('ofv','-')} / {final.get('aic','-')}",
            f"  VPC coverage        : {vpc.get('pct_observations_within_90_pi','-')}% within 90% PI",
        ]
        if truth:
            lines.append(
                f"  (validation vs known truth: CL_ref={truth.get('CL_ref')}, "
                f"V_ref={truth.get('V_ref')}, Ka={truth.get('Ka')}, "
                f"allometric_CL={truth.get('allometric_CL')})")
        lines.append(
            "  Note: naive-pooled estimation folds between-subject variability "
            "into residual error and biases the covariate exponent in small "
            "samples; a human should confirm with a full NLME fit and GOF plots.")
        return "\n".join(lines)

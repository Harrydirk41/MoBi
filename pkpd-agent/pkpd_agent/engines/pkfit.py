"""A real, self-contained PK estimation engine (numpy/scipy).

This is a genuine engine — it fits compartmental PK models to data by maximum
likelihood, computes standard errors from the observed information matrix, and
supports real model comparison (OFV / AIC / likelihood-ratio tests). It runs
in-process with no NONMEM / nlmixr2 / R.

Scope & honesty: estimation is **naive-pooled** (typical parameters, one
residual-error model) — not full nonlinear mixed-effects. Full NLME with random
effects is exactly the job that needs an external backend (the pharmpy/NONMEM
world). What this gives you is a real estimator the decision loop can drive end
to end, with real convergence flags, real precision, and real model selection.

Structural models (closed form):
  * 1cpt_oral : CL, V, Ka          (one compartment, first-order absorption)
  * 2cpt_oral : CL, V1, Q, V2, Ka  (two compartments, first-order absorption)

Residual error: proportional,  DV = pred * (1 + eps),  eps ~ N(0, sigma^2).
Covariate model (optional): param_i = param_pop * (cov_i / ref) ** coef,
with ``coef`` estimated (so a likelihood-ratio test on 1 df is meaningful).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from scipy import optimize

# --------------------------------------------------------------------------- #
# Structural models
# --------------------------------------------------------------------------- #

MODELS: dict[str, list[str]] = {
    "1cpt_oral": ["CL", "V", "Ka"],
    "2cpt_oral": ["CL", "V1", "Q", "V2", "Ka"],
}


def _pred_1cpt_oral(p: dict[str, float], dose: np.ndarray, t: np.ndarray) -> np.ndarray:
    CL, V, Ka = p["CL"], p["V"], p["Ka"]
    ke = CL / V
    denom = Ka - ke
    # guard the flip-flop / Ka≈ke singularity with a first-order limit
    if abs(denom) < 1e-6:
        return dose * ke / V * t * np.exp(-ke * t)
    return dose * Ka / (V * denom) * (np.exp(-ke * t) - np.exp(-Ka * t))


def _pred_2cpt_oral(p: dict[str, float], dose: np.ndarray, t: np.ndarray) -> np.ndarray:
    CL, V1, Q, V2, Ka = p["CL"], p["V1"], p["Q"], p["V2"], p["Ka"]
    k10 = CL / V1
    k12 = Q / V1
    k21 = Q / V2
    s = k10 + k12 + k21
    disc = math.sqrt(max(s * s - 4.0 * k10 * k21, 0.0))
    alpha = 0.5 * (s + disc)
    beta = 0.5 * (s - disc)
    eps = 1e-9

    def term(root: float) -> np.ndarray:
        others = [r for r in (alpha, beta, Ka) if r is not root]
        # coefficient for exp(-root t): (k21-root)/prod(other-root) with Ka scaling
        # Standard first-order-input two-compartment solution.
        num = (k21 - root)
        den = 1.0
        for r in others:
            den *= (r - root)
        return num / (den + eps) * np.exp(-root * t)

    coef = dose * Ka / V1
    return coef * (term(alpha) + term(beta) + term(Ka))


_PREDICTORS: dict[str, Callable[..., np.ndarray]] = {
    "1cpt_oral": _pred_1cpt_oral,
    "2cpt_oral": _pred_2cpt_oral,
}


# --------------------------------------------------------------------------- #
# Dataset
# --------------------------------------------------------------------------- #

@dataclass
class Dataset:
    """Row-aligned observation data (naive-pooled)."""
    subject: np.ndarray          # subject id per row
    time: np.ndarray             # hours
    dv: np.ndarray               # observed concentration
    dose: np.ndarray             # dose amount (mg) per row (subject's dose)
    covariates: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def n_obs(self) -> int:
        return int(self.dv.size)

    @property
    def n_subjects(self) -> int:
        return int(np.unique(self.subject).size)

    def summary(self) -> dict[str, Any]:
        return {
            "n_subjects": self.n_subjects,
            "n_obs": self.n_obs,
            "time_range_h": [float(self.time.min()), float(self.time.max())],
            "dv_range": [float(self.dv.min()), float(self.dv.max())],
            "covariates": {k: [float(v.min()), float(v.max())] for k, v in self.covariates.items()},
        }


# --------------------------------------------------------------------------- #
# Parameter packing (log scale for positive structural params)
# --------------------------------------------------------------------------- #

@dataclass
class _Spec:
    """Describes the free-parameter vector for one fit."""
    model: str
    struct_names: list[str]
    cov_param: str | None = None
    cov_name: str | None = None
    cov_ref: float = 1.0

    @property
    def names(self) -> list[str]:
        names = list(self.struct_names) + ["sigma_prop"]
        if self.cov_param:
            names.append(f"{self.cov_name}_on_{self.cov_param}")
        return names

    def unpack(self, theta: np.ndarray) -> tuple[dict[str, float], float, float | None]:
        n = len(self.struct_names)
        struct = {name: math.exp(theta[i]) for i, name in enumerate(self.struct_names)}
        sigma = math.exp(theta[n])
        coef = float(theta[n + 1]) if self.cov_param else None
        return struct, sigma, coef


# --------------------------------------------------------------------------- #
# Prediction with optional covariate
# --------------------------------------------------------------------------- #

def _predict(spec: _Spec, struct: dict[str, float], coef: float | None,
             data: Dataset) -> np.ndarray:
    predictor = _PREDICTORS[spec.model]
    if spec.cov_param and coef is not None:
        cov = data.covariates[spec.cov_name]
        scale = (cov / spec.cov_ref) ** coef
        rows = np.empty_like(data.time, dtype=float)
        # scale only the covariate-affected parameter, per row
        base = dict(struct)
        # vectorized over unique scale values would be faster; datasets are small
        for i in range(data.time.size):
            p = dict(base)
            p[spec.cov_param] = base[spec.cov_param] * scale[i]
            rows[i] = predictor(p, data.dose[i], np.array([data.time[i]]))[0]
        return rows
    return predictor(struct, data.dose, data.time)


# --------------------------------------------------------------------------- #
# Objective (-2 log likelihood, proportional error)
# --------------------------------------------------------------------------- #

def _ofv(theta: np.ndarray, spec: _Spec, data: Dataset) -> float:
    struct, sigma, coef = spec.unpack(theta)
    pred = _predict(spec, struct, coef, data)
    pred = np.clip(pred, 1e-8, None)
    sd = sigma * pred
    resid = data.dv - pred
    ll = -0.5 * (np.log(2 * np.pi * sd * sd) + (resid * resid) / (sd * sd))
    val = -2.0 * float(np.sum(ll))
    if not math.isfinite(val):
        return 1e12
    return val


def _num_hessian(f: Callable[[np.ndarray], float], x: np.ndarray,
                 eps: float = 1e-4) -> np.ndarray:
    n = x.size
    H = np.zeros((n, n))
    fx = f(x)
    for i in range(n):
        for j in range(i, n):
            xi = x.copy(); xi[i] += eps; xi[j] += eps
            xj = x.copy(); xj[i] += eps; xj[j] -= eps
            xk = x.copy(); xk[i] -= eps; xk[j] += eps
            xl = x.copy(); xl[i] -= eps; xl[j] -= eps
            H[i, j] = H[j, i] = (f(xi) - f(xj) - f(xk) + f(xl)) / (4 * eps * eps)
    return H


# --------------------------------------------------------------------------- #
# Public engine
# --------------------------------------------------------------------------- #

class PKFitEngine:
    """Real MLE estimation, NCA, and Monte-Carlo VPC."""

    # -- initial estimates from the data (so the optimizer starts sane) -- #
    def _initial_struct(self, model: str, data: Dataset) -> dict[str, float]:
        cmax = float(np.max(data.dv))
        dose = float(np.median(data.dose))
        # very rough NCA-flavored guesses
        auc = self._auc_pooled(data)
        cl = dose / max(auc, 1e-6)
        v = dose / max(cmax, 1e-6) * 2.0
        if model == "1cpt_oral":
            return {"CL": cl, "V": v, "Ka": 1.0}
        return {"CL": cl, "V1": v, "Q": cl, "V2": v, "Ka": 1.0}

    def _auc_pooled(self, data: Dataset) -> float:
        # AUC of the median profile across nominal times
        times = np.unique(data.time)
        med = np.array([np.median(data.dv[data.time == t]) for t in times])
        return float(np.trapezoid(med, times))

    # -- NCA ------------------------------------------------------------- #
    def nca(self, data: Dataset) -> dict[str, Any]:
        times = np.unique(data.time)
        med = np.array([np.median(data.dv[data.time == t]) for t in times])
        cmax = float(med.max())
        tmax = float(times[int(np.argmax(med))])
        auc = float(np.trapezoid(med, times))
        thalf = None
        if med.size >= 3 and med[-1] > 0 and med[-2] > med[-1]:
            k = (math.log(med[-3]) - math.log(med[-1])) / (times[-1] - times[-3])
            if k > 0:
                thalf = round(math.log(2) / k, 3)
        return {
            "c_max": round(cmax, 4), "t_max": round(tmax, 4),
            "auc": round(auc, 4), "t_half_terminal": thalf,
            "apparent_clearance_F": round(float(np.median(data.dose)) / max(auc, 1e-9), 4),
            "n_points": int(times.size), "source": "pkfit-nca",
        }

    # -- fit ------------------------------------------------------------- #
    def fit(self, data: Dataset, model: str = "1cpt_oral",
            covariate: dict[str, Any] | None = None) -> dict[str, Any]:
        if model not in MODELS:
            raise ValueError(f"unknown model {model!r}; choose {list(MODELS)}")
        struct_names = MODELS[model]
        spec = _Spec(model=model, struct_names=struct_names)
        if covariate:
            spec.cov_param = covariate["param"]
            spec.cov_name = covariate["cov"]
            spec.cov_ref = float(covariate.get("ref", 1.0))
            if spec.cov_param not in struct_names:
                raise ValueError(f"covariate parameter {spec.cov_param} not in {struct_names}")
            if spec.cov_name not in data.covariates:
                raise ValueError(f"covariate {spec.cov_name} not in dataset")

        init_struct = self._initial_struct(model, data)
        theta0 = [math.log(init_struct[n]) for n in struct_names] + [math.log(0.2)]
        if spec.cov_param:
            theta0.append(0.5)
        theta0 = np.array(theta0)

        obj = lambda th: _ofv(th, spec, data)  # noqa: E731
        # robust two-stage: Nelder-Mead to get close, then polish
        res = optimize.minimize(obj, theta0, method="Nelder-Mead",
                                options={"maxiter": 4000, "xatol": 1e-5, "fatol": 1e-5})
        res = optimize.minimize(obj, res.x, method="Nelder-Mead",
                                options={"maxiter": 2000, "xatol": 1e-6, "fatol": 1e-6})
        theta = res.x
        ofv = float(res.fun)

        struct, sigma, coef = spec.unpack(theta)
        estimates: dict[str, float] = {k: round(v, 5) for k, v in struct.items()}
        estimates["sigma_prop"] = round(sigma, 5)
        if spec.cov_param:
            estimates[f"{spec.cov_name}_on_{spec.cov_param}"] = round(coef, 5)

        # covariance from the observed information matrix
        rse: dict[str, float] = {}
        cond = float("inf")
        pos_def = False
        try:
            H = _num_hessian(obj, theta)
            cov = np.linalg.inv(0.5 * H)  # OFV = -2logL  ->  Fisher = 0.5*H
            eig = np.linalg.eigvalsh(0.5 * H)
            pos_def = bool(np.all(eig > 0))
            se = np.sqrt(np.clip(np.diag(cov), 0, None))
            cond = float(np.linalg.cond(cov))
            for i, name in enumerate(spec.names):
                # structural params are on log scale -> SE(log) ≈ RSE on natural scale
                if name in struct_names:
                    rse[name] = round(float(se[i]), 4)
                elif name == "sigma_prop":
                    rse[name] = round(float(se[i]), 4)
                else:  # covariate coefficient (natural scale)
                    rse[name] = round(float(se[i]) / (abs(coef) + 1e-9), 4)
        except np.linalg.LinAlgError:
            pass

        n_par = len(spec.names)
        n_obs = data.n_obs
        aic = ofv + 2 * n_par
        bic = ofv + n_par * math.log(n_obs)
        well_conditioned = math.isfinite(cond) and cond < 1e6
        successful = bool(
            res.success and pos_def and rse
            and all(math.isfinite(v) for v in rse.values())
            and well_conditioned
        )

        return {
            "model": model,
            "covariate": f"{spec.cov_name}_on_{spec.cov_param}" if spec.cov_param else None,
            "ofv": round(ofv, 3),
            "aic": round(aic, 3),
            "bic": round(bic, 3),
            "n_parameters": n_par,
            "n_obs": n_obs,
            "parameter_estimates": estimates,
            "relative_standard_errors": rse,
            "condition_number": round(cond, 1) if math.isfinite(cond) else None,
            "minimization_successful": successful,
            "source": "pkfit",
            "_theta": theta.tolist(),   # kept for VPC re-simulation
        }

    # -- VPC (Monte-Carlo) ---------------------------------------------- #
    def vpc(self, data: Dataset, fit: dict[str, Any], n_sim: int = 500) -> dict[str, Any]:
        spec = _Spec(model=fit["model"], struct_names=MODELS[fit["model"]])
        cov = fit.get("covariate")
        if cov:
            param, _, name = cov.partition("_on_")[0], None, cov.partition("_on_")[2]
            # covariate string is "<cov>_on_<param>"
            spec.cov_name = cov.split("_on_")[0]
            spec.cov_param = cov.split("_on_")[1]
            spec.cov_ref = self._infer_ref(data, spec.cov_name)
        theta = np.array(fit["_theta"])
        struct, sigma, coef = spec.unpack(theta)
        pred = _predict(spec, struct, coef, data)
        pred = np.clip(pred, 1e-8, None)
        rng = np.random.default_rng(12345)
        sims = pred[None, :] * (1 + rng.normal(0, sigma, size=(n_sim, pred.size)))
        lo = np.percentile(sims, 5, axis=0)
        hi = np.percentile(sims, 95, axis=0)
        inside = np.mean((data.dv >= lo) & (data.dv <= hi)) * 100.0
        return {
            "model": fit["model"],
            "pct_observations_within_90_pi": round(float(inside), 1),
            "n_sim": n_sim, "source": "pkfit-vpc",
        }

    @staticmethod
    def _infer_ref(data: Dataset, cov_name: str) -> float:
        return float(np.median(data.covariates[cov_name]))


# --------------------------------------------------------------------------- #
# Synthetic-but-real dataset (simulated from a KNOWN truth so recovery can be
# checked). 1-compartment oral truth with allometric WT on CL.
# --------------------------------------------------------------------------- #

def simulate_dataset(n_subjects: int = 12, seed: int = 7) -> tuple[Dataset, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    nominal_t = np.array([0.25, 0.5, 1.0, 2.0, 3.5, 5.0, 7.0, 9.0, 12.0, 24.0])
    # WT enters ONLY through CL (allometric 0.75); V and Ka are WT-independent,
    # so a WT-on-CL covariate model is correctly specified and should recover
    # the 0.75 exponent.
    truth = {"CL_ref": 3.5, "V_ref": 45.0, "Ka": 1.2, "WT_ref": 70.0,
             "allometric_CL": 0.75, "allometric_V": 0.0,
             "iiv_sd": 0.20, "prop_err": 0.15}

    subj, time, dv, dose, wt = [], [], [], [], []
    for s in range(1, n_subjects + 1):
        w = float(rng.uniform(55, 95))
        dose_amt = 320.0  # fixed mg dose (so WT enters only through allometry)
        cl_i = truth["CL_ref"] * (w / truth["WT_ref"]) ** truth["allometric_CL"] \
            * math.exp(rng.normal(0, truth["iiv_sd"]))
        v_i = truth["V_ref"] * (w / truth["WT_ref"]) ** truth["allometric_V"] \
            * math.exp(rng.normal(0, truth["iiv_sd"]))
        ka_i = truth["Ka"] * math.exp(rng.normal(0, truth["iiv_sd"]))
        p = {"CL": cl_i, "V": v_i, "Ka": ka_i}
        clean = _pred_1cpt_oral(p, np.full_like(nominal_t, dose_amt), nominal_t)
        noisy = clean * (1 + rng.normal(0, truth["prop_err"], size=nominal_t.size))
        noisy = np.clip(noisy, 1e-4, None)
        for t, c in zip(nominal_t, noisy):
            subj.append(s); time.append(t); dv.append(float(c))
            dose.append(dose_amt); wt.append(w)

    data = Dataset(
        subject=np.array(subj), time=np.array(time), dv=np.array(dv),
        dose=np.array(dose), covariates={"WT": np.array(wt)},
    )
    return data, truth

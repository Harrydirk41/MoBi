#!/usr/bin/env Rscript
# nlmixr2_fit.R  -  real NLME population PK fit, driven by the Python agent.
#
#   Rscript nlmixr2_fit.R <data.csv> <out.json> [model] [est] [covparam]
#
# data.csv is NONMEM-style: columns ID, TIME, DV, AMT, EVID, CMT
#   dosing rows  : EVID=1, AMT=<dose>, CMT=1 (depot), DV=0
#   observation  : EVID=0, DV=<conc>,  CMT=2 (central)
#   optional     : WTR = covariate/reference (baseline, per subject) - present
#                  when a covariate is requested; used as (WTR)^coef.
# model    : "1cpt_oral" (default) | "2cpt_oral"
# est      : "focei" (default) | "saem"
# covparam : "none" (default) | "CL" | "V"   -> allometric covariate on that
#            parameter, estimated inside the NLME fit.
#
# Emits JSON: { ok, model, covariate, ofv, aic, bic, n_obs,
#               parameter_estimates, relative_standard_errors,
#               iiv_cv_percent, shrinkage_percent,
#               pct_observations_within_90_pi (best-effort),
#               minimization_successful, source }
#
# NOTE: validate on your machine. The VPC coverage is a model-based
# Monte-Carlo VPC computed here (closed-form 1-cpt oral) using the estimated
# fixed effects + IIV + residual; it is guarded so a failure never blocks the
# fit.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: nlmixr2_fit.R <data.csv> <out.json> [model] [est] [covparam]")
data_csv   <- args[[1]]; out_json <- args[[2]]
model_name <- ifelse(length(args) >= 3, args[[3]], "1cpt_oral")
est_method <- ifelse(length(args) >= 4, args[[4]], "focei")
covparam   <- ifelse(length(args) >= 5, args[[5]], "none")

emit_error <- function(msg) {
  writeLines(jsonlite::toJSON(list(ok = FALSE, message = msg, source = "nlmixr2"),
                              auto_unbox = TRUE), out_json)
  quit(status = 0)
}

ok <- suppressWarnings(suppressMessages(
  requireNamespace("nlmixr2", quietly = TRUE) &&
  requireNamespace("jsonlite", quietly = TRUE)))
if (!ok) emit_error("nlmixr2 or jsonlite not installed")

suppressMessages({ library(nlmixr2); library(jsonlite) })
dat <- read.csv(data_csv)

# ---- model factory (optional allometric covariate on CL or V) -------------
mk_1cpt <- function(cov) {
  if (cov == "CL") {
    function() {
      ini({ tka <- log(1.0); tcl <- log(3.0); tv <- log(45)
            eta.cl ~ 0.09; eta.v ~ 0.09; cov_cl <- 0.75; prop.err <- 0.15 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl) * WTR^cov_cl
              v <- exp(tv + eta.v); linCmt() ~ prop(prop.err) })
    }
  } else if (cov == "V") {
    function() {
      ini({ tka <- log(1.0); tcl <- log(3.0); tv <- log(45)
            eta.cl ~ 0.09; eta.v ~ 0.09; cov_v <- 1.0; prop.err <- 0.15 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl)
              v <- exp(tv + eta.v) * WTR^cov_v; linCmt() ~ prop(prop.err) })
    }
  } else {
    function() {
      ini({ tka <- log(1.0); tcl <- log(3.0); tv <- log(45)
            eta.cl ~ 0.09; eta.v ~ 0.09; prop.err <- 0.15 })
      model({ ka <- exp(tka); cl <- exp(tcl + eta.cl)
              v <- exp(tv + eta.v); linCmt() ~ prop(prop.err) })
    }
  }
}

mk_2cpt <- function() {
  function() {
    ini({ tka <- log(1.0); tcl <- log(3.0); tv <- log(45)
          tq <- log(2.0); tv2 <- log(60)
          eta.cl ~ 0.09; eta.v ~ 0.09; prop.err <- 0.15 })
    model({ ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
            q <- exp(tq); v2 <- exp(tv2); linCmt() ~ prop(prop.err) })
  }
}

mdl <- if (model_name == "2cpt_oral") mk_2cpt() else mk_1cpt(covparam)

fit <- tryCatch(
  suppressWarnings(nlmixr2(mdl, dat, est = est_method, control = list(print = 0))),
  error = function(e) emit_error(paste("fit failed:", conditionMessage(e))))

# ---- fixed effects + precision -------------------------------------------
pf  <- as.data.frame(fit$parFixedDf)
est <- as.list(stats::setNames(pf$Estimate, rownames(pf)))
rse <- as.list(stats::setNames(pf[["%RSE"]] / 100, rownames(pf)))

objf <- tryCatch(fit$objDf$OBJF[1], error = function(e) NA)
aic  <- tryCatch(fit$objDf$AIC[1],  error = function(e) NA)
bic  <- tryCatch(fit$objDf$BIC[1],  error = function(e) NA)

# ---- IIV (CV%) and shrinkage (guarded) -----------------------------------
iiv <- tryCatch({
  om <- diag(as.matrix(fit$omega))
  as.list(stats::setNames(round(sqrt(exp(om) - 1) * 100, 1), names(om)))
}, error = function(e) NULL)

shrink <- tryCatch({
  sh <- fit$shrink
  if (is.null(sh)) NULL else as.list(round(as.numeric(sh[nrow(sh), ]), 1))
}, error = function(e) NULL)

# ---- model-based Monte-Carlo VPC (closed-form 1-cpt oral; guarded) --------
vpc_cov <- tryCatch({
  if (model_name != "1cpt_oral") stop("vpc only for 1cpt here")
  set.seed(20240101)
  ids   <- unique(dat$ID)
  om_cl <- as.matrix(fit$omega)["eta.cl", "eta.cl"]
  om_v  <- as.matrix(fit$omega)["eta.v",  "eta.v"]
  tcl <- est[["tcl"]]; tv <- est[["tv"]]; tka <- est[["tka"]]
  sig <- est[["prop.err"]]
  cov_coef <- if (covparam == "CL") est[["cov_cl"]]
              else if (covparam == "V") est[["cov_v"]] else 0
  pred1 <- function(dose, ka, ke, v, t) {
    d <- ka - ke
    if (abs(d) < 1e-6) dose * ke / v * t * exp(-ke * t)
    else dose * ka / (v * d) * (exp(-ke * t) - exp(-ka * t))
  }
  nsim <- 300
  times <- sort(unique(dat$TIME[dat$EVID == 0]))
  simvals <- setNames(vector("list", length(times)), as.character(times))
  for (s in seq_len(nsim)) {
    for (id in ids) {
      sub  <- dat[dat$ID == id, ]
      dose <- sub$AMT[sub$EVID == 1][1]
      wtr  <- if ("WTR" %in% names(sub)) sub$WTR[1] else 1
      ecl  <- rnorm(1, 0, sqrt(om_cl)); ev <- rnorm(1, 0, sqrt(om_v))
      cl <- exp(tcl + ecl) * (if (covparam == "CL") wtr^cov_coef else 1)
      v  <- exp(tv + ev)  * (if (covparam == "V")  wtr^cov_coef else 1)
      ka <- exp(tka); ke <- cl / v
      ot <- sub$TIME[sub$EVID == 0]
      p  <- pred1(dose, ka, ke, v, ot)
      p  <- pmax(p, 1e-8) * (1 + rnorm(length(p), 0, sig))
      for (j in seq_along(ot)) {
        k <- as.character(ot[j])
        simvals[[k]] <- c(simvals[[k]], p[j])
      }
    }
  }
  inside <- 0; total <- 0
  for (k in names(simvals)) {
    q <- quantile(simvals[[k]], c(0.05, 0.95), names = FALSE)
    obs <- dat$DV[dat$EVID == 0 & as.character(dat$TIME) == k]
    inside <- inside + sum(obs >= q[1] & obs <= q[2]); total <- total + length(obs)
  }
  round(100 * inside / total, 1)
}, error = function(e) NULL)

out <- list(
  ok = TRUE,
  model = model_name,
  covariate = if (covparam == "none") NULL else paste0("WT_on_", covparam),
  est = est_method,
  ofv = objf, aic = aic, bic = bic,
  n_obs = sum(dat$EVID == 0),
  parameter_estimates = est,
  relative_standard_errors = rse,
  iiv_cv_percent = iiv,
  shrinkage_percent = shrink,
  pct_observations_within_90_pi = vpc_cov,
  minimization_successful = isTRUE(all(is.finite(unlist(rse)))),
  source = "nlmixr2"
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, na = "null"), out_json)

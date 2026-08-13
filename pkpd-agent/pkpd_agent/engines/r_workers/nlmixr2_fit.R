#!/usr/bin/env Rscript
# nlmixr2_fit.R  -  real NLME population PK fit, driven by the Python agent.
#
#   Rscript nlmixr2_fit.R <data.csv> <out.json> [model] [est]
#
# data.csv is NONMEM-style: columns ID, TIME, DV, AMT, EVID, CMT
#   dosing rows  : EVID=1, AMT=<dose>, CMT=1 (depot), DV=0
#   observation  : EVID=0, DV=<conc>,  CMT=2 (central)
# model : "1cpt_oral" (default) | "2cpt_oral"
# est   : "focei" (default) | "saem"
#
# Emits JSON: { ok, ofv, aic, bic, n_obs, parameter_estimates,
#               relative_standard_errors, minimization_successful, source }
#
# NOTE: validate/tune this on your machine after installing nlmixr2. Model
# parameterization and result-field names follow nlmixr2 conventions but may
# need small edits for your data/columns.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: nlmixr2_fit.R <data.csv> <out.json> [model] [est]")
data_csv <- args[[1]]; out_json <- args[[2]]
model_name <- ifelse(length(args) >= 3, args[[3]], "1cpt_oral")
est_method <- ifelse(length(args) >= 4, args[[4]], "focei")

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

one_cmt_oral <- function() {
  ini({
    tka <- log(1.0); tcl <- log(3.0); tv <- log(45)
    eta.cl ~ 0.09; eta.v ~ 0.09
    prop.err <- 0.15
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv + eta.v)
    linCmt() ~ prop(prop.err)
  })
}

two_cmt_oral <- function() {
  ini({
    tka <- log(1.0); tcl <- log(3.0); tv <- log(45)
    tq  <- log(2.0); tv2 <- log(60)
    eta.cl ~ 0.09; eta.v ~ 0.09
    prop.err <- 0.15
  })
  model({
    ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    q  <- exp(tq);  v2 <- exp(tv2)
    linCmt() ~ prop(prop.err)
  })
}

mdl <- if (model_name == "2cpt_oral") two_cmt_oral else one_cmt_oral

fit <- tryCatch(
  suppressWarnings(nlmixr2(mdl, dat, est = est_method,
                          control = list(print = 0))),
  error = function(e) emit_error(paste("fit failed:", conditionMessage(e))))

pf <- as.data.frame(fit$parFixedDf)     # fixed-effect estimates + precision
est <- as.list(stats::setNames(pf$Estimate, rownames(pf)))
rse <- as.list(stats::setNames(pf[["%RSE"]] / 100, rownames(pf)))  # fraction
objf <- tryCatch(fit$objDf$OBJF[1], error = function(e) NA)
aic  <- tryCatch(fit$objDf$AIC[1],  error = function(e) NA)
bic  <- tryCatch(fit$objDf$BIC[1],  error = function(e) NA)

out <- list(
  ok = TRUE,
  model = model_name,
  est = est_method,
  ofv = objf, aic = aic, bic = bic,
  n_obs = sum(dat$EVID == 0),
  parameter_estimates = est,
  relative_standard_errors = rse,
  minimization_successful = isTRUE(all(is.finite(unlist(rse)))),
  source = "nlmixr2"
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, na = "null"), out_json)

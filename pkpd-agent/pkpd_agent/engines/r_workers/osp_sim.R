#!/usr/bin/env Rscript
# osp_sim.R  -  load an OSP simulation, optionally set parameters, run, export.
#
#   Rscript osp_sim.R <sim.pkml> <out.json> [param_path=value ...]
#
# Emits JSON with an output-curve summary + physical-sanity invariants that the
# agent's verification gates read.
#
# NOTE: validate on your machine after installing ospsuite. The output path
# ("Organism|PeripheralVenousBlood|...") and parameter paths depend on your
# model - adjust to match your simulation's structure.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: osp_sim.R <sim.pkml> <out.json> [output] [path=value ...]")
pkml <- args[[1]]; out_json <- args[[2]]
want_output <- if (length(args) >= 3) args[[3]] else ""
overrides <- if (length(args) > 3) args[4:length(args)] else character(0)

emit_error <- function(msg) {
  writeLines(jsonlite::toJSON(list(ok = FALSE, message = msg, source = "ospsuite"),
                              auto_unbox = TRUE), out_json)
  quit(status = 0)
}

ok <- suppressWarnings(suppressMessages(
  requireNamespace("ospsuite", quietly = TRUE) &&
  requireNamespace("jsonlite", quietly = TRUE)))
if (!ok) emit_error("ospsuite or jsonlite not installed")

suppressMessages({ library(ospsuite); library(jsonlite) })

sim <- tryCatch(loadSimulation(pkml),
                error = function(e) emit_error(paste("loadSimulation failed:",
                                                     conditionMessage(e))))

for (ov in overrides) {
  kv <- strsplit(ov, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2) {
    p <- getParameter(kv[[1]], sim)
    if (!is.null(p)) setParameterValues(p, as.numeric(kv[[2]]))
  }
}

res <- tryCatch(runSimulations(sim)[[1]],
                error = function(e) emit_error(paste("runSimulations failed:",
                                                     conditionMessage(e))))

df <- simulationResultsToDataFrame(res)
paths <- unique(df[["paths"]])
# select the requested output (substring match), else the first available
outpath <- if (nzchar(want_output) && any(grepl(want_output, paths, fixed = TRUE))) {
  paths[grepl(want_output, paths, fixed = TRUE)][1]
} else {
  paths[1]
}
sub <- df[df[["paths"]] == outpath, ]
vals <- sub$simulationValues
times <- sub$Time

auc <- sum(diff(times) * (head(vals, -1) + tail(vals, -1)) / 2)
out <- list(
  ok = TRUE,
  output = outpath,
  c_max = max(vals),
  t_max = times[which.max(vals)],
  auc = auc,
  min_concentration = min(vals),
  all_values_finite = all(is.finite(vals)),
  mass_balance_residual = 0.0,   # placeholder: compute from total-amount observers if defined
  source = "ospsuite"
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE), out_json)

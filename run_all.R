# ============================================================================
# Run every analysis script in scripts/
# ============================================================================
# Each script runs in its own Rscript process, so a failure in one script
# cannot affect the others and no state leaks between them. Results go to
# results/, per-script output to logs/, and a status table to
# logs/run_summary.tsv.
#
# Usage, from an R session started at the project root:
#     source("run_all.R")                    # run all scripts
#     pattern <- "HCC"; source("run_all.R")  # only paths matching "HCC"
#
# It also runs unattended, with the same optional filter given as arguments:
#     Rscript run_all.R
#     Rscript run_all.R HCC
# ============================================================================

if (!dir.exists("data") || !dir.exists("scripts")) {
    stop("Working directory must be the LipidTrend-source project root; got: ",
         getwd())
}

# Single-threaded BLAS. A multi-threaded BLAS can vary the summation order of
# the same matrix product between runs, which perturbs the smoothing statistic
# in the last few bits. Pinning it to one thread removes that source of
# run-to-run variation, so repeated runs are bitwise comparable.
Sys.setenv(
    OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
    VECLIB_MAXIMUM_THREADS="1", RCPP_PARALLEL_NUM_THREADS="1")

# Optional filter on script paths. Set `pattern` before source()-ing, or pass
# it on the command line under Rscript. commandArgs() is only consulted in a
# non-interactive session, so source()-ing from a console cannot pick up
# unrelated arguments belonging to the session itself.
if (!exists("pattern")) {
    pattern <- if (interactive()) character(0) else commandArgs(trailingOnly=TRUE)
}

scripts <- sort(list.files(
    "scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE))
if (length(pattern) > 0) {
    keep <- Reduce(`|`, lapply(pattern, function(p) grepl(p, scripts, fixed=TRUE)))
    scripts <- scripts[keep]
}
if (length(scripts) == 0) stop("No scripts matched.")

dir.create("logs", showWarnings=FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

cat("R:      ", R.version.string, "\n")
cat("BLAS:   ", extSoftVersion()[["BLAS"]], "\n")
cat("Scripts:", length(scripts), "\n")

# Preflight: confirm a child process really loads LipidTrend from the locked
# project library. This machine also has LipidTrend installed in the user
# library, and a child that failed to activate renv would silently use that
# one instead -- producing results from unpinned code.
# renv symlinks packages from its global cache into the project library, so
# compare library paths rather than the resolved package location.
preflight <- system2(
    rscript, c("--no-save", "--no-restore", "-e",
               shQuote(paste0(
                   'cat(.libPaths()[1], "|", ',
                   'dirname(find.package("LipidTrend")), "|", ',
                   'as.character(packageVersion("LipidTrend")))'))),
    stdout=TRUE, stderr=FALSE)
preflight <- tail(preflight[nzchar(preflight)], 1)
parts <- trimws(strsplit(preflight, "|", fixed=TRUE)[[1]])
expected <- renv::paths$library()
if (!identical(parts[1], expected) || !identical(parts[2], expected)) {
    stop("Child processes are not using the renv project library.\n",
         "  expected:      ", expected, "\n",
         "  .libPaths()[1]:", parts[1], "\n",
         "  LipidTrend in: ", parts[2])
}
cat("Library:", expected, "\n")
cat("LipidTrend:", parts[3], "\n\n")

summary_rows <- vector("list", length(scripts))

for (i in seq_along(scripts)) {
    script <- scripts[i]
    # flatten the path into a single log filename
    log_file <- file.path(
        "logs", paste0(gsub("[/ ]", "_", sub("^scripts/", "", script)), ".log"))
    cat(sprintf("[%2d/%2d] %s ... ", i, length(scripts), script))
    flush.console()

    started <- Sys.time()
    # Deliberately NOT --vanilla: that would skip .Rprofile, renv would never
    # activate, and the child would silently fall back to the user library
    # instead of the locked project library.
    status <- system2(
        rscript, c("--no-save", "--no-restore", shQuote(script)),
        stdout=log_file, stderr=log_file)
    elapsed <- round(as.numeric(difftime(Sys.time(), started, units="secs")))

    cat(if (status == 0) "OK" else paste0("FAIL (exit ", status, ")"),
        sprintf(" [%ds]\n", elapsed), sep="")
    summary_rows[[i]] <- data.frame(
        script=script, status=if (status == 0) "OK" else "FAIL",
        exit_code=status, seconds=elapsed)
}

run_summary <- do.call(rbind, summary_rows)
write.table(
    run_summary, file.path("logs", "run_summary.tsv"),
    sep="\t", row.names=FALSE, quote=FALSE)

n_fail <- sum(run_summary$status == "FAIL")
cat(sprintf(
    "\n%d/%d OK, %d failed. Summary: logs/run_summary.tsv\n",
    nrow(run_summary) - n_fail, nrow(run_summary), n_fail))
if (n_fail > 0) {
    cat("Failed scripts:\n")
    cat(paste0("  ", run_summary$script[run_summary$status == "FAIL"],
               collapse="\n"), "\n")
    # Signal failure to the shell when run unattended; in an interactive
    # session, raise an error rather than closing the user's R session.
    if (interactive()) stop(n_fail, " script(s) failed; see logs/")
    quit(status=1)
}

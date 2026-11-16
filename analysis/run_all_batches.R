# =============================================================================
# Scaling up: many trait pairs x many conditioning traits, in batches
#
# Set-up of the 8-CVD paper:
#   * 8 "core" traits: the CVD GWAS whose pairs are analysed
#   * 288 conditioning traits (risk factors, biomarkers, metabolites)
#   * 15 CVD pairs, each conditioned on every other CVD and on every
#     conditioning trait
#
# Steps
#   1. LDSC in batches: every batch = all core traits + up to 50 conditioning traits
#   2. per batch: every pair x every conditioning trait in that batch, both orientations
#      -> one CSV per batch and pair
#   3. combine all CSVs and apply the decision rules (summarise_conditional)
#
# Edit the CONFIG block, then run from the repository root:
#   Rscript analysis/run_all_batches.R        # everything, batch after batch
#   Rscript analysis/run_all_batches.R 3      # only batch 3 (e.g. one cluster job per batch)
# Finished batches and pairs are skipped when you re-run, so the script can be
# restarted after an interruption. After per-batch jobs, run it once without a
# number to combine and summarise.
# =============================================================================
source("R/conditional_rg.R")

# one thread per process: parallelism comes from running models side by side
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
only_batch <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
if (is.na(only_batch)) only_batch <- NULL

# ---- CONFIG -----------------------------------------------------------------
if (!exists("cfg")) cfg <- list(
  # traits whose pairs you analyse (included in every LDSC batch)
  core_files = Sys.glob("sumstats/core/*.sumstats.gz"),
  # conditioning traits (risk factors, metabolites, ...)
  extra_files = Sys.glob("sumstats/conditioning/*.sumstats.gz"),
  ld = "ld/eur_w_ld_chr/",
  # pairs as "A~B", using trait names (file names without .sumstats.gz)
  pairs = c("UKB_AAA~UKB_VTE", "UKB_VTE~UKB_IHD", "UKB_VTE~UKB_AF", "UKB_VTE~UKB_HF", "UKB_VTE~UKB_VV",
            "UKB_AS~UKB_IHD", "UKB_AS~UKB_HF", "UKB_VV~UKB_HF", "UKB_VV~UKB_AF", "UKB_IHD~UKB_AF",
            "UKB_IHD~UKB_HF", "UKB_IHD~UKB_CD", "UKB_HF~UKB_AF", "UKB_HF~UKB_CD", "UKB_AF~UKB_CD"),
  # core traits are also used as conditioning traits for the pairs they are not part of
  core_as_conditioning = TRUE,
  alpha = 0.05,                                 # significance threshold for conditioning_p
  batch_size = 50,                              # conditioning traits per LDSC batch
  workers = max(1, parallel::detectCores() - 1),  # parallel model fits
  out_dir = "results"
)

# ---- 1. LDSC in batches ---------------------------------------------------------
ldsc_files <- run_ldsc_batches(core = cfg$core_files, extra = cfg$extra_files,
                               ld = cfg$ld, batch_size = cfg$batch_size,
                               out_dir = file.path(cfg$out_dir, "ldsc_batches"), only = only_batch)

core_names <- trait_names_from_files(cfg$core_files)
pair_table <- do.call(rbind, lapply(strsplit(cfg$pairs, "~", fixed = TRUE),
                                    function(ab) data.frame(A = ab[1], B = ab[2])))
missing <- setdiff(c(pair_table$A, pair_table$B), core_names)
if (length(missing)) stop("Pair traits not among the core traits: ", paste(missing, collapse = ", "))
message(nrow(pair_table), " pairs: ", paste0(pair_table$A, "~", pair_table$B, collapse = ", "))

# ---- 2. conditional models, batch by batch ---------------------------------------
pair_dir <- file.path(cfg$out_dir, "conditional")
dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

for (b in if (is.null(only_batch)) seq_along(ldsc_files) else only_batch) {
  LDSC  <- load_ldsc(ldsc_files[b])
  extra_here <- setdiff(colnames(LDSC$S), core_names)
  for (i in seq_len(nrow(pair_table))) {
    A <- pair_table$A[i]; B <- pair_table$B[i]
    out_csv <- file.path(pair_dir, sprintf("batch%02d_%s~%s.csv", b, A, B))
    if (file.exists(out_csv)) next
    # conditioning traits: the other core traits (once, in batch 1, as they are
    # in every batch) + the conditioning traits of this batch
    other_core <- if (cfg$core_as_conditioning && b == 1) setdiff(core_names, c(A, B))
    cond <- c(other_core, extra_here)
    message(sprintf("[batch %d/%d] %s~%s: %d conditioning traits", b, length(ldsc_files), A, B, length(cond)))
    res <- conditional_rg(LDSC, A, B, cond, workers = cfg$workers)
    if (!is.null(res)) write.csv(cbind(batch = b, res), out_csv, row.names = FALSE)
  }
}

# ---- 3. combine and decide ----------------------------------------------------------
if (!is.null(only_batch)) {
  message("Batch ", only_batch, " done. Run without a batch number to combine all batches.")
  quit(save = "no")
}
all_res <- do.call(rbind, lapply(list.files(pair_dir, "\\.csv$", full.names = TRUE), utils::read.csv))
summ <- summarise_conditional(all_res, alpha = cfg$alpha)
write.csv(summ$all, file.path(cfg$out_dir, "all_conditional_results.csv"), row.names = FALSE)
write.csv(summ$significant, file.path(cfg$out_dir, "significant_conditioning_traits.csv"), row.names = FALSE)
message(sprintf("Done: %d tests, %d significant pair x conditioning trait combinations (%s)",
                nrow(summ$all), NROW(summ$significant), cfg$out_dir))

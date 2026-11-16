# =============================================================================
# Conditional genetic correlations with Genomic SEM
# -----------------------------------------------------------------------------
# Question: does the genetic correlation between traits A and B change when we
# condition on the genetic liability of a third trait C (the conditioning trait)?
#
#   marginal rg        rg(A,B), from multivariable LDSC
#   direct             the A-B path when C is in the model (conditional rg)
#   indirect           marginal - direct = (C ~ B path) x (A ~ C path)
#   conditioning p     p-value of the A ~ C path in the model with both B and C.
#                      When C is genetically correlated with B, this is the test
#                      of marginal = direct (Romero et al. 2025, Biol Psychiatry,
#                      Supplementary Note 7).
#
# The model is fitted in both orientations (A on B, and B on A) because the
# direction of the A-B relationship is unknown. The smaller indirect effect of
# the two is reported as a conservative estimate.
#
# Source this file, then see example/run_example.R and
# analysis/run_all_batches.R.
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicSEM)
})

# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

#' Trait names from munged file names: "path/UKB_AF.sumstats.gz" -> "UKB_AF"
trait_names_from_files <- function(files) {
  sanitize_names(sub("\\.sumstats(\\.gz)?$", "", basename(files)))
}

#' lavaan-safe names (no "-" or other special characters)
sanitize_names <- function(x) {
  x <- gsub("-", "_", x, fixed = TRUE)
  gsub("[^A-Za-z0-9_.]", "_", x)
}

stop_if_missing <- function(paths, label = "input") {
  miss <- paths[!file.exists(paths)]
  if (length(miss)) stop(sprintf("Missing %s files:\n- %s", label, paste(miss, collapse = "\n- ")))
  invisible(TRUE)
}

#' Load an LDSC result saved with saveRDS() (.rds) or save() (.RData / .R)
load_ldsc <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj)) {
    e <- new.env(parent = emptyenv())
    objs <- load(path, envir = e)
    obj <- get(if ("res" %in% objs) "res" else objs[1], envir = e)
  }
  stopifnot(is.list(obj), "S" %in% names(obj))
  obj
}

#' Genetic correlation matrix from an LDSC object
rg_matrix <- function(LDSC) {
  r <- suppressWarnings(stats::cov2cor(LDSC$S))
  dimnames(r) <- list(colnames(LDSC$S), colnames(LDSC$S))
  r
}

# -----------------------------------------------------------------------------
# Step 1: multivariable LDSC
# -----------------------------------------------------------------------------

#' Run GenomicSEM::ldsc on munged summary statistics and save the result.
#'
#' @param files       munged .sumstats.gz files (GenomicSEM::munge output)
#' @param ld, wld     folder with LD scores (e.g. eur_w_ld_chr/)
#' @param sample.prev, population.prev  NA gives observed-scale estimates. All
#'                    estimates used here are standardised, so the scale does
#'                    not change them.
#' @param out_file    where to save the LDSC object (.rds)
run_ldsc <- function(files, ld, wld = ld, trait.names = trait_names_from_files(files),
                     sample.prev = NA, population.prev = NA, out_file = NULL, ...) {
  stop_if_missing(files, "summary statistics")
  k <- length(files)
  args <- list(traits = files, ld = ld, wld = wld,
               sample.prev = rep_len(sample.prev, k), population.prev = rep_len(population.prev, k),
               trait.names = trait.names, ...)
  # write GenomicSEM's log next to the saved object instead of the working directory
  if (!is.null(out_file) && is.null(args$ldsc.log)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    args$ldsc.log <- sub("\\.rds$", "", out_file)
  }
  res <- do.call(ldsc, args)
  colnames(res$S) <- rownames(res$S) <- sanitize_names(colnames(res$S))
  if (!is.null(out_file)) {
    saveRDS(res, out_file)
    message("Saved LDSC object: ", out_file)
  }
  invisible(res)
}

#' Many conditioning traits: run LDSC in batches.
#'
#' Every batch contains all `core` traits (the traits whose pairs you analyse)
#' plus up to `batch_size` of the `extra` traits (conditioning traits). The
#' sampling covariance matrix of LDSC grows with the 4th power of the number of
#' traits, so one run with hundreds of traits is not feasible; batches of ~50
#' conditioning traits keep each run small.
#'
#' @param only   run only these batch numbers (e.g. one batch per cluster job)
#' @return paths of the saved LDSC objects, one per batch
run_ldsc_batches <- function(core, extra, ld, wld = ld, batch_size = 50L,
                             out_dir = "ldsc_batches", prefix = "ldsc_batch_", only = NULL, ...) {
  stop_if_missing(c(core, extra), "summary statistics")
  extra <- sort(setdiff(extra, core))
  batches <- split(extra, ceiling(seq_along(extra) / batch_size))
  message(sprintf("%d core + %d conditioning traits -> %d LDSC batches",
                  length(core), length(extra), length(batches)))
  out <- character(length(batches))
  for (i in seq_along(batches)) {
    out[i] <- file.path(out_dir, sprintf("%s%02d.rds", prefix, i))
    if (!is.null(only) && !i %in% only) next
    if (file.exists(out[i])) { message("Exists, skipping: ", out[i]); next }
    message(sprintf("Batch %d/%d: %d traits", i, length(batches), length(core) + length(batches[[i]])))
    run_ldsc(c(core, batches[[i]]), ld = ld, wld = wld, out_file = out[i], ...)
  }
  out
}

# -----------------------------------------------------------------------------
# Step 2: conditional models
# -----------------------------------------------------------------------------

# one regression path (standardised) from usermodel() output
.path <- function(res, lhs, rhs) {
  r <- res[res$op == "~" & res$lhs == lhs & res$rhs == rhs, , drop = FALSE]
  if (nrow(r) != 1) return(c(est = NA, se = NA, p = NA))
  c(est = as.numeric(r$STD_Genotype), se = as.numeric(r$STD_Genotype_SE), p = as.numeric(r$p_value))
}

#' Fit one orientation: outcome ~ predictor + C, C ~ predictor
#'
#' Output columns:
#'   direct                      outcome ~ predictor, adjusted for C (conditional rg)
#'   conditioning_effect, _p     outcome ~ C, adjusted for the predictor.
#'                               conditioning_p tests whether direct differs from
#'                               the marginal rg
#'   conditioning_predictor_effect, _p   C ~ predictor (must be non-zero for the test)
#'   indirect                    conditioning_predictor_effect x conditioning_effect
#'                               (= marginal rg - direct)
#'   total                       direct + indirect (= marginal rg)
#'   ratio                       indirect / total: proportion of the rg shared with C
fit_orientation <- function(LDSC, outcome, predictor, cond, estimation = "DWLS") {
  model <- paste0(outcome, " ~ ", predictor, " + ", cond, "\n", cond, " ~ ", predictor)
  out <- tryCatch({
    utils::capture.output(fit <- suppressWarnings(suppressMessages(
      usermodel(LDSC, estimation = estimation, model = model, CFIcalc = FALSE))))
    fit$results
  }, error = function(e) {
    message(sprintf("[ERR] %s ~ %s + %s: %s", outcome, predictor, cond, conditionMessage(e)))
    NULL
  })
  if (is.null(out)) return(NULL)
  d  <- .path(out, outcome, predictor)   # direct (conditional) path
  oc <- .path(out, outcome, cond)        # conditioning trait -> outcome
  cp <- .path(out, cond, predictor)      # predictor -> conditioning trait
  indirect    <- unname(cp["est"] * oc["est"])
  indirect_se <- unname(sqrt(oc["est"]^2 * cp["se"]^2 + cp["est"]^2 * oc["se"]^2))   # delta method
  total <- unname(d["est"] + indirect)
  data.frame(
    outcome = outcome, predictor = predictor, conditioning_trait = cond,
    direct = d["est"], direct_se = d["se"], direct_p = d["p"],
    conditioning_effect = oc["est"], conditioning_se = oc["se"], conditioning_p = oc["p"],
    conditioning_predictor_effect = cp["est"], conditioning_predictor_se = cp["se"],
    conditioning_predictor_p = cp["p"],
    indirect = indirect, indirect_se = indirect_se,
    total = total, ratio = indirect / total,
    row.names = NULL
  )
}

#' Conditional rg of pair A-B for one or more conditioning traits, in both orientations.
#'
#' @param LDSC                 output of run_ldsc()/GenomicSEM::ldsc (or a path to it)
#' @param A, B                 trait names (as in colnames(LDSC$S)); pair is written "A~B"
#' @param conditioning_traits  one or more conditioning traits
#' @param workers              >1 fits conditioning traits in parallel (future.apply)
#' @return one row per conditioning trait x orientation. Orientation "A_on_B" has
#'         A as outcome and B as predictor; "B_on_A" the reverse.
conditional_rg <- function(LDSC, A, B, conditioning_traits, workers = 1, estimation = "DWLS") {
  if (is.character(LDSC)) LDSC <- load_ldsc(LDSC)
  traits <- colnames(LDSC$S)
  need <- c(A, B, conditioning_traits)
  if (any(!need %in% traits)) stop("Not in LDSC object: ", paste(setdiff(need, traits), collapse = ", "))
  rg <- rg_matrix(LDSC)
  one <- function(cond) {
    rbind(cbind(direction = "A_on_B", fit_orientation(LDSC, A, B, cond, estimation)),
          cbind(direction = "B_on_A", fit_orientation(LDSC, B, A, cond, estimation)))
  }
  res <- if (workers > 1) {
    oplan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
    options(future.globals.maxSize = 8 * 1024^3)
    future.apply::future_lapply(conditioning_traits, one, future.seed = TRUE)
  } else lapply(conditioning_traits, one)
  res <- do.call(rbind, res)
  if (is.null(res)) return(NULL)
  res$pair <- paste0(A, "~", B)
  res$marginal_rg <- rg[A, B]
  first <- c("pair", "direction", "conditioning_trait", "marginal_rg")
  res[, c(first, setdiff(names(res), first))]
}

# -----------------------------------------------------------------------------
# Step 3: decide which conditioning traits change a pair's genetic correlation
# -----------------------------------------------------------------------------

#' Apply the decision rules used in the papers.
#'
#'  1. Eligible conditioning trait: genetically correlated with both traits of
#'     the pair, i.e. conditioning_predictor_p < 0.01 in both orientations.
#'     Without this, conditioning_p does not test the change in rg
#'     (Supplementary Note 7).
#'  2. Significant: conditioning_p < alpha in BOTH orientations.
#'  3. Report the orientation with the smaller indirect effect (conservative
#'     minimum, as the direction between A and B is unknown).
#'  4. Keep conditioning traits whose overlap reduces the correlation
#'     (ratio > 0). A negative ratio means conditioning makes the A-B rg
#'     stronger (a suppressor), not weaker.
#'
#' @param alpha  significance threshold for conditioning_p
#' @return list(all = every test with eligibility and significance,
#'              significant = one row per pair x significant conditioning trait)
summarise_conditional <- function(res, alpha = 0.05, predictor_p = 0.01, positive_only = TRUE) {
  res <- as.data.frame(res)
  key <- paste(res$pair, res$conditioning_trait)
  eligible <- tapply(res$conditioning_predictor_p < predictor_p, key,
                     function(x) length(x) == 2 && all(x %in% TRUE))
  res$eligible <- as.vector(eligible[key])
  res$significant <- res$eligible & !is.na(res$conditioning_p) & res$conditioning_p < alpha
  both <- tapply(res$significant, key, function(x) length(x) == 2 && all(x))
  cand <- res[as.vector(both[key]), , drop = FALSE]
  sig <- do.call(rbind, lapply(split(cand, paste(cand$pair, cand$conditioning_trait)),
                               function(d) d[which.min(d$indirect), , drop = FALSE]))
  if (!is.null(sig)) {
    sig$suppressor <- sig$ratio < 0
    if (positive_only) sig <- sig[!sig$suppressor, , drop = FALSE]
    sig <- sig[order(sig$pair, -sig$ratio), , drop = FALSE]
    rownames(sig) <- NULL
  }
  list(all = res, significant = sig)
}

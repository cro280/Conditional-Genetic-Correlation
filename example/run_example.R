# =============================================================================
# Test run: does the genetic overlap between IHD and AF also overlap with
#   (1) heart failure (HF)?
#   (2) height?
#
# Run from the repository root:   Rscript example/run_example.R
# Takes about a minute (most of it is LDSC).
# =============================================================================
source("R/conditional_rg.R")

# ---- inputs -----------------------------------------------------------------
# Munged summary statistics (GenomicSEM::munge output, one file per trait)
files <- c("example_data/UKB_IHD.sumstats.gz",
           "example_data/UKB_AF.sumstats.gz",
           "example_data/UKB_HF.sumstats.gz",
           "example_data/height.sumstats.gz")

# European LD scores (eur_w_ld_chr/, see README)
ld <- "ld/eur_w_ld_chr/"

dir.create("example_output", showWarnings = FALSE)
options(width = 200)
show <- function(x, cols, digits = 3) print(x[, cols], digits = digits, row.names = FALSE)

# ---- step 1: multivariable LDSC on all four traits ----------------------------
ldsc_file <- "example_output/ldsc_IHD_AF_HF_height.rds"
LDSC <- if (file.exists(ldsc_file)) load_ldsc(ldsc_file) else
  run_ldsc(files, ld = ld, out_file = ldsc_file)

cat("\n== Step 1. Marginal genetic correlations ==\n")
print(round(rg_matrix(LDSC), 2))

# ---- step 2: IHD-AF conditioned on HF and on height, both orientations -------
res <- conditional_rg(LDSC, A = "UKB_IHD", B = "UKB_AF", conditioning_traits = c("UKB_HF", "height"))
write.csv(res, "example_output/IHD_AF_conditional.csv", row.names = FALSE)

cat("\n== Step 2. Conditional models (A_on_B: IHD ~ AF + C, C ~ AF;  B_on_A: AF ~ IHD + C, C ~ IHD) ==\n")
show(res, c("direction", "conditioning_trait", "marginal_rg", "direct", "direct_se", "direct_p",
            "conditioning_p", "conditioning_predictor_p"))
show(res, c("direction", "conditioning_trait", "indirect", "indirect_se", "total", "ratio"))

# ---- step 3: decisions --------------------------------------------------------
s <- summarise_conditional(res, alpha = 0.05)
write.csv(s$all, "example_output/IHD_AF_conditional_flagged.csv", row.names = FALSE)

cat("\n== Step 3. Eligible (correlated with both traits) and significant ==\n")
show(s$all, c("direction", "conditioning_trait", "conditioning_predictor_p", "eligible",
              "conditioning_p", "significant"))

cat("\n== Result: significant conditioning traits, conservative orientation (smallest indirect effect) ==\n")
if (!is.null(s$significant) && nrow(s$significant)) {
  show(s$significant, c("pair", "conditioning_trait", "direction", "marginal_rg", "direct", "direct_p",
                        "indirect", "ratio", "conditioning_p"))
  write.csv(s$significant, "example_output/IHD_AF_significant.csv", row.names = FALSE)
}
excluded <- setdiff(unique(res$conditioning_trait), s$significant$conditioning_trait)
if (length(excluded)) cat("\nNot reported:", paste(excluded, collapse = ", "),
                         "(see README for why)\n")

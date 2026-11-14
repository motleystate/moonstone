# Reference for tests/analysis/test_permanova_model.py: vegan adonis2 and betadisper on synthetic data.
# Run from this folder: Rscript permanova_adonis2.R   (R and vegan from Ubuntu: apt install r-base r-cran-vegan)
# Writes the synthetic inputs (permanova_counts.csv, permanova_metadata.csv) and permanova_adonis2_expected.csv.
library(vegan)
set.seed(2026)
n <- 60
meta <- data.frame(
  group = rep(c("a", "b", "c"), each = 20),
  sex = sample(c("f", "m"), n, replace = TRUE),
  age = round(rnorm(n, 35, 8)),
  bmi = round(rnorm(n, 25, 4), 1),
  batch = rep(c("b1", "b2", "b3", "b4"), times = n / 4),
  row.names = sprintf("S%02d", 1:n)
)
# 12 taxa; group a/b/c, age and batch shift the abundances
effect <- cbind(group = rep(c(0, 1, 2), each = 20), age = (meta$age - 35) / 8)
logmean <- outer(seq(6, 1, length.out = 12), rep(1, n)) +
  outer(rnorm(12), effect[, "group"]) * 0.5 + outer(rnorm(12), effect[, "age"]) * 0.3 +
  matrix(rnorm(12 * n, 0, 0.6), 12, n)
counts <- sapply(1:n, function(j) rmultinom(1, 5000, exp(logmean[, j])))
dimnames(counts) <- list(sprintf("taxon%02d", 1:12), rownames(meta))
write.csv(counts, "permanova_counts.csv")
write.csv(meta, "permanova_metadata.csv")

d <- vegdist(t(counts), method = "bray")
rows <- list()
for (strata in c("none", "batch")) {
  perm <- if (strata == "none") how(nperm = 9999) else how(nperm = 9999, blocks = meta$batch)
  for (by in c("terms", "margin")) {
    set.seed(1)
    fit <- adonis2(d ~ group + sex + age + bmi, data = meta, by = by, permutations = perm)
    fit <- fit[!rownames(fit) %in% c("Residual", "Total"), ]
    rows[[length(rows) + 1]] <- data.frame(
      strata = strata, by = by, term = rownames(fit), df = fit$Df, R2 = fit$R2, F = fit$F, p = fit[["Pr(>F)"]])
  }
}
write.csv(do.call(rbind, rows), "permanova_adonis2_expected.csv", row.names = FALSE)

# betadisper + permutest, one factor at a time
disp <- list()
for (term in c("group", "sex")) {
  set.seed(1)
  bd <- betadisper(d, meta[[term]], type = "median")
  pt <- permutest(bd, permutations = 9999)$tab
  disp[[term]] <- data.frame(term = term, F = pt$F[1], p = pt[["Pr(>F)"]][1])
}
write.csv(do.call(rbind, disp), "permanova_betadisper_expected.csv", row.names = FALSE)

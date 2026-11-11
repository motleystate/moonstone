# Reference for tests/analysis/test_association.py: R lm() on CLR values of synthetic data.
# The counts have no zeros, so the CLR needs no zero replacement and is the same in R and Python.
# Run from this folder: Rscript association_lm.R   (R 4.3.3, base packages only)
counts <- as.matrix(read.csv("association_counts.csv", row.names = 1, check.names = FALSE))
metadata <- read.csv("association_metadata.csv", row.names = 1)
metadata <- metadata[colnames(counts), ]
metadata$disease <- relevel(factor(metadata$disease), ref = "healthy")
metadata$smoker <- relevel(factor(metadata$smoker), ref = "no")
logs <- log(counts)
clr <- sweep(logs, 2, colMeans(logs))
rows <- list()
for (feature in rownames(clr)) {
  fit <- lm(clr[feature, ] ~ disease + age + smoker, data = metadata)
  table <- summary(fit)$coefficients
  ci <- confint(fit)
  for (term in c("diseasebv", "diseasecandida", "age")) {
    rows[[length(rows) + 1]] <- data.frame(
      feature = feature, term = term, coef = table[term, 1], std_error = table[term, 2],
      t_value = table[term, 3], p_value = table[term, 4], ci_low = ci[term, 1], ci_high = ci[term, 2]
    )
  }
}
write.csv(do.call(rbind, rows), "association_lm_expected.csv", row.names = FALSE)

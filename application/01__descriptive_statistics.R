rm(list = ls())
library(zanicc)
library(xtable)

FF <- function(x, digits = 3, width = 2) formatC(x, digits = digits, width = width, format = "f")
format_numbers <- function(x1, x2) paste0(FF(x1), " (", FF(x2), ")")

local_dir <- "./application/data"

# Import modern and fossil pollen counts
data("pollen_climate", package = "zanicc")
Y <- pollen_climate$Y
Y_fossil <- readRDS(file = file.path(local_dir, "Y_fossil.rds"))

taxa_names <- colnames(Y)
all.equal(colnames(Y_fossil), taxa_names)

# Modern
tab <- data.frame(
  mean = format_numbers(colMeans(Y), colMeans(Y_fossil)),
  di = format_numbers(apply(Y, 2, var) / colMeans(Y),
                      apply(Y_fossil, 2, var) / colMeans(Y_fossil)),
  prop0 = format_numbers(100*colMeans(Y == 0), 100*colMeans(Y_fossil == 0)),
  zi = format_numbers(apply(Y, 2, zi_binomial, N = rowSums(Y)),
                      apply(Y_fossil, 2, zi_binomial, N = rowSums(Y_fossil))))
rownames(tab) <- taxa_names
tab <- tab[sort(taxa_names), ]
tab$taxa <- paste0("\\textit{", rownames(tab), "}")
tab <- tab[, c("taxa", "mean", "di", "prop0", "zi")]
colnames(tab) <- c("Taxa", "Abundance", "$\\operatorname{DI}\\lbrack Y_j\\rbrack$",
                   "\\% zeros", "$\\operatorname{ZI}\\lbrack Y_j\\rbrack$")
print(xtable(tab, digits = c(rep(0, 6))), include.rownames = FALSE,
      sanitize.text.function = force)

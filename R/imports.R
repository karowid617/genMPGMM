# Package-level imports and global variable declarations.
#
# stats functions used directly (rnorm, median) must be declared so
# R CMD check does not flag them as undefined globals.
#
# .data is the rlang/ggplot2 tidy-evaluation pronoun; it is declared via
# globalVariables so check does not flag it as an unbound variable in
# plot helpers that live behind a requireNamespace() guard.

#' @importFrom MASS mvrnorm
#' @importFrom mclust adjustedRandIndex
#' @importFrom stats median rnorm
NULL

utils::globalVariables(".data")

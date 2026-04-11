\name{dkge_between_permute}
\alias{dkge_between_permute}
\title{Permutation tests for between-subject DKGE RRR terms}
\usage{
dkge_between_permute(
  object,
  terms = NULL,
  method = c("freedman_lane"),
  B = 999L,
  blocks = NULL,
  seed = NULL,
  adjust = "none",
  statistic = c("frob"),
  scope = c("global", "features", "both"),
  feature_adjust = c("none", "fdr", "maxT")
)
}
\arguments{
\item{object}{A \code{dkge_between_rrr} object.}
\item{terms}{Character vector of terms or model-matrix columns to test.}
\item{method}{Permutation method. Currently \code{"freedman_lane"}.}
\item{B}{Number of permutations.}
\item{blocks}{Optional exchangeability blocks of length \code{n_subjects}.}
\item{seed}{Optional random seed.}
\item{adjust}{P-value adjustment method passed to \code{stats::p.adjust()}.}
\item{statistic}{Test statistic. Currently \code{"frob"}.}
\item{scope}{Whether to compute only global term tests, only featurewise term
tests, or both.}
\item{feature_adjust}{Multiplicity adjustment for featurewise tests: none,
Benjamini-Hochberg FDR, or maxT.}
}
\value{Object of class \code{dkge_between_permutation}.}
\description{
Performs term-specific global tests for a fitted \code{dkge_between_rrr()}
model. The default uses Freedman-Lane residual permutations from a reduced
model that excludes the tested term, then compares the reduced and full
reduced-rank residual sums of squares. When \code{scope} includes featurewise
inference, the object also stores observed term maps and per-feature
permutation p-values derived from the term-map norm.
}

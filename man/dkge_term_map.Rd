\name{dkge_term_map}
\alias{dkge_term_map}
\title{Extract a term-specific coefficient map}
\usage{
dkge_term_map(object, term, contrast = NULL, drop = TRUE)
}
\arguments{
\item{object}{A \code{dkge_between_rrr} object.}
\item{term}{Formula term or model-matrix column name.}
\item{contrast}{Optional numeric contrast over model-matrix columns.}
\item{drop}{Logical; drop single-row results to a vector.}
}
\value{Numeric vector or matrix of coefficients in target feature space.}
\description{
Extracts a coefficient map for a model term, model-matrix column, or explicit
contrast from a fitted between-subject reduced-rank model.
}

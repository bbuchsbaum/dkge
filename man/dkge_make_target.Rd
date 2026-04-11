\name{dkge_make_target}
\alias{dkge_make_target}
\title{Build a subject-by-feature target for between-subject DKGE models}
\usage{
dkge_make_target(
  fit = NULL,
  type = c("matrix", "component_scores", "transported_maps"),
  Y = NULL,
  contrast = NULL,
  contrast_obj = NULL,
  transport = NULL,
  values = NULL,
  loadings = NULL,
  centroids = NULL,
  sizes = NULL,
  medoid = NULL,
  mapper = NULL,
  crossfit = c("none", "analytic", "loso", "kfold"),
  feature_ids = NULL,
  subject_ids = NULL,
  coverage = NULL,
  precision = NULL,
  subject_weights = NULL,
  feature_weights = NULL,
  geometry = NULL,
  provenance = NULL,
  ...
)
}
\arguments{
\item{fit}{Optional \code{dkge} object used to derive DKGE targets.}
\item{type}{Target type.}
\item{Y}{Optional subject-by-feature matrix.}
\item{contrast, contrast_obj}{Contrast inputs for transported-map targets.}
\item{transport}{Optional \code{dkge_transport_spec}.}
\item{values}{Optional subject values as a matrix or list.}
\item{loadings}{Optional subject loading matrices used by transport.}
\item{centroids, sizes, medoid, mapper}{Transport inputs.}
\item{crossfit}{Cross-fitting method used when computing \code{contrast}.}
\item{feature_ids, subject_ids}{Optional identifiers.}
\item{coverage, precision}{Optional target reliability metadata.}
\item{subject_weights, feature_weights}{Optional reliability weights.}
\item{geometry, provenance}{Optional metadata.}
\item{...}{Additional arguments forwarded to contrast/transport helpers.}
}
\value{Object of class \code{dkge_target}.}
\description{
Creates the explicit boundary object between the DKGE representation stage and
between-subject multivariate models. The target stores a subject-by-feature
matrix plus optional coverage, precision, and weighting metadata.
}

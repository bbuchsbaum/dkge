\name{dkge_between_rrr}
\alias{dkge_between_rrr}
\alias{predict.dkge_between_rrr}
\alias{fitted.dkge_between_rrr}
\alias{residuals.dkge_between_rrr}
\title{Reduced-rank regression on a DKGE subject target}
\usage{
dkge_between_rrr(
  target,
  design,
  rank = NULL,
  weights = c("none", "target"),
  feature_mask = NULL,
  tol = 1e-08
)
}
\arguments{
\item{target}{A \code{dkge_target} from \code{dkge_make_target()} or a numeric matrix.}
\item{design}{A \code{dkge_subject_model} from \code{dkge_subject_model()} or a numeric model matrix.}
\item{rank}{Reduced rank. Defaults to the maximum estimable rank.}
\item{weights}{\code{"none"}, \code{"target"}, or a list with optional subject and feature weights.}
\item{feature_mask}{Optional logical vector selecting target features.}
\item{tol}{Numerical tolerance for rank checks.}
}
\value{Object of class \code{dkge_between_rrr}.}
\description{
Fits a directional multivariate second-level model
\deqn{Y = X B + E,\quad rank(B) <= r}
using a compact SVD of the fitted design-space response. This is intended for
subject-level factors, traits, and covariates after DKGE has produced a
subject-by-feature target.
}
\details{
\code{predict()} supports \code{type = "response"} for predicted subject-by-feature
values and \code{type = "scores"} for latent subject scores. When \code{newdata}
is a data frame and the fitted object retains a formula-based
\code{dkge_subject_model}, the model matrix is rebuilt automatically.
}

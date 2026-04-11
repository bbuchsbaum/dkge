\name{dkge_subject_model}
\alias{dkge_subject_model}
\title{Build a between-subject model matrix}
\usage{
dkge_subject_model(
  formula,
  data,
  subject_ids = NULL,
  nuisance = NULL,
  contrasts.arg = NULL,
  na.action = stats::na.fail
)
}
\arguments{
\item{formula}{A model formula such as \code{~ group * trait + age + sex}.}
\item{data}{Data frame containing subject-level variables.}
\item{subject_ids}{Optional subject identifiers.}
\item{nuisance}{Optional character vector of nuisance term labels.}
\item{contrasts.arg}{Optional contrasts passed to \code{stats::model.matrix()}.}
\item{na.action}{NA handler for the model frame.}
}
\value{Object of class \code{dkge_subject_model}.}
\description{
Wraps formula handling for subject-level DKGE analyses. The result owns the
model matrix, term metadata, and subject identifiers used by second-level
models such as \code{dkge_between_rrr()}.
}

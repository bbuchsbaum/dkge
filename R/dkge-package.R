#' dkge: Design-Kernel Group Embedding
#'
#' Core routines for design-aware group embedding, cross-validated contrasts, and transport utilities.
#'
#' @keywords internal
#' @importFrom stats aggregate contr.helmert contr.sum p.adjust pt setNames
#' @importFrom utils head
"_PACKAGE"

#' DKGE glossary
#'
#' Definitions for the terms the DKGE documentation uses most. Each one states
#' what the package does, not what the word means elsewhere.
#'
#' @section Effect space:
#' The shared coordinate system in which DKGE works. Each subject's beta matrix
#' has one row per named design effect, so the rows mean the same thing across
#' subjects even when the columns (clusters) do not. That row space is the
#' effect space, and its dimension `q` is the number of design effects. The
#' design kernel `K` supplies the metric on it.
#'
#' @section Design kernel:
#' A `q` by `q` positive semidefinite matrix stating which design effects should
#' be treated as similar. It expresses a scientific belief you can name, such as
#' adjacent ordinal levels being related. It is not a tuning matrix, and
#' `diag(q)` is the honest choice when the structure is uncertain. See
#' [design_kernel()].
#'
#' @section Salience:
#' The design-weighted component coordinates, \eqn{K U}. Rows are effects or
#' design cells and columns are latent components, so reading down a column
#' shows which effects a component combines. See [dkge_component_saliences()].
#'
#' @section Cross-fitting:
#' Evaluating a subject through a basis their own data did not shape. The
#' leave-one-subject-out path removes subject `s` from the pooled compressed
#' covariance, re-derives the basis by eigendecomposition of what remains, and
#' reads that subject's values through it. K-fold cross-fitting does the same
#' with folds instead of single subjects. It limits basis-reuse optimism; it
#' does not supply a population p-value.
#'
#' @section Medoid:
#' The reference subject whose clusters the other subjects are mapped onto.
#' Nothing computes it: it is an index you supply, defaulting to subject 1, and
#' the medoid subject's own values pass through the transport unchanged. A
#' reader who knows the clustering sense of the word should note that the
#' package does not select a most-representative subject.
#'
#' @section Transport:
#' Moving per-cluster values from each subject's own parcellation onto a common
#' set of locations, so that fields defined on different clusterings can be
#' compared or averaged. The plan is built from cluster centroids by one of
#' several mappers (`"sinkhorn"`, `"ridge"`, `"ols"`, or a k-nearest-neighbour
#' variant). Transport is a spatial alignment step and it is separate from the
#' effect-space alignment that `K` performs.
#'
#' @section LOSO:
#' Leave-one-subject-out. See Cross-fitting.
#'
#' @name dkge-glossary
#' @keywords documentation
NULL

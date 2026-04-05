# dkge-classify.R
# This file has been split into:
#   dkge-classify-core.R    — public API, target prep, S3 methods
#   dkge-classify-cv.R      — cross-validation loops (cell and delta targets)
#   dkge-classify-backends.R — classifier backends (LDA, logistic regression)

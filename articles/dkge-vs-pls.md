# DKGE versus Partial Least Squares

This document provides a comparative analysis of Design-Kernel Group
Embedding (DKGE) and the well-established partial least squares (PLS)
framework for neuroimaging analysis. To establish a foundation for
comparison, we first examine the classical PLS approach as described in
McIntosh & Lobaugh (2004), following their established terminology and
conceptual framework.

------------------------------------------------------------------------

## What classical PLS does

Classical PLS neuroimaging analysis operates by constructing latent
variables that capture the strongest relationships between brain
activity and experimental or behavioral variables. The **primary
objective** is to build latent variables (LVs) that maximise the
covariance between an exogenous block (such as task design or behavior)
and a brain-data block (elements × time).

The framework encompasses several **methodological variants**, each
tailored to different research questions: task PLS analyzes condition
differences, behavior PLS examines brain–behavior coupling, seed PLS
investigates functional connectivity patterns, and spatiotemporal PLS
treats space × time jointly for fMRI/ERP/MEG data.

The **computational workflow** follows a systematic sequence of
steps: 1. Data are arranged as a single matrix
$M \in {\mathbb{R}}^{{(nk)} \times {(mt)}}$ with observations nested
inside conditions. 2. Cross-block covariance is formed with an
orthonormal design matrix, followed by singular value decomposition
$C^{\top}M = USV^{\top}$. 3. Element/time saliences (singular images),
design saliences, and singular values are extracted from the
decomposition. 4. Brain scores $B = MU$ and design scores $D = CV$ are
computed to characterize patterns in each domain. 5. For behavior PLS
specifically, the design block is replaced with behavior matrices and
correlated with $M$ before applying SVD.

**Statistical inference** in PLS relies on resampling approaches:
permutation tests assess LV significance, bootstrap procedures evaluate
voxel salience reliability, and Procrustes alignment stabilizes
resampled LVs across iterations.

This methodology emphasizes *whole-pattern* effects and time-resolved
couplings, leveraging linear algebra and resampling techniques to
identify distributed, reliable neural patterns that relate to
experimental or behavioral variables.

## How DKGE relates

DKGE builds upon the fundamental concept of using a compact latent space
to explain patterns in neuroimaging data, but extends this approach in
ways specifically tailored for modern multi-subject neuroimaging
analysis. While PLS operates on raw brain data to find latent variables,
DKGE works with subject-wise GLM outputs and incorporates design-aware
alignment and spatial transport mechanisms. The comprehensive comparison
table below summarizes both the conceptual commonalities between these
approaches and the key methodological advances that DKGE introduces.

| Aspect                    | Partial Least Squares                                                                                        | DKGE                                                                                                                                  |
|---------------------------|--------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| Latent-space construction | SVD on cross-block covariance between design/behaviour and brain data; columns of $U$ and $V$ are saliences. | Eigen-decomposition of a design-kernel-weighted covariance, producing orthonormal components $U$.                                     |
| Design information        | Implicit via orthonormal contrasts; no way to encode graded similarity between effects.                      | Explicit design kernel $K$ encodes factorial structure, smoothness, or prior correlations among effects.                              |
| Data normalisation        | Conditions averaged or mean-centred before SVD; each voxel treated equally.                                  | Row standardisation of subject betas using the pooled design Cholesky factor; optional spatial/reliability weights $\Omega_{s}$.      |
| Cross-validation          | Permutation for LV significance, bootstrap for salience stability (no LOSO cross-fitting).                   | LOSO / K-fold cross-fitting (`dkge_contrast`), analytic approximations, parametric or bootstrap inference with cached transports.     |
| Transport / alignment     | Outputs latent scores; spatial interpretation relies on the original voxel grid.                             | Provides barycentric kNN and Sinkhorn transports, anchor graphs, and voxel decoders for consistent spatial maps across parcellations. |
| Reliability weighting     | All observations weighted equally; stability assessed post hoc via bootstrap ratios.                         | Subject- and cluster-level reliabilities enter directly (e.g. sizes, inverse variances), influencing fits and transport.              |
| Spatiotemporal support    | ST-PLS handles time by stacking features.                                                                    | DKGE works on any GLM-derived beta blocks; temporal modelling is delegated to the design matrix and optional kernels.                 |
| Implementation focus      | Exploratory LVs; complementary to other analyses.                                                            | Integrated workflow for group GLM analysis, transport, and inference tailored to fMRI/ERP pipelines.                                  |

## Similarities worth noting

Despite their methodological differences, PLS and DKGE share several
fundamental characteristics that reflect their common mathematical
foundations. Both approaches rely on SVD or eigendecomposition
techniques to obtain orthogonal latent patterns and their associated
scores, ensuring that the derived components capture independent sources
of variation in the data.

Resampling methodology plays a central role in both frameworks, though
implemented differently: PLS employs permutation tests and bootstrap
procedures for statistical inference, while DKGE provides a broader
toolkit including analytic approximations, leave-one-subject-out (LOSO)
cross-fitting, bootstrap procedures, and transport-aware resampling
utilities.

Both methods require careful interpretation that involves examining
latent loadings or saliences in conjunction with subject scores to
properly understand how experimental conditions or behavioral variables
relate to the underlying neural patterns.

## Why DKGE can replace or complement PLS in modern pipelines

DKGE offers several methodological advantages that make it particularly
well-suited for contemporary neuroimaging analysis pipelines, addressing
limitations that can arise with classical PLS approaches:

1.  **Enhanced design control** — The design kernel framework provides a
    principled mechanism for encoding factorial relations, smoothness
    constraints, or hierarchical structure among experimental
    conditions, capabilities that are difficult to achieve within the
    classical PLS framework.

2.  **Bias-aware statistical contrasts** — Cross-fitting procedures
    systematically avoid the optimistic bias that can occur when
    estimating condition or behavior effects on the same data used for
    model fitting, while built-in inference tools provide rigorous
    statistical testing without requiring separate implementations.

3.  **Sophisticated spatial alignment** — Transport operators enable
    mapping of component or contrast values onto consistent anchor or
    voxel grids, facilitating meaningful group-level interpretation even
    when working with heterogeneous parcellation schemes across
    subjects.

4.  **Principled reliability weighting** — Subject-level and
    cluster-level reliability weights ensure that noisy measurements
    contribute proportionally less to the final estimates, improving
    overall stability compared to the uniform weighting typically
    employed in PLS analyses.

5.  **Seamless workflow integration** — DKGE integrates directly into
    downstream analysis pipelines including rendering, bootstrapping,
    and component analysis without requiring reconstruction of transport
    operators or design alignment logic at each stage.

In summary, DKGE can be conceptualized as a design-aware,
transport-enabled extension of fundamental PLS concepts. It preserves
the interpretability and whole-pattern emphasis that make latent
variable approaches valuable while providing the inference capabilities,
reliability weighting, and spatial alignment mechanisms that modern
multi-subject neuroimaging studies require.

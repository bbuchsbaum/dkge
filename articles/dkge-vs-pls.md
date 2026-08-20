# DKGE versus Partial Least Squares

This page is a conceptual orientation, not an accuracy or power
benchmark. “PLS” names several related methods, and implementation
choices matter. The summary follows the neuroimaging accounts of
[McIntosh and Lobaugh
(2004)](https://doi.org/10.1016/j.neuroimage.2004.07.020) and [Krishnan
et al. (2011)](https://doi.org/10.1016/j.neuroimage.2010.07.034). It
compares their classical task/behavior PLS formulations with DKGE’s
current subject-beta workflow so you can decide which estimand and input
contract match your study.

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
$`M \in \mathbb{R}^{(n k) \times (m t)}`$ with observations nested
inside conditions. 2. Cross-block covariance is formed with an
orthonormal design matrix, followed by singular value decomposition
$`C^\top M = U S V^\top`$. 3. Element/time saliences (singular images),
design saliences, and singular values are extracted from the
decomposition. 4. Brain scores $`B = M U`$ and design scores $`D = C V`$
are computed to characterize patterns in each domain. 5. For behavior
PLS specifically, the design block is replaced with behavior matrices
and correlated with $`M`$ before applying SVD.

**Statistical inference** in PLS relies on resampling approaches:
permutation tests assess LV significance, bootstrap procedures evaluate
voxel salience reliability, and Procrustes alignment stabilizes
resampled LVs across iterations.

This methodology emphasizes *whole-pattern* effects and time-resolved
couplings, leveraging linear algebra and resampling techniques to
identify distributed, reliable neural patterns that relate to
experimental or behavioral variables.

## How DKGE relates

Both methods construct compact latent representations, but they begin
from different objects. Classical task or behavior PLS constructs a
cross-block matrix linking a design/behavior block to brain
measurements. DKGE pools subject-level GLM effect moments and applies an
explicit effect-space kernel. The table records implemented differences;
it does not establish that one method is generally more accurate,
stable, or powerful.

| Aspect | Partial Least Squares | DKGE |
|----|----|----|
| Latent-space construction | SVD on cross-block covariance between design/behaviour and brain data; columns of $`U`$ and $`V`$ are saliences. | Eigen-decomposition of a design-kernel-weighted covariance, producing orthonormal components $`U`$. |
| Design information | Specified through the exogenous/design block; the exact coding and available structure depend on the PLS variant. | Explicit design kernel $`K`$ encodes factorial structure, smoothness, or prior correlations among effects. |
| Data normalisation | Conditions averaged or mean-centred before SVD; each voxel treated equally. | Row standardisation of subject betas using the pooled design Cholesky factor; optional spatial/reliability weights $`\Omega_s`$. |
| Cross-validation | Classical neuroimaging workflows commonly use permutation for LV significance and bootstrap for salience stability; predictive cross-validation can be added but is not the same estimand. | LOSO / K-fold cross-fitting (`dkge_contrast`), analytic approximations, parametric or bootstrap inference with cached transports. |
| Transport / alignment | Outputs latent scores; spatial interpretation relies on the original voxel grid. | Provides barycentric kNN and Sinkhorn transports, anchor graphs, and voxel decoders for consistent spatial maps across parcellations. |
| Reliability weighting | Baseline formulations often use uniform observation weights; weighted, sparse, and regularized PLS variants also exist. | Subject- and cluster-level reliabilities enter directly (e.g. sizes, inverse variances), influencing fits and transport. |
| Spatiotemporal support | ST-PLS handles time by stacking features. | DKGE works on any GLM-derived beta blocks; temporal modelling is delegated to the design matrix and optional kernels. |
| Implementation focus | Exploratory LVs; complementary to other analyses. | Integrated workflow for group GLM analysis, transport, and inference tailored to fMRI/ERP pipelines. |

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

## Which method matches the question?

Choose a PLS formulation when the target quantity is cross-block
covariance between a brain block and an experimental or behavioral block
and the chosen PLS implementation supplies the resampling procedure you
need. Choose DKGE when the inputs are already subject-level GLM beta
matrices, the scientific model requires an explicit effect-space metric,
or subject-specific parcellations must be transported before
feature-level summaries.

These are workflow distinctions, not demonstrated performance
advantages. DKGE’s cross-fitting addresses basis reuse only when the
selected mode actually refits the basis; reliability weights help only
when their construction is defensible; and spatial transport is useful
only when its features and mapper support plausible correspondence.
Conversely, the table does not cover the many sparse, regularized,
predictive, or multilevel PLS variants.

For an empirical comparison, prespecify one estimand, fit both methods
to the same analysis units, use identical outer resampling splits, and
compare a held-out metric with uncertainty. Do not infer superiority
from the fact that one package exposes a helper that the other workflow
leaves to the analyst.

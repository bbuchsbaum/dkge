Below is a concise white paper for engineers implementing Design‑Kernel Group Embedding (DKGE). It specifies data shapes, the mathematics, core algorithms, complexity, and practical engineering guidance.

⸻

DKGE: Design‑Kernel Group Embedding for Cluster‑level fMRI

Goal. Learn a shared, design‑respecting latent basis from per‑subject GLM coefficients, then compute unbiased contrasts and display group maps without voxel‑wise averaging.
	•	Compression: all heavy LA is in q×q, where q = #design effects (often 10–100).
	•	Robustness: LOSO (leave‑one‑subject‑out) cross‑fitting prevents circularity.
	•	Display: transport cluster values to a medoid parcellation (mass‑preserving OT).
	•	Design structure: encoded by a kernel K over effects (main, interactions, ordinal smoothness).

⸻

1) Data model & notation

Symbol	Shape	Meaning
Y_s	T_s\times P_s	Subject‑s cluster time‑series (columns are clusters)
X_s	T_s\times q	Subject‑s design matrix (effects / regressors)
\widehat B_s	q\times P_s	GLM betas per cluster: (X_s^\top X_s)^{-1} X_s^\top Y_s
\Omega_s	P_s\times P_s (diag)	Optional cluster reliability/size weights
K	q\times q	Design kernel (PSD) in effect space
R	q\times q	Shared “ruler”, s.t. R^\top R = \sum_s X_s^\top X_s
\widetilde B_s	q\times P_s	Row‑standardized betas: \widetilde B_s=R^\top\widehat B_s
K^{\pm 1/2}	q\times q	Kernel (inverse) square root
\widehat C	q\times q	Compressed covariance in K^{1/2} metric
U	q\times r	Group latent basis, U^\top K U = I_r
A_s	P_s\times r	Subject‑s cluster loadings: A_s=\widetilde B_s^\top K U
c	q\times 1	Contrast in original design basis
\alpha	r\times 1	Contrast coords: \alpha=U^\top K R^{-1}c
v_s	P_s\times 1	Cluster‑wise contrast values for subject s

Design kernel K: encode similarity between design effects by building per-factor kernels (nominal = identity, ordinal/circular/continuous = RBF with tunable length-scales), combine terms via Kronecker products with weights, add a small ridge, then optionally map to effect space with contrasts T (K = T^\top K_{\text{cell}} T).
Use `design_kernel()` to generate K in either cell or effect basis (with metadata on contrast blocks).

⸻

2) Core mathematics
	1.	Shared ruler (coding/scale invariance).
G_{\text{pool}}=\sum_s X_s^\top X_s,\; R^\top R=G_{\text{pool}}, then \widetilde B_s=R^\top \widehat B_s.
	2.	Compressed covariance in the design metric.
\widehat C \;=\; \sum_s w_s\;K^{1/2}\,\widetilde B_s\,\Omega_s\,\widetilde B_s^\top\,K^{1/2}.
Subject weights w_s (optional, MFA‑style) stabilize across subjects.
	3.	Generalized eigen in K-metric.
Eigen‑decompose \widehat C = V\Lambda V^\top; set U = K^{-1/2}V_{[:,1:r]} → U^\top K U=I_r.
	4.	Projections & contrasts.
A_s = \widetilde B_s^\top K U; for any contrast c, compute \alpha = U^\top K R^{-1}c, then v_s=A_s\alpha.
	5.	Unbiased estimation.
For subject s, re‑learn U^{(-s)} using \widehat C^{(-s)} (remove s’s contribution) and evaluate v_s with U^{(-s)}.

⸻

3) Algorithm 1 — DKGE fit (tiny q\times q)

Inputs: {B_s ∈ ℝ^{q×P_s}, X_s ∈ ℝ^{T_s×q}, Ω_s diag (optional)}_{s=1..S}, K ∈ ℝ^{q×q}, rank r, ridge ≥ 0
Output: U ∈ ℝ^{q×r} (K-orthonormal), {B̃_s}, K^{±1/2}, Chat, {per-subject contributions}

1. G_pool ← Σ_s X_sᵀ X_s;  R ← chol(G_pool)  # upper: Rᵀ R = G_pool
2. For each s:  B̃_s ← Rᵀ B_s
3. Eigendecompose K → K = V diag(κ) Vᵀ, κ_i ≥ ε  → K^{±1/2} = V diag(κ_i^{±1/2}) Vᵀ
4. Weights (optional):
     w_s ← 1 / σ₁(K^{1/2} B̃_s Ω_s^{1/2})²         # shrink toward 1
5. Accumulate tiny matrix:
     Chat ← 0
     For each s:
        right_s ← B̃_s Ω_s B̃_sᵀ   (column-scaling if Ω diagonal)
        S_s ← K^{1/2} right_s K^{1/2}
        Chat ← Chat + w_s S_s
6. If ridge > 0: Chat ← Chat + ridge·I
7. Eig(Chat) → V Λ Vᵀ;  U ← K^{-1/2} V[:,1:r]
Return U, {B̃_s}, K^{±1/2}, Chat, {S_s}, R, weights

Complexity. Per subject: O(qP_s) for B̃_s Ω_s B̃_s^\top (Ω diagonal), plus O(q^2) update; eigen O(q^3).
Memory. O(q^2) (store K, K^{±1/2}, \widehat C); cluster data never forms a large square.

⸻

4) Algorithm 2 — LOSO cross‑fitted contrasts

Inputs: fit (U, R, K, K^{±1/2}, Chat, {S_s}, weights, {B̃_s}), c ∈ ℝ^q
Output: {v_s}  (unbiased per-subject cluster contrasts)

Precompute: c̃ = R^{-1} c
For each subject s:
  Chat^{(-s)} ← Chat − w_s S_s
  Eig(Chat^{(-s)}) → V Λ Vᵀ
  U^{(-s)} ← K^{-1/2} V[:,1:r]
  α ← U^{(-s)ᵀ} K c̃
  A_s ← B̃_sᵀ K U^{(-s)}
  v_s ← A_s α

Note. No time split is required; LOSO uses subject‑held‑out geometry and avoids circularity.

⸻

5) Algorithm 3 — Transport to medoid (entropic OT)

Cost. For subject s cluster i and medoid cluster j:
C_{ij}=\lambda_{\text{emb}}\|A_{s,i\cdot}-A_{m,j\cdot}\|^2 + \lambda_{\text{spa}}\|x_{s,i}-x_{m,j}\|^2/\sigma^2.

Inputs: {v_s ∈ ℝ^{P_s}}, {A_s ∈ ℝ^{P_s×r}}, centroids {X_s ∈ ℝ^{P_s×3}},
        medoid index m, masses {μ_s}, ν_m, cost weights, ε>0
Output: Group values on medoid clusters (length Q)

For each subject s:
  If s = m: y_s ← v_s; continue
  Build cost C_s (embedding + spatial)
  Solve entropic OT (Sinkhorn) with ε: plan T_s ≈ argmin ⟨C_s, T⟩ + ε·KL(T||μ_s⊗ν_m)
  Pushforward: (y_s)_j ← (Σ_i v_{s,i} T_{s,ij}) / (Σ_i T_{s,ij})
Aggregate: group_j ← median_s (y_s)_j    # or robust mean

Complexity. Sinkhorn ~ O(P_s Q) per iteration; use log‑domain updates to avoid underflow.
Note. When memory is tight, compute pushforward without storing T_s (accumulate numerators/denominators).

⸻

6) Optional: CPCA “inside‑span” filter (q‑space)

To isolate design‑explainable structure before eigen:
	•	With row subspace T\in\mathbb{R}^{q\times q_0}:
\widehat P = K^{1/2} T (T^\top K T)^{-1} T^\top K^{1/2}.
Split \widehat C into \widehat C_{\text{design}}=\widehat P \widehat C \widehat P and
\widehat C_{\text{resid}}=(I-\widehat P)\widehat C (I-\widehat P); eigendecompose each → U_{\text{design}}, U_{\text{resid}}.

This keeps parts inside span and orthogonal in the stacked data, improving interpretability. (Toggleable; zero change in memory profile.)

⸻

7) Out‑of‑sample prediction (frozen basis)

Given a frozen model \{U, K, R\}:
	•	Loadings: A_s = \widetilde B_s^\top K U, \widetilde B_s=R^\top B_s.
	•	Contrasts: \alpha = U^\top K R^{-1}c, v_s = A_s \alpha.
	•	Streaming: apply per subject; no retraining needed.

⸻

8) Engineering guidance

Shapes & typical sizes. q ~ 10–100; P_s ~ 1–10k clusters; S ~ 10–200.
Precision. Use double precision; add small ridge in K eigenspectrum and in \widehat C if eigen‑gaps are small.
Symmetry. Symmetrize tiny matrices before eigen: A\leftarrow(A+A^\top)/2.
Parallelism. Accumulate \widehat C across subjects (OpenMP; per‑thread accumulators) → single reduction.
Streaming. Two‑pass design:
	1.	accumulate G_{\text{pool}} → R;
	2.	stream subjects to update \widehat C and optionally cache S_s on disk for fast LOSO.
Memory. Store only q\times q matrices; never build large P_s\times P_s objects.
GLM adapters. Compute \widehat B_s from any GLM engine (OLS/GLS); the DKGE core expects only \widehat B_s and X_s.
Design kernel. Modular builder: main effects, interactions, ordinal/circular smoothness; validate PSD; normalize (unit trace/Frobenius).

⸻

9) Model selection & inference
	•	Rank r. LOSO explained variance (in K^{1/2} metric) with one‑SE rule.
	•	Kernel K tuning. Grid over term weights / length‑scales; pick by LOSO EV or held‑out contrast accuracy.
	•	Group inference.
	•	Fast: sign‑flip max‑T on transported subject maps (FWER across clusters).
	•	Heavy: Freedman–Lane (permute residuals at time‑series level; refit GLMs; re‑run DKGE) to reflect design‑fitting uncertainty.

⸻

10) Algorithm 4 — Streaming DKGE fit (two‑pass)

Inputs: loader with n(), X(s), B(s), Ω(s); K; rank r; ridge; cache_dir (optional)
Output: dkge_stream object: U, Chat, R, K^{±1/2}, weights, {S_s or paths}

Pass 1: G_pool ← Σ_s X(s)ᵀ X(s);  R ← chol(G_pool)
Compute K^{±1/2}

Pass 2: Initialize Chat ← 0
For each subject s:
  B̃_s ← Rᵀ B(s)
  w_s  ← subject weight (optional MFA)
  S_s  ← K^{1/2} (B̃_s Ω(s) B̃_sᵀ) K^{1/2}
  Chat ← Chat + w_s S_s
  If cache_dir: save S_s to disk for LOSO
If ridge>0: Chat += ridge·I
Eig(Chat) → V;  U ← K^{-1/2} V[:,1:r]

LOSO (streamed). For subject s, read S_s from cache, form \widehat C^{(-s)}, eig, then project B_s on the fly to get A_s, v_s.

⸻

11) Minimal API to implement
	•	Kernel builder: design_kernel(factors, terms, rho, basis="effect", contrasts=...) → K (+ K^{±1/2}).
	•	Fit: dkge_fit(B_list, X_list, K, Omega_list, ridge, rank, w_method, w_tau) (+ streamed variant).
	•	LOSO: dkge_loso_contrast(fit, s, c) and dkge_all_loso_contrasts(fit, c).
	•	Transport: dkge_transport_to_medoid_sinkhorn(...) (and a CPU‑fast log‑domain variant).
	•	Inference: dkge_signflip_maxT(Y, B) (and an optional FL scaffold).
	•	Predict: dkge_freeze(fit), dkge_predict(model, B_new_list, contrasts).
	•	Utilities: Procrustes alignment, consensus basis, CPCA inside‑span split.

⸻

12) Sanity checks & tests
	•	Unit tests
	•	U^\top K U = I to tolerance.
	•	LOSO vs full: basis changes but unbiased v_s behave as expected.
	•	Transport: OT marginals correct; pushforward equals explicit plan.
	•	Streaming vs batch: same U, \widehat C within FP tolerance.
	•	Inference: sign‑flip controls FWER on null simulations.
	•	Diagnostics
	•	Spectrum of \widehat C; eigen‑gaps.
	•	Weight distribution w_s.
	•	Alignment of bases across folds (Procrustes singular values).
	•	Coverage/dispersion per medoid cluster.

⸻

Appendix: Typical numeric safeguards
	•	Add \epsilon I to K eigenvalues before forming K^{±1/2}.
	•	Symmetrize K, \widehat C, and \widehat C^{(-s)}.
	•	Use backsolve with Cholesky R for R^{-1}c.
	•	Small ridge in \widehat C when eigen‑gaps are small.
	•	Normalize kernels (unit trace) to stabilize tuning.

⸻

Deliverables
Implement the algorithms above with:
	•	C++ kernels (Rcpp/Armadillo): accumulation of \widehat C, LOSO basis, log‑domain Sinkhorn, pairwise distances.
	•	Parallelism: subject‑level OpenMP; optional GPU for OT if desired.
	•	Streaming adapters: neuroim2 for on‑disk series; fmrireg for GLMs.

This document is sufficient to implement DKGE end‑to‑end in any scientific stack while preserving its small‑matrix, design‑faithful, and honest analysis guarantees.
# Predeclared q = 60 performance gate

Status: frozen before benchmark execution for Mote ticket
`bd-01M03KTXQK53KJD0HF1K8Y39CM`.

## Command and host contract

- Canonical command: `Rscript dev/benchmark-dkge-q60.R`.
- Three fresh-process repetitions are run per workload with seeds 19421:19423.
- Results are written to `inst/extdata/dkge-q60-performance.csv`.
- The benchmark records elapsed time, peak resident set size reported by
  `/usr/bin/time -l`, total allocations and largest allocation from `Rprofmem`,
  returned-object size, and retained design size separately.
- Thresholds below apply to this macOS reference host and are release smoke
  gates, not universal performance guarantees.

## Workloads

1. `trial_chunks`: `T = 240`, `q = 60`, `P = 100,000`, feature chunks of
   5,000. The full trial response would occupy 192 MB but is generated one
   block at a time. The returned beta matrix is 48 MB.
2. `group_legacy`: 18 subjects, `q = 60`, `P = 10,000` per subject, strongly
   unequal count precision and 15% structural missingness, ordinary pooled fit.
3. `group_advanced`: the same 18-subject workload with explicit precision,
   pair normalization, and analytic debiasing.
4. `fold_refit`: one leakage-free training refit from the advanced workload.
5. `poisson_bootstrap`: 20 literal-multiplicity q-space bootstrap repools on
   18 subjects with `P = 2,000` per subject.

The group workloads use precomputed subject-level sufficient statistics; they
do not retain trialwise `Y`. Trial reduction retains the shared `T x q` design
because it is part of the public subject record. Its size is reported.

## Complexity contract

| Stage | Dominant time | Retained memory | Forbidden allocation |
|---|---|---|---|
| chunked IID trial reduction | O(T q P + q^2 P) | O(q P + T q + P) | P x P |
| raw/analytic moment | O(q P + q^2 P) | O(q P + q^2 + P) | P x P |
| precision pair pooling | O(S q^2) after moments | O(S q^2) | P x P |
| pooled eigensolve | O(q^3) | O(q^2) | P x P |
| fold refit | same as one training fit | subject records plus O(q^2) | P x P |
| literal Poisson repool | O(B S q^2 P) currently | input records plus O(q^2) | P x P |

The last row is intentionally conservative: literal refits guarantee correct
random-effects and ruler semantics but are not advertised as a constant-time
bootstrap shortcut.

## Predeclared thresholds

Using the median of three runs:

- `trial_chunks`: elapsed <= 20 s, peak RSS <= 1.25 GB, returned object <=
  100 MB, retained design <= 0.25 MB, largest allocation <= 100 MB.
- `group_legacy`: elapsed <= 20 s and peak RSS <= 2.0 GB.
- `group_advanced`: elapsed <= 30 s and peak RSS <= 2.0 GB.
- `fold_refit`: elapsed <= 30 s and peak RSS <= 2.0 GB.
- `poisson_bootstrap`: elapsed <= 60 s and peak RSS <= 1.5 GB.
- Across repeated runs, elapsed IQR must not exceed the median for any workload.
- Source and allocation audits must find no `P x P` object. For the
  `trial_chunks` imaging workload, the largest allocation must remain below the
  80 GB required by a 100,000 x 100,000 double matrix by at least three orders
  of magnitude.

Any failed threshold is retained and reported without post-hoc revision.

# Functional-alignment v9 certification

This directory freezes certification for one source-tree-addressed candidate.
It does not replace the historical court artifacts. The byte-exact v2 protocol,
its fail-closed invalidation receipt, and the extracted Sinkhorn regression
fixture are retained beside the active protocol. The corresponding v2 source
test manifest, results, session receipt, and failed-court budget remain hashed
inputs under `data-raw/functional-alignment-certification/`. V2 materialized no
statistical outputs: a required solve ended just above the fixed marginal
tolerance at its 5,000-iteration ceiling. V3 changed only that ceiling to
50,000; correcting the output-directory label from v2 to v3 was explicitly
metadata-only. Before any v3 executor ran, a concurrent spatial-regularization
edit changed the package candidate. The byte-frozen `protocol-v3.json` and its
supersession receipt are therefore retained as immutable history. V4 bound a
new candidate and completed the formal computation, but its source tree changed
during the run. The analyzer detected that mismatch before reading either raw
or summarized statistical output. Its terminal invalidation receipt binds every
materialized V4 byte as inadmissible history. V5 was superseded before freeze
approval. V6's provisional freeze was subsequently revoked. Neither protocol
was executed after later candidate hardening changed the source tree. Their
frozen protocols and supersession receipts are also retained. V7 completed its
source suite, court, and known-warp study, but the exact built-tarball check then
exposed documentation and test-harness defects before collection. Those V7
results and the failed artifact candidates remain immutable, inspected history;
they were not reused as V8 evidence. V8 then completed a fresh source suite,
court, known-warp study, closed-world tarball build, and exact package check.
The separately required pkgdown render failed before collection because all 19
vignette version guards called the non-exported `utils::package_version()`.
Every V8 result is retained as historical evidence and is not reused as V9
evidence. V9 changes only those 19 calls to `base::package_version()` and adds a
validation-only guard regression; its certified runtime source hash is unchanged.
V4 through V9 make no numerical or statistical change. Every epsilon, tolerance,
seed, cell, arm, estimand, sign action, and gate is inherited, and a second
numerical amendment is forbidden. The active protocol's
`candidate_source_tree_sha256` binds the runtime candidate exactly, while the
source-suite, tarball, and documentation receipts bind the full package bytes.

Every executable uses two-sided loaded-byte binding. It snapshots the package
source and its complete test, court, or power dependency bundle; reloads every
helper and, where applicable, the package with compilation enabled; then
requires the same snapshot before computation begins. The same snapshot is
checked again before a successful receipt is written. A concurrent edit thus
creates failed, inadmissible evidence rather than allowing loaded code and
recorded bytes to have mixed provenance.

For V9, run the complete source suite and exact tarball/check/pkgdown gate first.
Each subsequent court or power runner writes into a directory whose name contains the live package
source-tree hash and the exact validation-harness bundle hash. It refuses to
reuse any already claimed run directory. Before computation, the court writes
its immutable harness and a `running` receipt; it atomically records either
`failed` or `formal_complete`, and analysis alone can transition that receipt
to `analyzed_complete`. The power runner likewise records `running` before any
cohort, then `failed` or `power_complete`; determinism verification first checks
the recorded output hashes and alone transitions to `determinism_complete`.
Thus a numerical failure or partial evidence directory cannot be silently
retried. Both runners record fixed seeds and configuration and retain adverse outcomes.
Power is descriptive: it cannot promote an alignment mode or override Type-I,
provenance, or eligibility gates. The source-suite runner
additionally addresses the exact tests and validation scripts, so test-only
repairs cannot silently reuse older evidence.

Run the complete source suite only after package code and validation helpers
are frozen:

```sh
Rscript inst/validation/functional-alignment-certification/run-source-tests.R
```

Build the closed-world candidate and check that exact tarball, retaining the
checker's `00_pkg_src` directory:

```sh
R CMD build --no-build-vignettes /path/to/dkge
R CMD check --no-manual --no-clean --ignore-vignettes \
  /path/to/dkge_version.tar.gz
```

The collector requires `--no-manual` and `--ignore-vignettes` in the check log
and rejects every option outside the declared allowlist, including any option
that weakens installation, tests, examples, codoc, or subdirectory checks.
Vignettes are excluded from the check only because the exact-tarball pkgdown
render is the separately collected executable-documentation evidence. Build
pkgdown from the same closed-world tarball through the
receipt-producing wrapper:

```sh
LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 \
Rscript inst/validation/functional-alignment-certification/run-pkgdown.R \
  --tarball=/path/to/dkge_version.tar.gz \
  --dest-dir=/path/to/pkgdown-site
```

The wrapper fails before rendering unless R reports a UTF-8 locale. Its receipt
records both that gate and the full locale string, and the collector rejects a
missing or false locale claim. This prevents a hash-valid receipt from
certifying encoding-corrupted article or figure text.

Only after that receipt validates, run the fresh held-out efficacy study and
formal court from the repository root:

```sh
Rscript inst/validation/functional-alignment-certification/run-known-warp-power.R
Rscript inst/validation/functional-alignment-certification/run-court.R
```

After the full source suite, built-tarball check, pkgdown build, and independent
review of the fresh court and power evidence are complete, assemble the
exact-candidate record with:

```sh
Rscript inst/validation/functional-alignment-certification/collect-certification.R \
  --tarball=/path/to/dkge_version.tar.gz \
  --check-dir=/path/to/dkge.Rcheck \
  --pkgdown-dir=/path/to/pkgdown-site \
  --review=data-raw/functional-alignment-certification/independent-review.json
```

`R CMD build` rewrites `DESCRIPTION` and normally omits `_pkgdown.yml` plus the
`pkgdown/` assets. The wrapper therefore verifies an explicit build projection:
all retained DESCRIPTION fields must be semantically identical, every packaged
source and documentation file must be byte-identical, and the omitted site
configuration is applied as a separately hashed overlay. Both the raw tarball
hash and the certified projection/overlay hashes are recorded in the receipt.
New certification stages also radix-sort source paths so the fingerprint is
independent of the process collation locale while reproducing the frozen
court's canonical source hash.

The collector unpacks the tarball and independently verifies that projection,
the exact checked source tree, the pkgdown site receipt, every court/power
helper, the complete V2-V7 invalidation and supersession chain, and the
source-test tree. The independent review must name these same hashes, the
analyzed court run-state, the aggregate evidence-output bundle, and the
collector hash, and it must contain no open blockers. The collector rebuilds
the exact frozen court and known-warp schedules from the protocol and recomputes
both summaries from raw rows. Collection fails if any binding differs or if the
frozen court is represented as permitting inferential promotion. Final
publication is staged and renamed atomically only after every gate passes.

# DKGE public-beta baseline freeze

Frozen on 2026-08-20 for P0 epic `bd-01M0GADWESTKW3J0AGVN6YKQR1`.
This record separates the audited remote source, the user-owned dirty overlay,
and the implementation candidate. None of them is evidence of a shipped beta.

## Source identities

| surface | identity | state at freeze |
|---|---|---|
| local checkout `main` | `30f5b5db2a9ff4d5422f2fb0db0f55fa88141826` | 0 ahead, 2 behind its tracking ref; broad dirty overlay |
| local tracking ref `origin/main` | `1fdacb37ef2320a942e08f09eb5078033169478c` | matched the live remote query |
| live `refs/heads/main` | `1fdacb37ef2320a942e08f09eb5078033169478c` | verified with `git ls-remote origin refs/heads/main` |
| isolated implementation branch | `codex/public-beta-p0` from `1fdacb37ef2320a942e08f09eb5078033169478c` | worktree `/private/tmp/dkge-beta-p0`; not published |
| dirty-overlay patch | SHA-256 `fe1ebd7e0f1e861a2c7e636e09ec14a90025652ff9b264f8260a71b5c706d4f6` | 141,156 bytes; copied only into the isolated worktree |

`data-raw/dkge-public-beta-dirty-manifest.csv` records every modified or
untracked path in the original checkout, including its status, byte size, file
SHA-256, and source HEAD. The original checkout is not edited by this epic.

## Exact reproduction sources

The R-level reproductions in `artifacts/review-verification/results.txt` were
run by `reproduce-review.R` against an archived copy of the live remote R
source. The frozen SHA is `1fdacb37ef2320a942e08f09eb5078033169478c`.
The harness loaded the checkout only to reuse compiled native code; `src/` had
no difference between local HEAD, dirty overlay, and that remote SHA. The
native warm-start reproduction used the package built from the same archived
remote source. The focused strength tests and both `R CMD check` logs also used
that archive. The regression files named in
`data-raw/dkge-public-beta-traceability.csv` are run against a pristine archive
of this SHA before any P0 implementation edits.

That pre-fix run completed with 32 failures, 1 warning, 0 skips, and 10 oracle
passes; `data-raw/dkge-public-beta-baseline-failures.md` records the scope. The
release-contract fixture separately freezes the dependency and workflow gaps.

## Reconciliation with earlier work

The dirty overlay contains the completed August 17/20 new-path remediation,
documentation, generated `.Rd` updates, and associated tests. It is retained
as user-owned work, not treated as released code. The live remote adds the
random-effects effect-weight implementation in `569f68e`/PR 8; that source and
its tests are preserved in the isolated candidate before applying the overlay.

The current P0 findings are not closure evidence for the older work. They add
package-wide beta blockers around PSD kernel geometry, scale-stable rank,
`dkge_infer()` compatibility, data-dependent transport, classification model
selection, fold and validation policy, Sinkhorn/native boundaries, zero-signal
aggregation, and release resolution. Overlapping trialwise concerns remain
nested under the open trust epic `bd-01M03KRRJ77XHPH0DDQ9V0YK6F`; the existing
final trust gate now waits on this P0 epic.

## Dependency and check state

At the audited SHA, `DESCRIPTION` used mutable GitHub `Remotes`, required
non-mainstream `fmridesign` and `fmrireg`, and suggested unavailable
`T4transport`. The full forced-Suggests check stopped at dependency resolution.
The no-Suggests check reached package tests but ended with 12 failures, 10
warnings, 51 skips, 2,361 passes, plus codoc, S3 usage, visible binding, and
undeclared `future`/`withr` findings. The hosted matrix has the desired five
OS/R lanes but currently fails during dependency resolution before package
tests; the latest verified route is `fmrihrf` through `fmridesign`. These are
release-infrastructure failures, not package-test success.

## Phase-1 kernel contract

The beta contract is genuine positive-semidefinite range-space support.
Kernel square roots preserve null eigenvalues as zero; inverse roots use the
Moore-Penrose inverse on the retained numerical range; null directions are
never jittered into metric dimensions. Applied absolute and relative spectral
tolerances, kernel rank/nullity, transformed-moment rank, fold-training rank,
and effective rank are diagnostic data. Requested ranks are capped by all
applicable ranks. A numerical rank of zero must return an explicit typed
zero-rank object where the API can represent it or fail with a precise
zero-rank error before arbitrary components are constructed.

## Editing and staging boundary

All P0 implementation work occurs in `/private/tmp/dkge-beta-p0` on
`codex/public-beta-p0`. Paths are reserved and staged per child mote; the broad
overlay remains unstaged until reconciled. Generated documentation is rebuilt
from roxygen source. No force push, beta tag, or publication action is implied
or authorized by local green evidence.

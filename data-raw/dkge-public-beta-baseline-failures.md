# Public-beta regression baseline

The four `test-beta-*` files were copied without implementation changes into a
pristine `git archive` of
`1fdacb37ef2320a942e08f09eb5078033169478c` and run with R 4.5.1 on macOS.
The run completed with **32 failures, 1 warning, 0 skips, and 10 passes**.

The ten passes are intentional proof/oracle assertions, such as demonstrating
that a loading-derived Sinkhorn operator changes under sign reversal. The 32
failures cover the required absent or violated postconditions: PSD range-space
roots and diagnostics, scale-stable and fold-local rank, inference
compatibility and one-target schema, sign-sensitive transport provenance,
classification selection policy, strict weight/fold/scalar contracts, RNG
preservation, zero-support Sinkhorn behavior, native warm-start validation,
rank-zero aggregate behavior, and the between-subject default.

The separately added `test-release-contracts.R` freezes the release-metadata
and workflow requirements. On the audited SHA it fails because mutable
`Remotes` and `T4transport` remain in core metadata and no sanitizer workflow
exists. This is expected pre-fix evidence, not a green package claim.

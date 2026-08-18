# T4transport differential validation

Date: 2026-08-16  
Mote ticket: `bd-01M0404VQD3AWFTSB3EXG86ZY4`

## Disposition

`T4transport` is **supported as an optional validation reference** on the
tested platform. DKGE does not call it in production code; the package's
native Sinkhorn implementation is compared against it when the suggested
dependency is available.

## Environment

- T4transport 0.1.8, installed from the current CRAN macOS arm64 binary into
  an isolated temporary library.
- R 4.5.1, `aarch64-apple-darwin20`.
- macOS Darwin 23.3.0 arm64.
- DKGE 0.0.0.9000 from the current uncommitted working tree at base commit
  `953d2a6`.

## Evidence

Command:

```sh
R_LIBS=/tmp/dkge-t4-lib Rscript dev/validate-dkge-t4transport.R
```

Both pre-existing fixed fixtures passed the `1e-6` differential and marginal
contracts:

| Fixture | Max plan difference | Objective difference | DKGE max marginal error |
|---|---:|---:|---:|
| Uniform Gaussian blobs | 5.42e-09 | 1.50e-08 | 8.71e-13 |
| Non-uniform 3-D weights | 1.54e-10 | 2.26e-10 | 1.11e-16 |

The existing `test-transport-sinkhorn.R` context also ran with the isolated
library: 16 assertions passed, none failed or skipped.

Machine-readable evidence is in
`inst/extdata/dkge-t4transport-validation.csv` (SHA-256
`03dd12bae2dd2590dbe78423d529734f806db9baa41b7687d402cac51683d163`).

## Scope

This validates agreement for small dense entropic-transport fixtures with
uniform and non-uniform marginals. It does not make T4transport a runtime
dependency, nor does it substitute for DKGE's larger native convergence and
cache tests.


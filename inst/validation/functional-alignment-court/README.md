# Functional-alignment Type-I court

This court measures whether data-adaptive subject-to-reference correspondence
distorts one-sample group inference. The protocol in `frozen-protocol.json` is
the authority: it fixes the null, raw sign action, arms, factors, estimands,
seeds, gates, stopping rule, and claim language before formal results.

Run the exact-oracle pilot first:

```sh
Rscript inst/validation/functional-alignment-court/court.R pilot
```

The pilot writes a runtime record and a mechanically selected formal budget to
`data-raw/functional-alignment-court/`. Only after those files exist may the
formal court run:

```sh
Rscript inst/validation/functional-alignment-court/court.R formal
```

Results are CSV plus a text manifest containing protocol/source hashes,
session information, seeds, worker count, elapsed time, and the predeclared
verdict. The source-tree hash covers the package code actually loaded from a
dirty checkout, rather than relying on Git HEAD alone. A
favorable same-data result never changes its `approximate` label; the
full-pipeline re-estimation arm is the exact null-action comparator for this
court. The formal `B = 31` randomization budget is Monte Carlo, not an
exhaustive sign enumeration.

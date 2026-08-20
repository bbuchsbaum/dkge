# Build an aggregate target from subject-level repeated-measures data

Constructs row-by-feature aggregate targets such as
\`group:task:measure\` cell means from subject-level repeated-measures
matrices. The returned object stores the source metadata required to
recompute aggregates under subject-label permutations or subject
bootstrap resampling. The unit of inference remains subjects; aggregate
rows are never resampled directly.

## Usage

``` r
dkge_aggregate_target(
  values,
  subject_data,
  cell_data = NULL,
  group_vars = NULL,
  cell_vars = NULL,
  subject_id_col = "subject_id",
  cell_id_col = NULL,
  weights = NULL,
  aggregate = c("mean", "weighted_mean"),
  row_sep = ":"
)

# S3 method for class 'dkge_aggregate_target'
print(x, ...)
```

## Arguments

- values:

  Named list of subject matrices. Each matrix is cells x features. When
  every subject matrix has rownames, rows are matched by those names (an
  error if the name sets differ or only some subjects have rownames).
  When subject 1 has colnames, every other subject must have the same
  column-name set; columns are reordered to subject 1.

- subject_data:

  Data frame with one row per subject.

- cell_data:

  Data frame with one row per repeated-measures cell. If \`NULL\`, a
  cell factor named \`cell\` is inferred from matrix row names (or
  \`cell1\`, \`cell2\`, ...) and is kept in the aggregate rows. Pass
  \`cell_vars = character(0)\` to collapse those cells.

- group_vars:

  Character vector of subject-level grouping variables.

- cell_vars:

  Character vector of cell-level variables. Defaults to the inferred
  \`cell\` column when \`cell_data\` is \`NULL\`, otherwise to all
  columns in \`cell_data\` except \`cell_id_col\`.

- subject_id_col:

  Subject identifier column in \`subject_data\`.

- cell_id_col:

  Optional cell identifier column in \`cell_data\`.

- weights:

  Optional subject or subject-by-cell weights.

- aggregate:

  Aggregation method. \`"mean"\` ignores \`weights\`;
  \`"weighted_mean"\` uses them.

- row_sep:

  Separator used to build aggregate row labels. No level of a
  \`group_vars\`/\`cell_vars\` variable may contain it, otherwise two
  distinct aggregate rows could be given the same row ID; violations are
  rejected.

- x:

  A \`dkge_aggregate_target\` object to print.

- ...:

  Unused; present for S3 method compatibility.

## Value

Object of class \`dkge_aggregate_target\`.

## Examples

``` r
set.seed(1)
values <- lapply(1:4, function(i) {
  M <- matrix(rnorm(2 * 3), 2, 3)
  dimnames(M) <- list(c("c1", "c2"), paste0("f", 1:3))
  M
})
names(values) <- paste0("s", 1:4)
subject_data <- data.frame(subject_id = names(values),
                           grp = c("A", "A", "B", "B"))
target <- dkge_aggregate_target(values, subject_data, group_vars = "grp")
nrow(target$Y)
#> [1] 4
```

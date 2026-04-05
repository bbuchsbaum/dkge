# DKGE minimal theme for ggplot2 outputs

Produces a light-weight theme used by the DKGE plotting helpers. Adjust
\`base_size\` / \`base_family\` to customise font size or typeface.

## Usage

``` r
theme_dkge(base_size = 12, base_family = "")
```

## Arguments

- base_size:

  Base font size.

- base_family:

  Base font family.

## Value

A ggplot2 theme object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point() +
    theme_dkge()
}
```

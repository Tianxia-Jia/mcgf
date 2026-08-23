# Calculate general stationary correlation.

Calculate general stationary correlation.

## Usage

``` r
.cor_stat(base, lagrangian, par_base, par_lagr, lambda, base_fixed = FALSE)
```

## Arguments

- base:

  Base model, `sep` or `fs` for now. Or correlation matrix/array.

- lagrangian:

  Lagrangian model, `none`, `lagr_tri`, `lagr_askey`, or `lagr_exp`.

- par_base:

  Parameters for the base model (symmetric), used only when
  `base_fixed = FALSE`.

- par_lagr:

  Parameters for the Lagrangian model. Used only when `lagrangian` is
  not `none`.

- lambda:

  Weight of the Lagrangian term, \\\lambda\in\[0, 1\]\\.

- base_fixed:

  Logical; if TRUE, `base` is the correlation.

## Value

Correlations for the general stationary model. Same dimension of `base`
if `base_fixed = FALSE`.

## Details

The general station model, a convex combination of a base model and a
Lagrangian model, has the form \$\$C(\mathbf{h},
u)=(1-\lambda)C\_{\text{Base}}(\mathbf{h}, u)+ \lambda
C\_{\text{Lagr}}(\mathbf{h}, u),\$\$ where \\\lambda\\ is the weight of
the Lagrangian term.

If `base_fixed = TRUE`, the correlation is of the form \$\$C(\mathbf{h},
u)=(1-\lambda)C\_{\text{Base}}+ \lambda C\_{\text{Lagr}}(\mathbf{h},
u),\$\$ where `base` is a correlation matrix/array and `par_base` and
`h` are not used.

When `lagrangian = "none"`, `lambda` must be 0.

## References

Gneiting, T., Genton, M., & Guttorp, P. (2006). Geostatistical
Space-Time Models, Stationarity, Separability, and Full Symmetry. In
C&H/CRC Monographs on Statistics & Applied Probability (pp. 151–175).
Chapman and Hall/CRC.

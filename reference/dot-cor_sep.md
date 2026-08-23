# Calculate correlation for separable model

Calculate correlation for separable model

## Usage

``` r
.cor_sep(spatial, temporal, par_s, par_t)
```

## Arguments

- spatial:

  Pure spatial model, `exp` or `cauchy` for now.

- temporal:

  Pure temporal model, `exp` or `cauchy` for now.

- par_s:

  Parameters for the pure spatial model. Nugget effect supported.

- par_t:

  Parameters for the pure temporal model.

## Value

Correlations for separable model.

## Details

The separable model is the product of a pure temporal model, \\C_T(u)\\,
and a pure spatial model, \\C_S(\mathbf{h})\\. It is of the form
\$\$C(\mathbf{h}, u)=C\_{T}(u)
\left\[(1-\text{nugget})C\_{S}(\mathbf{h})+\text{nugget}
\delta\_{\mathbf{h}=0}\right\],\$\$ where \\\delta\_{x=0}\\ is 1 when
\\x=0\\ and 0 otherwise. Here \\\mathbf{h}\in\mathbb{R}^2\\ and
\\u\in\mathbb{R}\\. Now only exponential and Cauchy correlation models
are available.

## References

Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
Space–Time Data, Journal of the American Statistical Association,
97:458, 590-600.

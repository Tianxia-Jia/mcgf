# Calculate correlation for fully symmetric model

Calculate correlation for fully symmetric model

## Usage

``` r
.cor_fs(nugget, c, gamma = 1/2, a, alpha, beta = 0, h, u)
```

## Arguments

- nugget:

  The nugget effect \\\in\[0, 1\]\\.

- c:

  Scale parameter of `cor_exp`, \\c\>0\\.

- gamma:

  Smooth parameter of `cor_exp`, \\\gamma\in(0, 1/2\]\\.

- a:

  Scale parameter of `cor_cauchy`, \\a\>0\\.

- alpha:

  Smooth parameter of `cor_cauchy`, \\\alpha\in(0, 1\]\\.

- beta:

  Interaction parameter, \\\beta\in\[0, 1\]\\.

- h:

  Euclidean distance matrix or array.

- u:

  Time lag, same dimension as `h`.

## Value

Correlations of the same dimension as `h` and `u`.

## Details

The fully symmetric correlation function with interaction parameter
\\\beta\\ has the form \$\$C(\mathbf{h},
u)=\dfrac{1}{(a\|u\|^{2\alpha} + 1)}
\left((1-\text{nugget})\exp\left(\dfrac{-c\\\mathbf{h}\\^{2\gamma}}
{(a\|u\|^{2\alpha}+1)^{\beta\gamma}}\right)+ \text{nugget}\cdot
\delta\_{\mathbf{h}=\boldsymbol 0}\right),\$\$ where \\\\\cdot\\\\ is
the Euclidean distance, and \\\delta\_{x=0}\\ is 1 when \\x=0\\ and 0
otherwise. Here \\\mathbf{h}\in\mathbb{R}^2\\ and \\u\in\mathbb{R}\\. By
default `beta = 0` and it reduces to the separable model.

## References

Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
Space–Time Data, Journal of the American Statistical Association,
97:458, 590-600.

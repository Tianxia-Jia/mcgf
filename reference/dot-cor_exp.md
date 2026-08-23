# Calculate exponential correlation

Calculate exponential correlation

## Usage

``` r
.cor_exp(x, c, gamma = 1/2, nugget = 0)
```

## Arguments

- x:

  A numeric vector, matrix, or array.

- c:

  Scale parameter, \\c\>0\\.

- gamma:

  Smooth parameter, \\\gamma\in(0, 1/2\]\\. Default is 1/2.

- nugget:

  The nugget effect \\\in\[0, 1\]\\.

## Value

Correlations of the same dimension as `x`.

## Details

The exponential correlation function with scale parameter \\c\\ and
smooth parameter \\\gamma\\ has the form
\$\$C(x)=(1-\text{nugget})\exp(-c\|x\|^{2\gamma})+\text{nugget}\cdot
\delta\_{x=0},\$\$ where \\\delta\_{x=0}\\ is 1 when \\x=0\\ and 0
otherwise.

## References

Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
Geostatistics. Journal of the Royal Statistical Society. Series C
(Applied Statistics), 47(3), 299–350.

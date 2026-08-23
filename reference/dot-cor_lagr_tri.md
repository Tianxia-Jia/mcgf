# Calculate Lagrangian correlation of the triangular form

Calculate Lagrangian correlation of the triangular form

## Usage

``` r
.cor_lagr_tri(v1, v2, k = 2, h1, h2, u)
```

## Arguments

- v1:

  Prevailing wind, u-component.

- v2:

  Prevailing wind, v-component.

- k:

  Scale parameter of \\\\\boldsymbol v\\\\, \\k\>0\\. Default is 2.

- h1:

  Horizontal distance matrix or array.

- h2:

  Vertical distance matrix or array, same dimension as `h1`.

- u:

  Time lag, same dimension as `h1`.

## Value

Correlations of the same dimension as `h1`.

## Details

The Lagrangian correlation function of the triangular form with
parameters \\\boldsymbol v = (v_1, v_2)^\top\in\mathbb{R}^2\\ has the
form \$\$C(\mathbf{h}, u)=\left(1-\dfrac{1}{k\\\boldsymbol v\\}
\left\|\dfrac{\mathbf{h}^\top\boldsymbol v}{\\\boldsymbol v\\}-
u\\\boldsymbol v\\\right\|\right)\_+,\$\$ where \\\\\cdot\\\\ is the
Euclidean distance, \\x\_+=\max(x,0), \mathbf{h} = (\mathrm{h}\_1,
\mathrm{h}\_2)^\top\in\mathbb{R}^2\\, and \\k \> 0\\ is the scale
parameter controlling the magnitude of asymmetry in correlation.

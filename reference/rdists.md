# Generate random distance matrices

Generate random distance matrices

## Usage

``` r
rdists(N, names, scale = 100)
```

## Arguments

- N:

  Number of locations.

- names:

  Names of locations.

- scale:

  Scale of the distance matrices. Default is 100.

## Value

List of signed distances.

## Details

This function generates random distance matrices using `rnorm`. `scale`
controls the scale of the distance matrices.

## Examples

``` r
set.seed(123)
rdists(3)
#> $h
#>           loc1      loc2     loc3
#> loc1   0.00000  33.54875 268.2444
#> loc2  33.54875   0.00000 239.0565
#> loc3 268.24442 239.05652   0.0000
#> 
#> $h1
#>           loc1      loc2      loc3
#> loc1   0.00000 -33.02982 -211.9184
#> loc2  33.02982   0.00000 -178.8886
#> loc3 211.91840 178.88858    0.0000
#> 
#> $h2
#>            loc1       loc2      loc3
#> loc1   0.000000  -5.877934 -164.4557
#> loc2   5.877934   0.000000 -158.5777
#> loc3 164.455660 158.577725    0.0000
#> 
rdists(3, scale = 1)
#> $h
#>          loc1     loc2     loc3
#> loc1 0.000000 2.401467 1.402200
#> loc2 2.401467 0.000000 1.039848
#> loc3 1.402200 1.039848 0.000000
#> 
#> $h1
#>           loc1      loc2       loc3
#> loc1  0.000000 1.7259774  1.1477691
#> loc2 -1.725977 0.0000000 -0.5782084
#> loc3 -1.147769 0.5782084  0.0000000
#> 
#> $h2
#>           loc1      loc2       loc3
#> loc1 0.0000000 -1.669744 -0.8054758
#> loc2 1.6697438  0.000000  0.8642680
#> loc3 0.8054758 -0.864268  0.0000000
#> 
rdists(3, names = LETTERS[1:3])
#> $h
#>          A        B        C
#> A   0.0000 132.1300 387.3512
#> B 132.1300   0.0000 255.3009
#> C 387.3512 255.3009   0.0000
#> 
#> $h1
#>           A         B        C
#> A   0.00000  29.00887 95.66126
#> B -29.00887   0.00000 66.65239
#> C -95.66126 -66.65239  0.00000
#> 
#> $h2
#>           A         B        C
#> A    0.0000  128.9063 375.3530
#> B -128.9063    0.0000 246.4468
#> C -375.3530 -246.4468   0.0000
#> 
```

cdsrmodels
================

This package contains modeling and biomarker analysis functions
created by the [Cancer Data Science](http://cancerdatascience.org/) team.

## Install

``` r
library(devtools)
devtools::install_github("broadinstitute/cdsr_models")
```
NOTE: If you have issues installing the WGCNA package which cdsr_models depends you can fix them by updating to R >4.0. If you don't want to update to R >4.0 contact William for a workaround.

The package can then be loaded by calling

``` r
library(cdsrmodels)
```

## Modeling functions

### discrete\_test

Compares binary features, such as lineage and mutation, running a two class
comparison on the difference in mean response between cell lines with the
feature and without it. Run on response vector `y` and feature matrix `X`

``` r
cdsrmodels::discrete_test(X, y)
```

### lin\_associations

Compares continuous features, such as gene expression, calculating
correlations between response and each feature. Run on feature matrix
`A`, response vector `y`, and an optional matrix of confounders `W`.
Other parameters can also be tuned and are explained in the function
documentation.

``` r
cdsrmodels::lin_associations(A, y, W=NULL)
```

### random\_forest

Fits a random forest to a feature matrix `X` and a response vector `y`
returning estimates of variable importance for each feature, as well as
model level statistics such as R-squared. Other parameters can also be
tuned and are explained in the function documentation.

``` r
cdsrmodels::random_forest(X, y)
```

### random\_forest\_gauss

Analogous to `random_forest` but uses [`gausscov`](https://cran.r-project.org/package=gausscov) for feature selection in each fold.

``` r
cdsrmodels::random_forest_gauss(X, y)
```

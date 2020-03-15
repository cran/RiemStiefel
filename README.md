
<!-- README.md is generated from README.Rmd. Please edit that file -->
Tools for Inference, Learning, and Optimization on Stiefel Manifold
===================================================================

<!-- badges: start -->
<!-- badges: end -->
Stiefel manifold is a set of orthonormal frames in Euclidean space. We provide algorithms for statistical inference, optimization, and learning over the Stiefel manifold.

Installation
------------

You can install the released version of RiemStiefel from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RiemStiefel")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kyoustat/RiemStiefel")
```

Available Functions
-------------------

| function    | description                                         |
|-------------|-----------------------------------------------------|
| `st.mean`   | Frechet mean computation.                           |
| `st.runif`  | random sample generation from uniform distribution. |
| `st.utestR` | test of uniformity.                                 |
| `stopt.SA`  | Simulated Annealing algorithm for optimization.     |

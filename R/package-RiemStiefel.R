#' Tools for Inference, Learning, and Optimization on Stiefel Manifold
#' 
#' Stiefel manifold \eqn{St(p,r)} is the set of all orthonormal \eqn{r}-frames in \eqn{R^p}, 
#' which is indeed a Riemannian manifold. For \eqn{X \in St(p,r)}, it is characterized as
#' \deqn{X^\top X = I_{r \times r}}. We provide algorithms for statistical inference, optimization, and learning over the Stiefel manifold. 
#' In our package, we use a convention to represent each data point on Stiefel manifold \eqn{St(p,r)} as \eqn{(p\times r)} matrix.
#' 
#' @docType package
#' @name RiemStiefel
#' @importFrom utils packageVersion
#' @importFrom RiemBase riemfactory
#' @importFrom RiemBaseExt rstat.frechet
#' @importFrom stats rnorm runif cov
#' @importFrom Rcpp evalCpp
#' @useDynLib RiemStiefel
NULL

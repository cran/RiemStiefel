#' Fréchet Mean on Stiefel Manifold
#' 
#' For manifold-valued data, Fréchet mean is the solution of following cost function,
#' \deqn{\textrm{min}_x \sum_{i=1}^n \rho^2 (x, x_i),\quad x\in\mathcal{M}}
#' for a given data \eqn{\{x_i\}_{i=1}^n} and \eqn{\rho(x,y)} is the geodesic distance 
#' between two points on manifold \eqn{\mathcal{M}}. It uses a gradient descent method 
#' with a backtracking search rule for updating. In the Stiefel manifold case, 
#' analytic formula is not known so we use numerical approximation scheme.
#' 
#' @param x either an array of size \eqn{(p\times r\times n)} or a list of length \eqn{n} whose elements are \eqn{(p\times r)} matrix on Stiefel manifold.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param eps stopping criterion for the norm of gradient.
#' @param parallel a flag for enabling parallel computation with OpenMP.
#' 
#' @return a named list containing
#' \describe{
#' \item{mu}{an estimated mean matrix of size \eqn{(p\times r)}.}
#' \item{variation}{Fréchet variation with the estimated mean.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #              Average Projection with 'iris' dataset
#' #-------------------------------------------------------------------
#' #  For PCA, take half of data from 'iris' data and repeat it 10 times.
#' #  We will compare naive PCA, intrinsic, and extrinsic mean.
#' 
#' data(iris)
#' label = iris$Species           # label information
#' idata = as.matrix(iris[,1:4])  # numeric data
#' ndata = nrow(idata)
#' 
#' # define a function for extracting pca projection
#' pcaproj <- function(X, p){
#'   return(eigen(stats::cov(X))$vectors[,1:p])
#' }
#' 
#' # extract embedding for random samples
#' proj10 = list()
#' for (i in 1:10){
#'    # index for random subsample
#'    rand.now    = base::sample(1:ndata, round(ndata/2)) 
#'    
#'    # PCA via personal tool
#'    proj10[[i]] = pcaproj(idata[rand.now,], p=2)
#' }
#' 
#' # compute intrinsic and extrinsic mean of projection matrices
#' mean.int = st.mean(proj10, type='intrinsic')$mu
#' mean.ext = st.mean(proj10, type='extrinsic')$mu
#' 
#' # compute 2-dimensional embeddings
#' fproj = pcaproj(idata, p=2)
#' f2 = idata%*%fproj     # PCA with full data
#' i2 = idata%*%mean.int  # projection with intrinsic mean
#' e2 = idata%*%mean.ext  # projection with extrinsic mean
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(f2, cex=0.5, pch=19, col=label, main='full PCA')
#' plot(i2, cex=0.5, pch=19, col=label, main='intrinsic projection')
#' plot(e2, cex=0.5, pch=19, col=label, main='extrinsic projection')
#' par(opar)
#' 
#' @references 
#' \insertRef{bhattacharya_nonparametric_2012}{RiemStiefel}
#' 
#' @export
st.mean <- function(x, type=c("intrinsic","extrinsic"), eps=1e-6, parallel=FALSE){
  ############################################################
  # Preprocessing
  x      = RiemBase::riemfactory(check_data(x), name="stiefel")
  mytype = match.arg(type)
  
  myeps      = as.double(eps)
  myparallel = as.logical(parallel)
  
  ############################################################
  # Computation
  output = RiemBaseExt::rstat.frechet(x, type=mytype, int.eps=myeps, parallel=myparallel)
  return(output)
}
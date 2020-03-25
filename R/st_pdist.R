#' Pairwise Distance for Data on Stiefel Manifold
#' 
#' For data on Stiefel manifold \eqn{x_1,x_2,\ldots,x_N \in St(r,p)}, compute pairwise distances \eqn{d(x_i,x_j)} via geodesic (\code{"intrinsic"}) distance 
#' or embedded Euclidean \code{"extrinsic"} metric. Since computing geodesic has no closed-form expressions, 
#' it relies on a numerical approximation which may incur heavier computational burden.
#' 
#' @param x either an array of size \eqn{(p\times r\times N)} or a list of length \eqn{N} whose elements are \eqn{(p\times r)} matrix on Stiefel manifold.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param as.dist a logical; \code{TRUE} to return a \code{\link[stats]{dist}} object or \code{FALSE} to return an \eqn{(N\times N)} symmetric matrix.
#' 
#' @return a \code{\link[stats]{dist}} object or \eqn{(N\times N)} symmetric matrix depending on \code{as.dist}.
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #      Generate a dataset with two types of Stiefel elements
#' #-------------------------------------------------------------------
#' #  group1 : first four columns of (8x8) identity matrix + noise
#' #  group2 : last  four columns of (8x8) identity matrix + noise
#' 
#' mydata = list()
#' sdval  = 0.05
#' diag8  = diag(8)
#' for (i in 1:10){
#'   mydata[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' for (i in 11:20){
#'   mydata[[i]] = qr.Q(qr(diag8[,5:8] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' 
#' ## compare 'intrinsic' and 'extrinsic' distances
#' dint = st.pdist(mydata, type="intrinsic", as.dist=FALSE)
#' dext = st.pdist(mydata, type="extrinsic", as.dist=FALSE)
#' 
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(dint[,20:1], main="intrinsic")
#' image(dext[,20:1], main="extrinsic")
#' par(opar)
#' 
#' @author Kisung You
#' @export
st.pdist <- function(x, type=c("intrinsic","extrinsic"), as.dist=TRUE){
  ############################################################
  # Preprocessing
  x       = RiemBase::riemfactory(check_data(x), name="stiefel")
  mytype  = match.arg(type)
  retdist = as.logical(as.dist)
  
  ############################################################
  # Computation
  output = RiemBaseExt::rstat.pdist(x, type=mytype, as.dist=retdist)
  
  ############################################################
  # Return
  if (retdist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}
#' Hierarchical Agglomerative Clustering on Stiefel Manifold
#' 
#' Given the \code{type} of distance measure and agglomeration scheme \code{method}, \code{gr.hclust} performs hierarchical clustering on 
#' Grassmann manifold using \pkg{fastcluster} package, which returns the same object as \pkg{stats} package's implementation while providing more efficient computation. 
#' See \code{\link[fastcluster]{hclust}} for more details.
#' 
#' @param x either an array of size \eqn{(p\times r\times N)} or a list of length \eqn{N} whose elements are \eqn{(p\times r)} matrix on Stiefel manifold.
#' @param type type of distance measure; \code{"intrinsic"} or \code{"extrinsic"}.
#' @param method he agglomeration method to be used. This must be (an unambiguous abbreviation of) one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details. 
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
#' ## try hierarchical clustering
#' #  compare 'intrinsic' and 'extrinsic' distance types
#' #  and use 'single' hclust option.
#' hint = st.hclust(mydata, type="intrinsic", method="single")
#' hext = st.hclust(mydata, type="extrinsic", method="single")
#' 
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(hint, main="intrinsic")
#' plot(hext, main="extrinsic")
#' par(opar)
#' 
#' @author Kisung You
#' @export
st.hclust <- function(x, type=c("intrinsic","extrinsic"),
                      method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
                                 "centroid", "median"),
                      members=NULL){
  ############################################################
  # Preprocessing
  x         = RiemBase::riemfactory(check_data(x), name="stiefel")
  mytype    = match.arg(type)
  mymethod  = match.arg(method)
  mymembers = members
  
  ############################################################
  # Compute Distance and Apply Hclust
  pdmat = RiemBaseExt::rstat.pdist(x, type=mytype, as.dist=TRUE)
  hcout = RiemBaseExt::rclust.hclust(pdmat,method=mymethod,members=mymembers)
  return(hcout)
}
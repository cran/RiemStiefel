#' Test of Uniformity via Rayleigh Statistic on Stiefel Manifold
#' 
#' This function is for hypothesis testing on Stiefel manifold \eqn{St(p,r)} whether the given data is uniformly distributed or not. 
#' We provide two options (original and modified) for Rayleigh-type statistics, which both follow Chi-squared distribution of degrees of 
#' freedom \eqn{pr}.
#' 
#' @param x either an array of size \eqn{(p\times r\times n)} or a list of length \eqn{n} whose elements are \eqn{(p\times r)} matrix on Stiefel manifold.
#' @param method \code{"original"} for conventional Rayleigh statistic or \code{"modified"} for better order of error.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' ## Test of Uniformity for 100 samples from St(10,5)
#' #  Data Generation
#' mydat = st.runif(n=100, p=10, r=5, rtype='list')
#' 
#' #  Run Tests using two methods
#' st.utestR(mydat, method='original')
#' st.utestR(mydat, method='modified')
#' 
#' \dontrun{
#' ## empirical Type 1 error using the same setting as above.
#' niter   = 10000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = st.runif(n=100, p=10, r=5, rtype='list')
#'   counter[i] = ifelse(st.utestR(X)$p.value < 0.05, 1, 0)
#'   print(paste0("iteration ",i,"/10000 complete..."))
#' }
#' 
#' ## print the result
#' print(paste0("* empirical Type 1 error for 'st.utestR': ",round(sum(counter/niter),5)))
#' }
#' 
#' @references 
#' \insertRef{mardia_directional_1999}{RiemStiefel} 
#' 
#' @export
st.utestR <- function(x, method=c("Original","Modified")){
  ##############################################################
  # Preprocessing
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  x = check_data(x)  # now 3d array
  p = dim(x)[1]      # r-frame in R^p with n observations
  r = dim(x)[2]
  n = dim(x)[3]
  method = tolower(method)
  alldip = c("original","modified")
  method = match.arg(method, alldip)
  
  ##############################################################
  # Preliminary Compute
  # xbar = cppaux_mean3d(x)
  xbar = array(0,c(p,r))
  for (i in 1:p){
    for (j in 1:r){
      xbar[i,j] = mean(as.vector(x[i,j,]))
    }
  }
  S    = p*n*sum(diag(t(xbar)%*%xbar))
  
  
  ##############################################################
  # Main Computation
  if (aux_strcmp(method, "original")){
    hname   = "Rayleigh Test of Uniformity of Stiefel Manifold"
    thestat = S
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  } else if (aux_strcmp(method, "modified")){
    hname   = "Modified Rayleigh Test of Uniformity of Stiefel Manifold"
    term1   = 1/(2*n)
    term2   = 1 - (S/((p*r) + 2))
    thestat = S*(1 - term1*term2)
    pvalue  = stats::pchisq(thestat, df=as.integer(p*r), lower.tail=FALSE)
  }
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  Ha      = paste("data is not uniformly distributed on St(",p,",",r,").",sep="")
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}



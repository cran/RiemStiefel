#' Generate Random Samples from Uniform Distribution on Stiefel Manifold
#' 
#' It generates \eqn{n} random samples from Stiefel manifold \eqn{St(p,r)} 
#' according to the procedure described in the reference.
#' 
#' @param n number of samples to be generated.
#' @param p original dimension (of the ambient space).
#' @param r dimension of the frame.
#' @param rtype return type; either 3d-array (\code{"array"}) or list (\code{"list"}).
#' 
#' @return a length \eqn{n} list or 3d array of size \eqn{(p,r,n)}.
#' 
#' @examples 
#' ## let's simply draw 3 times from St(10,5)
#' dat3 = st.runif(3, 10, 5, rtype="array")
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' image(dat3[,,1], main="sample 1")
#' image(dat3[,,2], main="sample 2")
#' image(dat3[,,3], main="sample 3")
#' par(opar)
#' 
#' @references 
#' \insertRef{chikuse_statistics_2003}{RiemStiefel} 
#' 
#' @export
st.runif <- function(n, p, r, rtype=c("list","array")){
  ##############################################################
  # Preprocessing
  n = round(n)
  p = round(p)
  r = round(r)
  
  ##############################################################
  # Main Computation
  rtype = tolower(rtype)
  alltype = c("list","array")
  rtype   = match.arg(rtype, alltype)
  if (aux_strcmp(rtype,"list")){
    output = list()
    for (i in 1:n){
      output[[i]] = st.runif.single(p,r)
    }
  } else if (aux_strcmp(rtype,"array")){
    output = array(0,c(p,r,n))
    for (i in 1:n){
      output[,,i] = st.runif.single(p,r)
    }
  }
  
  ##############################################################
  # Return
  return(output)
}

#------------------------------------------------
#' @keywords internal
#' @noRd
st.runif.single <- function(p,r){
  X = matrix(stats::rnorm(p*r),ncol=r)
  H = aux_halfinv((t(X)%*%X))
  return(X%*%H)
}

# source : https://math.stackexchange.com/questions/3097862/uniform-distribution-on-stiefel
#          Chikuse (2003) Statistics on Special Manifolds
#          Thm 2.2.1
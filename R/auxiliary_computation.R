## Auxiliary Computation
# (1) aux_strcmp  : MATLAB strcmp
# (2) aux_halfinv : A^(-1/2)





# (1) aux_strcmp ----------------------------------------------------------
#' @keywords internal
#' @noRd
aux_strcmp <- function(s1, s2) {
  if (!is.vector(s1, mode="character") || !is.vector(s1, mode="character"))
    stop("Arguments 's1' and 's2' must be character vectors.")
  
  if (length(s1) == length(s2)){
    return(all(s1 == s2))
  } else {
    return(FALSE)
  }
}

# (2) aux_halfinv ---------------------------------------------------------
#' @keywords internal
#' @noRd
aux_halfinv <- function(A){
  eigA = base::eigen(A)
  return((eigA$vectors %*% base::diag(1/sqrt(eigA$values)) %*% t(eigA$vectors)))
}

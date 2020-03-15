## checkers
## (1) check_data : return


# (1) check_data
#' @keywords internal
#' @noRd
check_data <- function(x){
  wow = RiemBase::riemfactory(x, name="stiefel")
  N = length(wow$data)
  n = wow$size[1]
  p = wow$size[2]
  
  output = array(0,c(n,p,N))
  for (i in 1:N){
    output[,,i] = wow$data[[i]]
  }
  return(output)
}
# x = list()
# for (i in 1:3){
#   xx = matrix(rnorm(5*3),ncol=3)
#   x[[i]] = qr.Q(qr(xx))
# }
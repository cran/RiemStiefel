#' Optimization over Stiefel Manifold with Simulated Annealing
#' 
#' Simulated Annealing is a black-box, derivative-free optimization algorithm 
#' that iterates via stochastic search in the neighborhood of current position.
#' 
#' @param func a function to be \emph{minimized}.
#' @param size a length-2 vector containing dimension information \eqn{(p,r)}.
#' @param n.start number of runs, i.e., algorithm is executed \code{n.start} times.
#' @param stepsize  size of random walk on each component.
#' @param maxiter maximum number of iterations for each run.
#' @param cooling triplet for cooling schedule. See the section for the usage.
#' @param init.val if \code{NULL}, starts from a random point. Otherwise, a Stiefel matrix of size \eqn{(p,r)} should be provided for fixed starting point.
#' @param print.progress a logical; \code{TRUE} to show 
#' 
#' @return a named list containing: \describe{
#' \item{cost}{minimized function value.}
#' \item{solution}{a \eqn{(p\times r)} matrix that attains the \code{cost}.}
#' \item{accfreq}{frequency of acceptance moves.}
#' }
#' 
#' @examples 
#' ## Optimization for eigen-decomposition
#' #  Let's find top-3 eigenvalues 
#' set.seed(121)                         # set seed
#' A = cov(matrix(rnorm(100*5), ncol=5)) # define covariance
#' myfunc <- function(p){                # cost function
#'   return(sum(-diag(t(p)%*%A%*%p)))
#' } 
#' 
#' #  Solve the optimization problem
#' Aout = stopt.SA(myfunc, size=c(5,3), n.start=40, maxiter=500)
#' 
#' #  Compute 3 Eigenvalues
#' #  1. use computed basis
#' abase   = Aout$solution
#' eig3sol = sort(diag(t(abase)%*%A%*%abase), decreasing=TRUE)
#'
#' #  2. use 'eigen' function
#' eig3dec = sort(eigen(A)$values, decreasing=TRUE)[1:3]
#' 
#' \donttest{
#' #   Visualize
#' opar <- par(no.readonly=TRUE)
#' yran = c(min(min(eig3sol),min(eig3dec))*0.95,
#'          max(max(eig3sol),max(eig3dec))*1.05)
#' plot(1:3, eig3sol, type="b", col="red",  pch=19, ylim=yran,
#'      xlab="index", ylab="eigenvalue", main="compare top 3 eigenvalues")
#' lines(1:3, eig3dec, type="b", col="blue", pch=19)
#' legend(1, 1, legend=c("optimization","decomposition"), col=c("red","blue"),
#'        lty=rep(1,2), pch=19)
#' par(opar)
#' }
#' 
#' @export
stopt.SA <- function(func, size, n.start=10, stepsize=0.1, maxiter=100, cooling=c("exponential",10,0.9), init.val=NULL, 
                     print.progress=FALSE){
  ##############################################################
  # Preprocessing
  # 1. function
  if (!is.function(func)){
    stop("* stopt.SA : an input 'func' should be a function.")
  }
  # 2. size
  if ((length(init.val)==0)&&(is.null(init.val))){ # not given
    if ((!is.vector(size))||(length(size)!=2)){
     stop("* stopt.SA : 'size' should be a vector of length 2.") 
    }  
    # for Stiefel Manifold Only
    n = round(size[1])
    p = round(size[2])
    initflag = FALSE
  } else { # if given, override size information
    n = nrow(init.val)
    p = ncol(init.val)
    initflag = TRUE 
  }
  # 3. other parameters
  my.nstart      = round(n.start)
  my.stepsize    = as.double(stepsize)
  my.temperature = sa_check_cooling(cooling, maxiter)
  
  ##############################################################
  # Main Run
  if (isTRUE(initflag)){
    out.now = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
    if (print.progress){
      print(paste0("* stopt.SA : iteration 1/",n.start," complete.."))
    }
  } else {
    init.val = base::qr.Q(base::qr(matrix(stats::rnorm(n*p),nrow=n)))
    out.now  = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
    if (print.progress){
      print(paste0("* stopt.SA : iteration 1/",n.start," complete.."))
    }
  }
  if (my.nstart > 1){
    for (it in 1:(my.nstart-1)){
      if (isTRUE(initflag)){
        out.tmp = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
      } else {
        init.val = base::qr.Q(base::qr(matrix(stats::rnorm(n*p),nrow=n)))
        out.tmp  = sa_engine_Stiefel(func, init.val, my.temperature, my.stepsize)
      }
      if (out.tmp$cost <= out.now$cost){ # update with a better one
        out.now = out.tmp
      }
      if (print.progress){
        print(paste0("* stopt.SA : iteration ",it+1,"/",n.start," complete.."))
      }
    }
  }
  
  ##############################################################
  # Report
  return(out.now)
}


# SA functions ------------------------------------------------------------
# (1) sa_check_cooling  : check the cooling schedule
#                         c("exponential", C, alpha)
#                         c("logarithmic", C, D)
#                         c("turnpike",    C, D)
# (2) sa_engine_Stiefel : working version of the function for Stiefel Manifold

#' @keywords internal
#' @noRd
sa_check_cooling <- function(coolschedule, itermax){
  if ((!is.vector(coolschedule))||(length(coolschedule)!=3)){
    stop("* stopt.SA : 'cooling' schedule must be a vector of length 3.")
  }
  themethod = coolschedule[1]
  itermax   = round(itermax)
  C         = as.double(coolschedule[2])
  if (C <= 0){
    stop("* stopt.SA : 'C' should be a nonnegative number.")
  }
  if (all(tolower(themethod)=="exponential")){
    alpha = as.double(coolschedule[3])
    if ((alpha <=0)||(alpha >=1)){
      stop("* stopt.SA : when the cooling schedule is 'exponential', 'alpha' should be a value in (0,1).")
    }
    outvec = C*(alpha^(1:itermax))
  } else if (all(tolower(themethod)=="logarithmic")){
    D = as.double(coolschedule[3])
    if (D <= 0){
      stop("* stopt.SA : when the cooling schedule is 'logarithmic', 'D' should be a positive real number.")
    }
    outvec = C/(log((1:itermax)+D))
  } else if (all(tolower(themethod)=="turnpike")){
    D = as.double(coolschedule[3])
    if (D <= 1){
      stop(" stopt.SA : when the cooling schedule is 'turnpike', 'D' should be a real number larger than 1.")
    }
    outvec = C*((D-1)/(log((1:itermax))))
  }
  
  if (any(is.infinite(outvec))){
    outvec[is.infinite(outvec)] = 2*max(outvec[!is.infinite(outvec)])
  }
  return(outvec)
}

# (2) sa_engine_Stiefel : working version of the function for Stiefel Manifold
#' @keywords internal
#' @noRd
sa_engine_Stiefel <- function(func, init.mat, temparature, stepsize){
  # initialization
  n = nrow(init.mat)
  p = ncol(init.mat)
  Eold    = func(init.mat)
  maxiter = length(temparature)
  
  # iteration
  sol.old = init.mat
  count   = 0
  for (k in 1:maxiter){
    # 1. generation unique to Stiefel 
    sol.tmp = sol.old + matrix(stats::rnorm(n*p, sd=stepsize), nrow=n)
    sol.tmp = base::qr.Q(base::qr(sol.tmp))
    # 2. energy evaluation & current temperature
    Etmp = func(sol.tmp)
    # 3. decision branching
    if (Etmp <= Eold){ # unconditional accept
      Enew    = Etmp
      sol.new = sol.tmp
      count   = count + 1
    } else {           # conditional accept
      Tk = temparature[k]
      tprob = exp((-(Etmp-Eold))/Tk)
      if (as.double(stats::runif(1)) < tprob){
        Enew    = Etmp
        sol.new = sol.tmp
        count   = count + 1
      } else {
        Enew    = Eold
        sol.new = sol.old
      }
    }
    # 4. update
    Eold    = Enew
    sol.old = sol.new
  }
  
  # report the run
  output = list()
  output$cost     = Eold
  output$solution = sol.old
  output$accfreq  = (count/maxiter)
  return(output)
}


  
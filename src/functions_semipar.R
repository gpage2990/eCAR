
# Function to create the b-spline basis from Brian
mybs <- function(x,L){
  library(splines)
  n  <- length(x)
  x  <- x/max(x)
  if(L==1){B<-matrix(1,n,1)}
  if(L==2){B<-cbind(1-x,x)}
  if(L==3){B<-cbind(1-x,1-4*(x-0.5)^2,x)}
  if(L>=4){B<-bs(x,df=L,intercept=TRUE)}
  return(B)}

# Function to create the b-spline basis (P-spline method, Eilers 1996)
bspline <- function(x, xl=min(x), xr=max(x), ndx, bdeg) {
  library(splines)
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
  B
}

rescale.row = function(A,vec){
  t(apply(A, 1, function(x, vec) x*vec, vec=vec))
}



## inla.rgenric.define model for the Gaussian case (eq 25, sec 4.3 overleaf paper)
myrgeneric.eigenLEROUX = function (cmd = c("graph", "Q", "mu", "initial",
                                           "log.norm.const", "log.prior", "quit"),
                                   theta = NULL)
{
  interpret.theta = function() {
    return(list(prec = exp(theta[1L]),
                r = exp(theta[2L])/(1+exp(theta[2L])),
                lam = exp(theta[3L])/(1+exp(theta[3L]))))
  }
  graph = function() {
    return(Q())
  }
  Q = function() {
    prec = interpret.theta()$prec
    r = interpret.theta()$r
    lam = interpret.theta()$lam
    Q = inla.as.sparse(diag(( r/(prec*(1-lam+lam*v)) + (1-r)/prec )^(-1)))
    return(Q)
  }
  mu = function() {
    return(numeric(0))
  }
  log.norm.const = function() {
    prec = interpret.theta()$prec
    r = interpret.theta()$r
    lam = interpret.theta()$lam
    a = ( r/(prec*(1-lam+lam*v)) + (1-r)/prec )^(-1)
    #return(sum(log(a)))
    val = -0.5*n * log(2*pi) + 0.5*sum(log(a))
    return(val)
    #    return(numeric(0))
  }
  log.prior = function() {
    prec = interpret.theta()$prec
    val = inla.pc.dprec(prec,u.prec,alpha.prec,log=T) + theta[1L] +
      theta[2L] - 2*log(exp(theta[2L]) + 1) +  # unif
      theta[3L] - 2*log(exp(theta[3L]) + 1)    # unif
    return(val)
  }
  initial = function() {
    return(c(0,1,1))
  }
  quit = function() {
    return(invisible())
  }
  if (is.null(theta)) theta = initial()
  val = do.call(match.arg(cmd), args = list())
  return(val)
}

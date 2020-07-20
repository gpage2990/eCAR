
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


### description of 'semipar.eCARglm.Leroux'
# E: offset (e.g. for covid data E=scale(pop,center=F))
# y: data vector (possibly including NA) of length n
# x: exposure data vector of length n
# W: matrix n x n with 1 if i~j otherwise 0
# X.conf: matrix with confounders, nrow(conf)=n

# TODO: include random effects
# X.random: indicator matrix for possible iid random effect, nrow(conf)=n
# recall this even if not needed; conf = tapply(X.conf,rep(1:ncol(X.conf),each=nrow(X.conf)),function(i) i)

semipar.eCARglm.Leroux = function(y, x, W, E, C=NA,
                               model="Poisson",
                               L=20, pcprior.sd=c(0.1,1), s2=2,
                               eval.fineGrid=FALSE,
                               verbose=FALSE){
  library(INLA)
  # compute spectral decomposition of the structure matrix 'R'
  n = length(y)
  M <- rowSums(W)
  R <- diag(M)-W
  Eigdec <- eigen(R)
  G <- Eigdec$vec
  v <- Eigdec$val

  X.conf <- C
  # compute cubic b-spline basis
  # Eilers basis:
  B = bspline(x=v, ndx=L-3, bdeg=3)  # Eilers basis
  if (eval.fineGrid) B.pred = bspline(x=seq(min(v),max(v),length.out=1000), ndx=L-3, bdeg=3)  else B.pred=B
  # Brian basis:
  #B = mybs(x=v, L=L)
  #if (eval.fineGrid) B.pred = mybs(x=seq(min(v), max(v), length.out=1000), L=L) else B.pred=B

  Z.tilde = matrix(NA, nrow = n, ncol = L)
  for(i in 1:L){
    Z.tilde[,i] =  rescale.row(G, B[,i]) %*% (t(G) %*% x)
  }

  if (all(is.na(X.conf))) {
    stk = inla.stack(
      data=list(y=y, E = E),
      A=list(1, Z.tilde, diag(n)),
      effects=list(intercept=rep(1,n),
                   id.beta=1:L,
                   id.z = 1:n))
    formula = y ~ -1 + intercept +
      f(id.beta,
        model = "rw1",
        constr=FALSE, scale.model = TRUE,
        hyper = list(prec=list(
          prior="pc.prec",
          param=c(pcprior.sd[1]/0.31, 0.01)))) +
      f(id.z, model = "generic1",
        Cmatrix = Diagonal(x=1, n=n)-inla.as.sparse(R),
        hyper = list(prec=list(
          prior="pc.prec",
          param=c(pcprior.sd[2]/0.31, 0.01))))
  } else {
      stk = inla.stack(
      data=list(y=y, E = E),
      A=list(1, Z.tilde, diag(n)),
      effects=list(conf= cbind(rep(1,n), X.conf),
                   id.beta=1:L,
                   id.z = 1:n))
      formula = y ~ -1 + conf +
        f(id.beta,
          model = "rw1",
          constr=FALSE, scale.model = TRUE,
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[1]/0.31, 0.01)))) +
        f(id.z, model = "generic1",
          Cmatrix = Diagonal(x=1, n=n)-inla.as.sparse(R),
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[2]/0.31, 0.01))))
  }

  #  likelihood Negative Binomial 'NegBin'
  if (model=="Negative Binomial"){
    r = inla(formula,
            data = inla.stack.data(stk),
            family='nbinomial', E=E,
            control.family = list(
              variant=0,
              hyper=list(theta= list(
                prior="normal", param=c(0,1/s2^2)))),
            control.predictor = list(A = inla.stack.A(stk), compute=TRUE),
            control.compute = list(dic=TRUE, config=TRUE),
            verbose=verbose)
  }


  #  likelihood Binomial 'Bin'
  if (model=="Binomial"){
    r = inla(formula,
             data = inla.stack.data(stk),
             family='binomial', Ntrials=E,
             control.predictor = list(A = inla.stack.A(stk), compute=TRUE),
             control.compute = list(dic=TRUE, config=TRUE),
             verbose=verbose)
  }

  # likelihood 'Poisson'
  if (model=="Poisson"){
    r = inla(formula,
             data = inla.stack.data(stk),
             family='poisson', E=E,
             control.predictor = list(A = inla.stack.A(stk), compute=TRUE),
             control.compute = list(dic=TRUE, config=TRUE),
             verbose=verbose)
  }

  # sample from the posterior and compute posterior mean curve
  sample.tmp = inla.posterior.sample(1000,r)
  splinecoefs = matrix(nrow=L, ncol=1000)
  ind.beta =  c(paste("id.beta:",1:L, sep=""))
  for(j in 1:1000){
    splinecoefs[,j] = (sample.tmp[[j]]$latent)[ind.beta,]
  }
  beta.mn = apply(exp(B.pred%*%splinecoefs),1,mean)
  beta.q025 = apply(exp(B.pred%*%splinecoefs),1,quantile, 0.025)
  beta.q975 = apply(exp(B.pred%*%splinecoefs),1,quantile, 0.975)

  # return results
  if (!eval.fineGrid) {
    out = list(omega=v, beta.mn=beta.mn, beta.q025=beta.q025, beta.q975=beta.q975)
  } else {
    out = list(omega=seq(min(v),max(v),length.out=1000),
               beta.mn=beta.mn, beta.q025=beta.q025, beta.q975=beta.q975)
  }
  return(out)
}




### TODO:
  # 1) add the Gaussian case in the same function 'semipar.eCAR.Leroux'
  # 2) try to fit the covid data with the confounders...see if it works the stack and inla call
        # SOLVE ISSUE ABOUT IND.LATENT.....





# Wrapper to fit the joint Leroux model
par.eCAR.Leroux <- function(y,x,W,E=NULL,C=NULL,model="Gaussian",
                    joint_prior_lamx_lamz = FALSE,
                    lamx.fix.val = NULL, sig2x.fix.val = NULL,
                    m=0,s2=10,alamx=1,blamx=1,alamz=1,blamz=1,
                    asig=1,bsig=1,atau=1,btau=1,asigx=1,bsigx=1,
                    mb0=0,s2b0=100,me=0,s2e=100,mx=0, s2x=100,
                    tau_cand_sd=1, sig2_cand_sd=1,
                    draws=10000,burn=5000,thin=5, verbose=TRUE){
  # W - is the neighborhood matrix
  # C - is the matrix of confounder covariates that are included not

  # Note that y and x are both measure from spatial domain.
  # That is, they haven't been projected to spectral domain.
  #
  # model specifies the data model with options being
  # 1 - Gaussian
  # 2 - Poisson
  # 3 - Binomial
  # 4 - Negative Binomial (Poisson with overdispersion)


  # HPD calculator (see TeachingDemos)
  emp.hpd <- function(x, conf=0.95){
    conf <- min(conf, 1-conf)
    n <- length(x)
    nn <- round( n*conf )
    x <- sort(x)
    xx <- x[ (n-nn+1):n ] - x[1:nn]
    m <- min(xx)
    nnn <- which(xx==m)[1]
    return( c( x[ nnn ], x[ n-nn+nnn ] ) )
  }


  cat("A", model, "model is being fit \n")
  out <- NULL

  nout <- (draws-burn)/thin

  nobs <- length(y)

  M <- diag(apply(W,1,sum))
  R <- M-W

  eigendec.R <- eigen(R)
  evals <- eigendec.R$values
  evecs <- eigendec.R$vectors

  ystar <- t(evecs) %*% y
  xstar <- t(evecs) %*% x

  ncov <- ncol(C)
  if(is.null(C)){
    ncov <- 0
    C <- 0
    Cstar <- 0
  } else {
    Cstar <- t(evecs) %*% as.matrix(C)
  }

  updateXparms <- FALSE
  if(is.null(lamx.fix.val) | is.null(sig2x.fix.val)){
    updateXparms <- TRUE
    sig2x.fix.val <- 1
    lamx.fix.val <- 0.5
  }


  modelPriors = c(m, s2, alamx, blamx, alamz, blamz, asig, bsig,
                  atau, btau, asigx, bsigx, mb0, s2b0, me, s2e, mx, s2x)
  MHsd <- c(tau_cand_sd, sig2_cand_sd)
  beta0 <- beta <- alpha <- tau <- sig2 <- rep(1, nout)
  sig2x <- rep(sig2x.fix.val,nout)
  lamx <- lamz <- rep(lamx.fix.val, nout)
  theta <- nb_r <- matrix(0, nrow=nout, ncol=nobs)
  eta <- matrix(0, nrow=nout, ncol=ncov)


  if(model=="Gaussian"){
    # I first center the X and the Y to set the intercept to zero
    ystar <- ystar - mean(ystar)
    xstar <- xstar - mean(xstar)

    C.out <- .C("mcmcloop_leroux_gauss",
                as.integer(draws), as.integer(burn), as.integer(thin),
                as.integer(nobs), as.double(ystar), as.double(xstar),
                as.double(evals), as.double(t(Cstar)), as.integer(ncov),
                as.double(modelPriors), as.double(MHsd),
                as.integer(verbose), as.integer(joint_prior_lamx_lamz),
                as.integer(updateXparms),
                beta.out=as.double(beta), alpha.out=as.double(alpha),
                tau.out=as.double(tau), sig2x.out=as.double(sig2x),
                lamx.out=as.double(lamx), lamz.out=as.double(lamz),
                sig2.out=as.double(sig2), eta.out=as.double(eta))


  }
  if(model!="Gaussian"){
    if(model == "Poisson") modelnum = 1
    if(model == "Binomial") modelnum = 2
    if(model == "Negative Binomial") modelnum = 3

    C.out <- .C("mcmcloop_leroux_GLM",
                as.integer(draws), as.integer(burn), as.integer(thin),
                as.integer(nobs), as.double(y), as.double(x), as.double(E),
                as.double(evals), as.double(t(evecs)),
                as.double(t(W)), as.double(t(C)), as.integer(ncov),
                as.integer(modelnum), as.double(modelPriors),
                as.double(MHsd), as.integer(verbose),
                as.integer(joint_prior_lamx_lamz),
                as.integer(updateXparms),
                beta.out=as.double(beta), alpha.out=as.double(alpha),
                tau.out=as.double(tau), sig2x.out=as.double(sig2x),
                lamx.out=as.double(lamx), lamz.out=as.double(lamz),
                theta.out=as.double(theta), beta0.out=as.double(beta0),
                eta.out=as.double(eta), r.out=as.double(nb_r))

  }



  out$beta <- matrix(C.out$beta.out, nrow=nout, byrow=TRUE)
  out$gamma <- matrix(C.out$alpha.out, nrow=nout, byrow=TRUE)
  out$tau <- matrix(C.out$tau.out, nrow=nout, byrow=TRUE)
  out$sig2x <- matrix(C.out$sig2x.out, nrow=nout, byrow=TRUE)
  out$lamx <- matrix(C.out$lamx.out, nrow=nout, byrow=TRUE)
  out$lamz <- matrix(C.out$lamz.out, nrow=nout, byrow=TRUE)
  if(model=="Gaussian") out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
  out$rho <- matrix((out$gamma*sqrt(out$sig2x))/sqrt(out$tau + out$gamma^2*out$sig2x), nrow=nout, byrow=TRUE)
  out$sig2z <- matrix(out$tau + out$gamma^2*out$sig2x, nrow=nout, byrow=TRUE)
  if(model!="Gaussian"){
    out$theta <- matrix(C.out$theta.out, nrow=nout, byrow=TRUE)
    out$beta0 <- matrix(C.out$beta0.out, nrow=nout, byrow=TRUE)
  }
  if(ncov>0) out$eta <- matrix(C.out$eta.out, nrow=nout, byrow=TRUE)
  if(model=="Negative Binomial") out$nb_r = matrix(C.out$r.out, nrow=nout, byrow=TRUE)




  # Produce beta as a function of eigen-value
  Dseq <- seq(min(evals), max(evals), length=1000)
  c.beta <- matrix(NA, nrow=nout, ncol=length(Dseq))
  for(t in 1:nout){
    c.beta[t,] <- out$beta[t] + out$gamma[t]*sqrt((1-out$lamx[t]+out$lamx[t]*Dseq)/(1-out$lamz[t]+out$lamz[t]*Dseq))
  }


  if(model=="Gaussian"){
    beta.mn <-  matrix(apply(c.beta,2,function(x) mean(x)),nrow=1000,byrow=TRUE)
    beta.q025 <-  matrix(apply(c.beta,2,function(x) emp.hpd(x))[1,],nrow=1000,byrow=TRUE)
    beta.q975 <-  matrix(apply(c.beta,2,function(x) emp.hpd(x))[2,],nrow=1000,byrow=TRUE)
  }
  if(model!="Gaussian"){
    beta.mn <-  matrix(apply(exp(c.beta),2,function(x) mean(x)),nrow=1000,byrow=TRUE)
    beta.q025 <-  matrix(apply(exp(c.beta),2,function(x) emp.hpd(x))[1,],nrow=1000,byrow=TRUE)
    beta.q975 <-  matrix(apply(exp(c.beta),2,function(x) emp.hpd(x))[2,],nrow=1000,byrow=TRUE)
  }
  omega <- Dseq
  result <- eCAR.out(data_model = model,
                  beta_omega = cbind(beta.mn, beta.q025, beta.q975, omega),
                  posterior_draws = out,
                  DIC=NULL,
                  regrcoef = NULL)

  return(result)
}

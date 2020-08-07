# E: offset (e.g. for covid data E=scale(pop,center=F))
# y: data vector (possibly including NA) of length n
# x: exposure data vector of length n
# W: matrix n x n with 1 if i~j otherwise 0
# C: matrix with confounders, nrow(conf)=n
# TO DO: include random effects? doable using the inla formula extension...
# recall this even if not strictly needed; conf = tapply(C,rep(1:ncol(C),each=nrow(C)),function(i) i)

semipar.eCAR.Leroux <- function(y, x, W, E, C=NA,
                               model="Gaussian",
                               L=20, pcprior.sd=c(0.1,1), s2=10,
                               eval.fineGrid=FALSE,
                               verbose=FALSE, ...){
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package \"INLA\" is needed for this function to work. To install it, please go to http://www.r-inla.org/download.",
         call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package \"Matrix\" is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package \"splines\" is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # compute spectral decomposition of the structure matrix 'R'
  n <- length(y)
  M <- rowSums(W)
  R <- diag(M)-W
  Eigdec <- eigen(R)
  G <- Eigdec$vec


  # first, 3 utilities functions:
  # bsplines, rescale.row, myrgeneric.eigenLEROUX

  # Function to create the b-spline basis (P-spline method, Eilers 1996)
  bspline <- function(x, xl=min(x), xr=max(x), ndx, bdeg) {
    dx <- (xr - xl) / ndx
    knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
    B <- splines::spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
    B
  }

  rescale.row <- function(A,vec){
    t(apply(A, 1, function(x, vec) x*vec, vec=vec))
  }

  # rgenric model for the Gaussian case (eq 25, sec 4.3 overleaf paper)
  eigenLEROUX <- function (cmd = c("graph",
                                   "Q",
                                   "mu",
                                   "initial",
                                   "log.norm.const",
                                   "log.prior",
                                   "quit"),
                           theta = NULL) {
    interpret.theta <- function() {
      return(list(prec = exp(theta[1L]),
                  r = exp(theta[2L])/(1+exp(theta[2L])),
                  lam = exp(theta[3L])/(1+exp(theta[3L]))))
    }
    graph <- function() {
      return(Q())
    }
    Q <- function() {
      prec <- interpret.theta()$prec
      r <- interpret.theta()$r
      lam <- interpret.theta()$lam
      Q <- INLA::inla.as.sparse(diag((r/(prec*(1 - lam + lam*Eigdec$val)) + (1-r)/prec )^(-1)))
      return(Q)
    }
    mu <- function() {
      return(numeric(0))
    }
    log.norm.const <- function() {
      prec <- interpret.theta()$prec
      r <- interpret.theta()$r
      lam <- interpret.theta()$lam
      a <- (r/(prec*(1 - lam + lam*Eigdec$val)) + (1 - r) / prec)^(-1)
      val <- -0.5*n * log(2*pi) + 0.5*sum(log(a))
      return(val)
    }
    log.prior <- function() {
      prec <- interpret.theta()$prec
      val <- INLA::inla.pc.dprec(prec, pcprior.sd[2]/0.31, 0.01, log=T) + theta[1L] +
        theta[2L] - 2*log(exp(theta[2L]) + 1) +  # unif
        theta[3L] - 2*log(exp(theta[3L]) + 1)    # unif
      return(val)
    }
    initial <- function() {
      return(c(0,1,1))
    }
    quit <- function() {
      return(invisible())
    }
    if (is.null(theta)) theta <- initial()
    val <- do.call(match.arg(cmd), args = list())
    return(val)
  }


  # compute cubic b-spline basis (Eilers basis):
  B <- bspline(x=Eigdec$val, ndx=L-3, bdeg=3)
  if (eval.fineGrid) {
    B.pred <- bspline(x=seq(min(Eigdec$val),max(Eigdec$val),length.out=1000), ndx=L-3, bdeg=3)
  } else {
    B.pred<-B
  }

  if (model=="Gaussian") {
    xstar <- as.vector(t(G)%*%x)
    ystar <- as.vector(t(G)%*%y)
    BXstar <- sweep(B,1, xstar,"*")

    # semipar model fit (INLA) for Gaussian case
    if (is.null(C)) {
      stk <- INLA::inla.stack(data=list(y=ystar),
                              A=list(matrix(1,n,1),
                                     BXstar,
                                     diag(n)),
                              effects=list(intercept=1,
                                           id.beta = 1:L,
                                           id.z = 1:n))
      r <- INLA::inla(y ~  -1 + intercept +
                 f(id.beta,
                   model = "rw1",
                   scale.model = TRUE,
                   constr=FALSE,
                   hyper = list(prec=list(
                     prior="pc.prec",
                     param=c(pcprior.sd[1]/0.31, 0.01))))+
                 f(id.z, model=INLA::inla.rgeneric.define(
                   model=eigenLEROUX)),
               data = INLA::inla.stack.data(stk),
               control.family = list(hyper=list(prec=list(initial=12, fixed=TRUE))),
               control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
               control.compute = list(config=TRUE, dic=TRUE),
               verbose=verbose)
    } else{
      stk <- INLA::inla.stack(data=list(y=ystar),
                              A=list(matrix(1,n,1),
                                     BXstar,
                                     diag(n)),
                              effects=list(conf= cbind(rep(1,n), C),
                                           id.beta = 1:L,
                                           id.z = 1:n))
      r <- INLA::inla(y ~  -1 + conf +
                 f(id.beta,
                   model = "rw1",
                   constr=FALSE,
                   scale.model = TRUE,
                   hyper = list(prec=list(
                     prior="pc.prec",
                     param=c(pcprior.sd[1]/0.31, 0.01))))+
                   f(id.z, model=INLA::inla.rgeneric.define(
                     model=eigenLEROUX)),
                 data = INLA::inla.stack.data(stk),
                 control.family = list(hyper=list(prec=list(initial=12, fixed=TRUE))),
                 control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
                 control.compute = list(config=TRUE, dic=TRUE),
                 verbose=verbose)
    }

  } else {

    # non Gaussian cases (Binomial, Negative Binomial, Poisson)
    Z.tilde <- matrix(NA, nrow = n, ncol = L)
    for(i in 1:L){
      Z.tilde[,i] <-  rescale.row(G, B[,i]) %*% (t(G) %*% x)
    }

    if (is.null(C)) {
      stk <- INLA::inla.stack(data=list(y=y, E = E),
                              A=list(1, Z.tilde, diag(n)),
                              effects=list(intercept=rep(1,n),
                                           id.beta=1:L,
                                           id.z = 1:n))
      formula <- y ~ -1 + intercept +
        f(id.beta,
          model = "rw1",
          constr=FALSE, scale.model = TRUE,
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[1]/0.31, 0.01)))) +
        f(id.z, model = "generic1",
          Cmatrix = Matrix::Diagonal(x=1, n=n)-INLA::inla.as.sparse(R),
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[2]/0.31, 0.01))))
    } else {
      stk <- INLA::inla.stack(data=list(y=y, E = E),
                              A=list(1, Z.tilde, diag(n)),
                              effects=list(conf= cbind(rep(1,n), C),
                                           id.beta=1:L,
                                           id.z = 1:n))
      formula <- y ~ -1 + conf +
        f(id.beta,
          model = "rw1",
          constr=FALSE, scale.model = TRUE,
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[1]/0.31, 0.01)))) +
        f(id.z, model = "generic1",
          Cmatrix = Matrix::Diagonal(x=1, n=n)-INLA::inla.as.sparse(R),
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[2]/0.31, 0.01))))
    }

    #  likelihood Negative Binomial
    if (model=="Negative Binomial"){
      r <- INLA::inla(formula,
               data = INLA::inla.stack.data(stk),
               family='nbinomial', E=E,
               control.family = list(
                 variant=0,
                 hyper=list(theta= list(
                   prior="normal", param=c(0,1/s2)))),
               control.predictor = list(A = INLA::inla.stack.A(stk),
                                        compute=TRUE),
               control.compute = list(dic=TRUE, config=TRUE),
               verbose=verbose)
    }


    #  likelihood Binomial
    if (model=="Binomial"){
      r <- INLA::inla(formula,
               data = INLA::inla.stack.data(stk),
               family='binomial', Ntrials=E,
               control.predictor = list(A = INLA::inla.stack.A(stk),
                                        compute=TRUE),
               control.compute = list(dic=TRUE, config=TRUE),
               verbose=verbose)
    }

    # likelihood 'Poisson'
    if (model=="Poisson"){
      r <- INLA::inla(formula,
               data = INLA::inla.stack.data(stk),
               family='poisson', E=E,
               control.predictor = list(A = INLA::inla.stack.A(stk),
                                        compute=TRUE),
               control.compute = list(dic=TRUE, config=TRUE),
               verbose=verbose)
    }

  }

  # sample from the posterior and compute posterior mean curve
  sample.tmp <- INLA::inla.posterior.sample(1000,r)
  splinecoefs <- matrix(nrow=L, ncol=1000)
  prec.beta <- matrix(nrow=1, ncol=1000)
  prec.z <- matrix(nrow=1, ncol=1000)
  lambda.z <- matrix(nrow=1, ncol=1000)
  ind.beta <- c(paste("id.beta:",1:L, sep=""))
  for(j in 1:1000){
    splinecoefs[,j] <- (sample.tmp[[j]]$latent)[ind.beta,]
    prec.beta[1,j] <- sample.tmp[[j]]$hyperpar["Precision for id.beta"]
    prec.z[1,j] <- sample.tmp[[j]]$hyperpar["Precision for id.z"]
    lambda.z[1,j] <- sample.tmp[[j]]$hyperpar["Beta for id.z"]
  }

  # ARRANGE OUTPUT, MAKE THE TWO CASES: exp(beta_omega); beta_omega
  # default for non-Gaussian case: exp(beta_omega)

  if (model=="Gaussian") {
    beta.mn <- apply(B.pred%*%splinecoefs, 1, base::mean)
    beta.q025 <- apply(B.pred%*%splinecoefs, 1, stats::quantile, 0.025)
    beta.q975 <- apply(B.pred%*%splinecoefs, 1, stats::quantile, 0.975)
  } else {
    beta.mn <- apply(exp(B.pred%*%splinecoefs), 1, base::mean)
    beta.q025 <- apply(exp(B.pred%*%splinecoefs), 1, stats::quantile, 0.025)
    beta.q975 <- apply(exp(B.pred%*%splinecoefs), 1, stats::quantile, 0.975)
  }

  # return results
  if (!eval.fineGrid) {
    out <- eCAR.out(data_model = model,
                    beta_omega = cbind(beta.mn, beta.q025, beta.q975,Eigdec$val),
                    posterior_draws = list(postsample.beta=splinecoefs,
                                       postsample.prec.beta=as.numeric(prec.beta),
                                       postsample.prec.z=as.numeric(prec.z),
                                       postsample.lambda.z=as.numeric(lambda.z),
                                       B.pred = B.pred))
  } else {
    out <- eCAR.out(data_model = model,
                    beta_omega = cbind(beta.mn,
                                   beta.q025,
                                   beta.q975,
                                   seq(min(Eigdec$val),max(Eigdec$val),length.out=1000)),
                    posterior_draws = list(postsample.beta=splinecoefs,
                                       postsample.prec.beta=as.numeric(prec.beta),
                                       postsample.prec.z=as.numeric(prec.z),
                                       postsample.lambda.z=as.numeric(lambda.z),
                                       B.pred = B.pred))
  }

  return(out)
}

# TO DO list:
# - produce a list as output including the post.sample;
# - add 'standard=TRUE' ? may be helpful to provide also the constant beta in the output
# - #CHECK THIS!!!! Do I need to do matrix(apply(G,2,sum),n,1) ? it should not matter



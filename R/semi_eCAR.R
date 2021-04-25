# E: offset (e.g. for covid data E=scale(pop,center=F))
# y: data vector (possibly including NA) of length n
# x: exposure data vector of length n
# W: matrix n x n with 1 if i~j otherwise 0
# C: matrix with confounders, nrow(conf)=n; NOT CURRENTLY WORKING PROPERLY...
# TO DO: include random effects? doable using the inla formula extension...
# recall this even if not strictly needed; conf = tapply(C,rep(1:ncol(C),each=nrow(C)),function(i) i)

semipar.eCAR.Leroux <- function(y, x, W, E, C=NULL,
                                names.covariates=NULL,
                                model="Gaussian",
                                L=10, pcprior.sd=c(0.1,1), s2=10,
                                method = "spectral",
                                num.threads.inla = NULL,
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
  if (!is.null(num.threads.inla)) INLA::inla.setOption(num.threads = num.threads.inla)

  n <- length(y)
  M <- rowSums(W)
  R <- diag(M)-W
  # compute spectral decomposition of the structure matrix 'R' (we do it also for the 'naive')
  Eigdec <- eigen(R)
  G <- Eigdec$vec
  omega <- Eigdec$val

  if (method=="spectral"){

    # define 3 utilities functions: bsplines(), rescale.row(), myrgeneric.eigenLEROUX()
    # may this be doine outside semipar.eCAR.Leroux?

    # create the b-spline basis (P-spline method, Eilers 1996)
    bspline <- function(x, xl=min(x), xr=max(x), ndx, bdeg) {
      dx <- (xr - xl) / ndx
      knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
      B <- splines::spline.des(knots, x, bdeg + 1, 0 * x, outer.ok = TRUE)$design
      B
    }

    rescale.row <- function(A,vec){
      t(apply(A, 1, function(x, vec) x*vec, vec=vec))
    }

    # 'rgeneric' model for the Gaussian case (eq 25, sec 4.3 overleaf paper)
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
        Q <- INLA::inla.as.sparse(diag((r/(prec*(1 - lam + lam*omega)) + (1-r)/prec )^(-1)))
        return(Q)
      }
      mu <- function() {
        return(numeric(0))
      }
      log.norm.const <- function() {
        prec <- interpret.theta()$prec
        r <- interpret.theta()$r
        lam <- interpret.theta()$lam
        a <- (r/(prec*(1 - lam + lam*omega)) + (1 - r) / prec)^(-1)
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

    # deal with the covariate matrix 'C'
    if (is.null(C)) {
      X.cov <- data.frame(intercept=rep(1,n))
      form.cov <- paste(" -1 + ", " intercept ")
      form1 <- stats::as.formula(paste("y", form.cov, sep=" ~" ))
      # NOTE: if we want to pass a formula as an input we need to do work here in 'form1'
    } else {
      X.cov <- data.frame(intercept=1, C)
      colnames(X.cov) <- c("intercept", paste("cov",1:ncol(C), sep = ""))
      if (!is.null(names.covariates) & length(names.covariates)==ncol(C)) colnames(X.cov) <- c("intercept", names.covariates)
      form.cov <- paste(" -1 + ", paste(colnames(X.cov), collapse = " + "))
      form1 <- stats::as.formula(paste("y", form.cov, sep=" ~" ))
      # NOTE: if we want to pass a formula as an input we need to do work here in 'form1'
    }

    # compute cubic b-spline basis (Eilers basis):
    B <- bspline(x=omega, ndx=L-3, bdeg=3)
    B.pred <- bspline(x=seq(min(omega),max(omega),length.out=1000), ndx=L-3, bdeg=3)


    ### Gaussian case ###
    if (model=="Gaussian") {
      # xstar <- as.vector(t(G)%*%x)
      # ystar <- as.vector(t(G)%*%y)
      xstar <- as.vector(crossprod(G,x))
      ystar <- as.vector(crossprod(G[!is.na(y),],y))
      BXstar <- sweep(B,1, xstar,"*")

      # inla call for method=="spectral" ; Gaussian case
      stk <- INLA::inla.stack(data=list(y=ystar),
                              A=list(1,
                                     BXstar,
                                     diag(n)),
                              effects=list(X.cov,
                                           #BBmisc::convertColsToList(X.cov),
                                           id.beta = 1:L,
                                           id.z = 1:n))
      form2 <- y ~ . +
        f(id.beta,
          model = "rw1",
          scale.model = TRUE,
          constr=FALSE,
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(
              pcprior.sd[1]/0.31,
              0.01)))) +
        f(id.z, model=INLA::inla.rgeneric.define(
          model=eigenLEROUX,
          omega=Eigdec$val,
          pcprior.sd=pcprior.sd,
          n=n))
      res <- INLA::inla(stats::update.formula(form1,form2),
                data = INLA::inla.stack.data(stk),
                control.family = list(hyper=list(prec=list(initial=12, fixed=TRUE))),
                control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
                control.compute = list(config=TRUE, dic=TRUE),
                verbose=verbose, ...)

    } else {

      ### non Gaussian cases (Binomial, Negative Binomial, Poisson) ###

      Z.tilde <- matrix(NA, nrow = n, ncol = L)
      for(i in 1:L){
        Z.tilde[,i] <-  rescale.row(G, B[,i]) %*% (t(G) %*% x)
      }

      # inla call for method=="spectral" ; non Gaussian cases
      stk <- INLA::inla.stack(data=list(y=y, E = E),
                              A=list(1, Z.tilde, diag(n)),
                              effects=list(X.cov,
                                           #BBmisc::convertColsToList(X.cov),
                                           id.beta=1:L,
                                           id.z = 1:n))
      form2 <- y ~ . +
        f(id.beta,
          model = "rw1",
          scale.model = TRUE,
          constr=FALSE,
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[1]/0.31, 0.01)))) +
        f(id.z, model = "generic1",
          Cmatrix = Matrix::Diagonal(x=1, n=n)-INLA::inla.as.sparse(R),
          hyper = list(prec=list(
            prior="pc.prec",
            param=c(pcprior.sd[2]/0.31, 0.01))))

      #  likelihood Negative Binomial
      if (model=="Negative Binomial"){
        res <- INLA::inla(stats::update.formula(form1,form2),
                  data = INLA::inla.stack.data(stk),
                  family='nbinomial', E=E,
                  control.family = list(
                    variant=0,
                    hyper=list(theta= list(
                      prior="normal", param=c(0,1/s2)))),
                  control.predictor = list(A = INLA::inla.stack.A(stk),
                                           compute=TRUE),
                  control.compute = list(dic=TRUE, config=TRUE),
                  verbose=verbose, ...)
      }

      #  likelihood Binomial
      if (model=="Binomial"){
        res <- INLA::inla(stats::update.formula(form1,form2),
                  data = INLA::inla.stack.data(stk),
                  family='binomial', Ntrials=E,
                  control.predictor = list(A = INLA::inla.stack.A(stk),
                                           compute=TRUE),
                  control.compute = list(dic=TRUE, config=TRUE),
                  verbose=verbose, ...)
      }

      # likelihood 'Poisson'
      if (model=="Poisson"){
        res <- INLA::inla(stats::update.formula(form1,form2),
                  data = INLA::inla.stack.data(stk),
                  family='poisson', E=E,
                  control.predictor = list(A = INLA::inla.stack.A(stk),
                                           compute=TRUE),
                  control.compute = list(dic=TRUE, config=TRUE),
                  verbose=verbose, ...)
      }
    }

    # sample from the posterior and compute beta_omega
    L.supported = sum(colSums(B) != 0)
    sample.tmp <- INLA::inla.posterior.sample(1000,res)
    splinecoefs <- matrix(nrow=L.supported, ncol=1000)
    prec.beta <- matrix(nrow=1, ncol=1000)
    prec.z <- matrix(nrow=1, ncol=1000)
    lambda.z <- matrix(nrow=1, ncol=1000)
    ind.beta <- c(paste("id.beta:",1:L.supported, sep=""))
    for(j in 1:1000){
      splinecoefs[,j] <- (sample.tmp[[j]]$latent)[ind.beta,]
      prec.beta[1,j] <- sample.tmp[[j]]$hyperpar["Precision for id.beta"]
      prec.z[1,j] <- sample.tmp[[j]]$hyperpar["Precision for id.z"]
      lambda.z[1,j] <- sample.tmp[[j]]$hyperpar["Beta for id.z"]
    }


    # ARRANGE OUTPUT, MAKE THE TWO CASES: exp(beta_omega); beta_omega
    # default for non-Gaussian case: exp(beta_omega)
    if (model=="Gaussian") {
      beta.mn <- apply(B.pred[,which(colSums(B) != 0)]%*%splinecoefs, 1, base::mean)
      beta.q025 <- apply(B.pred[,which(colSums(B) != 0)]%*%splinecoefs, 1, stats::quantile, 0.025)
      beta.q975 <- apply(B.pred[,which(colSums(B) != 0)]%*%splinecoefs, 1, stats::quantile, 0.975)
    } else {
      beta.mn <- apply(exp(B.pred[,which(colSums(B) != 0)]%*%splinecoefs), 1, base::mean)
      beta.q025 <- apply(exp(B.pred[,which(colSums(B) != 0)]%*%splinecoefs), 1, stats::quantile, 0.025)
      beta.q975 <- apply(exp(B.pred[,which(colSums(B) != 0)]%*%splinecoefs), 1, stats::quantile, 0.975)
    }

    # # regr coef estimates
    # if (model=="Gaussian") {
    #   fixed <- res$summary.fixed
    #   colnames(fixed) <- c("mean", "sd", "quant0.025", "quant0.25", "quant0.5", "quant0.75", "quant0.975")
    #   rownames(fixed) <- colnames(X.cov)
    # } else {
    #   fixed.list <- lapply(res$marginals.fixed, FUN=function(x) INLA::inla.zmarginal(INLA::inla.tmarginal(function(a) exp(a), x), silent=TRUE))
    #   fixed <- matrix(unlist(lapply(fixed.list, data.frame)), ncol(X.cov), 7, byrow=TRUE)
    #   colnames(fixed) <- c("mean", "sd", "quant0.025", "quant0.25", "quant0.5", "quant0.75", "quant0.975")
    #   rownames(fixed) <- c("intercept", paste("exp", names(fixed.list)[-1], sep="_"))
    #   fixed["intercept",] <- unlist(res$summary.fixed["intercept",])
    # }
    # regr coef est
    if (model=="Gaussian") {
      fixed <- res$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")]
      rownames(fixed) <- colnames(X.cov)
    } else {
      fixed.list <- lapply(res$marginals.fixed, FUN=function(x) INLA::inla.zmarginal(INLA::inla.tmarginal(function(a) exp(a), x), silent=TRUE))
      fixed.tmp <- matrix(unlist(lapply(fixed.list, data.frame)), nrow=length(fixed.list), ncol=7, byrow=TRUE)
      colnames(fixed.tmp) <- c("mean_exp", "sd", "0.025quant", "0.25quant", "0.5quant", "0.75quant", "0.975quant")
      rownames(fixed.tmp) <- colnames(X.cov)#c("intercept", paste("exp", names(fixed.list)[-1], sep="_"))
      fixed.tmp["intercept",] <- unlist(res$summary.fixed["intercept",])
      fixed <- fixed.tmp[, c("mean_exp", "sd", "0.025quant", "0.975quant")]
      #fixed
    }

    # produce 'out'
    out <- eCAR.out(data_model = model,
                    beta_omega = data.frame(beta.mn = beta.mn,
                                            beta.q025 = beta.q025,
                                            beta.q975 = beta.q975,
                                            omega = seq(min(omega),max(omega),length.out=1000)),
                    posterior_draws = list(postsample.beta=splinecoefs,
                                           postsample.prec.beta=as.numeric(prec.beta),
                                           postsample.prec.z=as.numeric(prec.z),
                                           postsample.lambda.z=as.numeric(lambda.z),
                                           B.pred = B.pred),
                    DIC=res$dic$dic,
                    regrcoef = as.data.frame(fixed))
  } else {
    # method == "naive"

    # deal with the covariate matrix 'C'
    if (is.null(C)) {
      X.cov <- data.frame(intercept=1, x=x)
      form.cov <- paste(" -1 + ", " intercept + ", " x ")
      form1 <- stats::as.formula(paste("y", form.cov, sep=" ~" ))
      # NOTE: if we want to pass a formula as an input we need to do add it here in 'form1'
    } else {
      X.cov <- data.frame(intercept=1, x=x)
      X.cov <- cbind(X.cov, C)
      colnames(X.cov) <- c("intercept", "x", paste("cov",1:ncol(C), sep = ""))
      if (!is.null(names.covariates) & length(names.covariates)==ncol(C)) colnames(X.cov) <- c("intercept", "x", names.covariates)
      form.cov <- paste(" -1 + ", paste(colnames(X.cov), collapse = " + "))
      form1 <- stats::as.formula(paste("y", form.cov, sep=" ~" ))
      # NOTE: if we want to pass a formula as an input we need to add it here in 'form1'
    }

    if (model=="Gaussian") {
      stk = INLA::inla.stack(
        data=list(y=y),
        A=list(1, diag(n)),
        effects=list(X.cov,
                     id.z=1:n))

      form2 <- y ~ . +  f(id.z, model = "generic1",
                          Cmatrix = Matrix::Diagonal(x=1, n=n) - INLA::inla.as.sparse(R),
                          hyper = list(prec=list(
                            prior="pc.prec",
                            param=c(pcprior.sd[2]/0.31, 0.01))))

      res <- INLA::inla(stats::update.formula(form1, form2),
                        data = INLA::inla.stack.data(stk),
                        control.family = list(hyper=list(prec=list(initial=12, fixed=TRUE))),
                        control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
                        control.compute = list(dic=TRUE, config=TRUE),
                        verbose=TRUE, ...)
    } else {
      stk = INLA::inla.stack(
        data=list(y=y, E=E),
        A=list(1, diag(n)),
        effects=list(X.cov,
                     id.z=1:n))

      form2 <- y ~ . +  f(id.z, model = "generic1",
                          Cmatrix = Matrix::Diagonal(x=1, n=n) - INLA::inla.as.sparse(R),
                          hyper = list(prec=list(
                            prior="pc.prec",
                            param=c(pcprior.sd[2]/0.31, 0.01))))


      if (model=="Poisson"){
        res <- INLA::inla(stats::update.formula(form1,form2),
                          data = INLA::inla.stack.data(stk),
                          family='poisson', E=E,
                          control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
                          control.compute = list(dic=TRUE, config=TRUE),
                          verbose=verbose, ...)
      }
      if (model=="Binomial"){
        res <- INLA::inla(stats::update.formula(form1,form2),
                          data = INLA::inla.stack.data(stk),
                          family='binomial', Ntrials=E,
                          control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
                          control.compute = list(dic=TRUE, config=TRUE),
                          verbose=verbose, ...)
      }
      if (model=="Negative Binomial"){
        res <- INLA::inla(stats::update.formula(form1, form2),
                          data = INLA::inla.stack.data(stk),
                          family='nbinomial', E=E,
                          control.family = list(variant=0,
                                                hyper=list(
                                                  theta=list(
                                                    prior="normal", param=c(0,1/s2)))),
                          control.predictor = list(A = INLA::inla.stack.A(stk), compute=TRUE),
                          control.compute = list(dic=TRUE, config=TRUE),
                          verbose=verbose, ...)
      }
    }

    # regr coef est
    if (model=="Gaussian") {
      fixed <- res$summary.fixed[,c("mean", "sd", "0.025quant", "0.975quant")]
      rownames(fixed) <- colnames(X.cov)
    } else {
      fixed.list <- lapply(res$marginals.fixed, FUN=function(x) INLA::inla.zmarginal(INLA::inla.tmarginal(function(a) exp(a), x), silent=TRUE))
      fixed.tmp <- matrix(unlist(lapply(fixed.list, data.frame)), nrow=ncol(X.cov), ncol=7, byrow=TRUE)
      colnames(fixed.tmp) <- c("mean_exp(beta)", "sd", "0.025quant", "0.25quant", "0.5quant", "0.75quant", "0.975quant")
      rownames(fixed.tmp) <- colnames(X.cov)#c("intercept", paste("exp", names(fixed.list)[-1], sep="_"))
      fixed.tmp["intercept",] <- unlist(res$summary.fixed["intercept",])
      fixed <- fixed.tmp[, c("mean_exp(beta)", "sd", "0.025quant", "0.975quant")]
      #fixed
    }

    # produce 'out'
    out <- eCAR.out(data_model = model,
                    beta_omega = data.frame(beta.mn = rep(fixed["x",1], 1000),
                                            beta.q025 = rep(fixed["x",3], 1000),
                                            beta.q975 = rep(fixed["x",4], 1000),
                                            omega = seq(min(omega),max(omega),length.out=1000)),
                    posterior_draws = list(postsample.beta=NULL,
                                           postsample.prec.beta=NULL,
                                           postsample.prec.z=NULL,
                                           postsample.lambda.z=NULL,
                                           B.pred = NULL),
                    DIC=res$dic$dic,
                    regrcoef = as.data.frame(fixed))
  }
  return(out)
}



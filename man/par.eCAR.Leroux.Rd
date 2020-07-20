\name{par.eCAR.Leroux}
\alias{par.eCAR.Leroux}
\title{R wrapper that accesses C code to fit parametric Leroux CAR spatial regression model that accounts for spatial confounding}
\usage{
par.eCAR.Leroux(y, x, W,
                E=NULL,
                C=NULL,
                model="Gaussian",
                m=0, s2=10,
                alamx=1, blamx=1,
                alamz=1, blamz=1,
                asig=1, bsig=1,
                atau=1, btau=1,
                asigx=1, bsigx=1,
                mb0=0,s2b0=100,
                me=0,s2e=100,
                mx=0,s2x=100,
                draws=10000, burn=5000, thin=5,
                verbose=TRUE)
}
\arguments{
\item{y}{response vector}

\item{x}{covariate vector for which casual effect is desired}

\item{W}{neighborhood matrix comprised of zeros and ones}

\item{E}{This argument is ignored if model is Gaussian.  For other models it takes on the following:

    * Poisson - E is vector that contains expected counts

    * Binomial - E is vector that contains number of trials

    * Negative Binomial - E is vector that contains an offset.
}

\item{C}{design matrix for the covariates that are included as controls}

\item{model}{Specifies the likelihood or data model.  Options are "Gaussian", "Poisson", "Binomial", "Negative Binomial"}

\item{m}{prior mean for (beta, alpha), where alpha = rho(sigz/sigx).  default is 0.}

\item{s2}{prior variance for (beta, alpha),where alpha = rho(sigz/sigx). default is 10}

\item{alamx}{prior shape1 parameter for lam.x, default is 1}

\item{blamx}{prior shape2 parameter for lam.x, default is 1}

\item{alamz}{prior shape1 parameter for lam.z, default is 1}

\item{blamz}{prior shape2 parameter for lam.z, default is 1}

\item{asig}{prior shape parameter for sigma2, default is 1}

\item{bsig}{prior scale parameter for sigma2, default is 1}

\item{atau}{prior shape parameter for tau, where tau = sigma2.z(1-rho^2). default is 1}

\item{btau}{prior scale parameter for tau, where tau = sigma2.z(1-rho^2).  default is 1}

\item{asigx}{prior shape parameter for sigma2.x, default is 1}

\item{bsixg}{prior scale parameter for sigma2.x, default is 1}

\item{mb0}{prior mean parameter for beta0, default is 0}

\item{s2b0}{prior variance parameter for beta0, default is 100}

\item{me}{prior mean parameter for eta, default is 0}

\item{s2e}{prior variance parameter for eta, default is 100}

\item{mx}{prior mean parameter for xi, default is 0. Only used for negative binomial model}

\item{s2x}{prior variance parameter for eta, default is 100. Only used for negative binomial model}

\item{draws}{number of MCMC iterates to be collected. default is 10000}

\item{burn}{number of MCMC iterates discared as burn-in. default is 5000}

\item{thin}{number by which the MCMC chain is thinne. default is 5}

\item{verbose}{If TRUE, then details associated with data being fit are printed to screen along with MCMC iterate counter}
}
\value{
A list containing arrays filled with MCMC iterates corresponding to model parameters. Letting T = (draws - burn)/thin be the number of MCMC iterates requested, the output list contains the following

\item{beta}{a T-dimensional array containing MCMC iterates associated with beta}

\item{alpha}{a T-dimensional array containing MCMC iterates associated with alpha}

\item{tau}{a T-dimensional array containing MCMC iterates associated with tau}

\item{sig2x}{a T-dimensional array containing MCMC iterates associated with sigma2.x}

\item{lamx}{a T-dimensional array containing MCMC iterates associated with lambda.x}

\item{lamz}{a T-dimensional array containing MCMC iterates associated with lambda.z}

\item{sig2}{a T-dimensional array containing MCMC iterates associated with sigma2}

\item{rho}{a T-dimensional array containing MCMC iterates associated with rho}

\item{sig2z}{a T-dimensional array containing MCMC iterates associated with sigma2.z}

\item{c.beta}{a T-dimensional array containing MCMC iterates associated with bias corrected beta}
}
\description{
\code{parametric.eCAR.Leroux} is the main function used to fit the parametric Leroux CAR model specified in the spectral domain.
}

\references{
}
\examples{

  # Our R-package
  library(eCAR)
  library(spdep) # to create a neigbhorhood matrix



  n1 <- n2 <- 20
  n <- n1*n2
  W <- nb2mat(cell2nb(n1, n2, type = "queen"), zero.policy=TRUE, style="B")
  M <- diag(apply(W,1,sum))
  R <- M-W
  e.dec <- eigen(R)
  e.val <- e.dec$values
  D.eigval = diag(e.val)

  # Generate data from the model specified in the spectral domain
  beta.true = 1.1
  sig2.true <- 0.1
  sigma.x = 1
  sigma.z = 1.5
  sigma2.x = sigma.x^2
  sigma2.z = sigma.z^2
  lam.x <- lam.x.true <- 0.1;
  lam.z <- lam.z.true <-  0.9
  rho <- rho.true <- 0.7
  alpha.true <- rho.true*(sigma.z/sigma.x)
  tau.true <- sigma2.z*(1-rho.true^2)


  # Function that computes matrix square root or its inverse
  matrix.sqrt <- function(X, inv=FALSE){
	  e.dec <- eigen(X, symmetric=TRUE)
	  if(inv){
		  out <- e.dec$vectors %*% diag(1/sqrt(e.dec$values)) %*% solve(e.dec$vectors)
	  }
	  if(!inv){
		  out <- e.dec$vectors %*% diag(sqrt(e.dec$values)) %*% solve(e.dec$vectors)
	  }
	  out
  }

  R_zz <- solve((1-lam.z.true)*diag(n) + lam.z.true*D.eigval)
  R_xx <- solve((1-lam.x.true)*diag(n) + lam.x.true*D.eigval)

  Rzz.sqrt = matrix.sqrt((1-lam.z.true)*diag(n) + lam.z.true*D.eigval, inv=1)
  Rxx.sqrt = matrix.sqrt((1-lam.x.true)*diag(n) + lam.x.true*D.eigval, inv=1)

  R = rbind(cbind(sigma2.x*R_xx,
                  rho*sqrt(sigma2.x*sigma2.z)*Rxx.sqrt%*%t(Rzz.sqrt)),
            cbind(rho*sqrt(sigma2.x*sigma2.z)*Rzz.sqrt%*%t(Rxx.sqrt),
                  sigma2.z*R_zz))

  # generate data
  set.seed(1)
  xz = rmvnorm(1,mean=rep(0,2*n),sigma = R)
  X.star = xz[1:n]
  Z.star = xz[(n+1):(2*n)]
  Y.star = beta.true * X.star + Z.star + rnorm(n, 0, sqrt(sig2.true))

  # Fit the model (see help file for default prior values)
  fit <- eCAR.Leroux(y=Y.star, x=X.star, W=W, verbose=FALSE, draws=100000, burn=50000, thin=50)

  plot(fit$beta, type='l'); abline(h=beta.true, col='red')
  plot(fit$c.beta, type='l'); abline(h=beta.true, col='red')
  sapply(fit, summary)


}

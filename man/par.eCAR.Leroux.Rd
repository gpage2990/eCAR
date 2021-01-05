\name{par.eCAR.Leroux}
\alias{par.eCAR.Leroux}
\title{R wrapper that accesses C code to fit parametric Leroux CAR spatial regression model that accounts for spatial confounding}
\description{
\code{par.eCAR.Leroux} is the main function used to fit the parametric Leroux CAR model specified in the spectral domain.
}

\usage{
par.eCAR.Leroux(y, x, W,
                E=NULL,
                C=NULL,
                model="Gaussian",
                joint_prior_lamx_lamz = FALSE,
                lamx.fix.val = NULL,
                sig2x.fix.val = NULL,
                m=0, s2=10,
                alamx=1, blamx=1,
                alamz=1, blamz=1,
                asig=1, bsig=1,
                atau=1, btau=1,
                asigx=1, bsigx=1,
                mb0=0,s2b0=100,
                me=0,s2e=100,
                mx=0,s2x=100,
                tau_cand_sd = 1,
                sig2_cand_sd = 1,
                draws=10000, burn=5000, thin=5,
                verbose=TRUE)
}
\arguments{
\item{y}{response vector}

\item{x}{covariate vector for which casual effect is desired}

\item{W}{neighborhood matrix comprised of zeros and ones}

\item{E}{This argument is ignored if model is Gaussian.  For other models it takes on the following:
\itemize{
    \item Poisson - E is vector that contains expected counts

    \item Binomial - E is vector that contains number of trials

    \item Negative Binomial - E is vector that contains an offset.
}
}
\item{C}{design matrix for the covariates that are included as controls}

\item{model}{Specifies the likelihood or data model.  Options are "Gaussian", "Poisson", "Binomial", "Negative Binomial"}


\item{joint_prior_lamx_lamz}{Logical.  If TRUE, then a uniform prior on space such that lamz > lamx.  If FALSE, independent beta priors are used.  }


\item{lamx.fix.val}{If a value is supplied then lambda_x is not updated in the MCMC algorithm, but rather treated as the fixed known supplied value}

\item{sig2x.fix.val}{If a value is supplied then sigma2_x is not updated in the MCMC algorithm, but rather treated as the fixed known supplied value}

\item{m}{prior mean for (beta, gamma), where gamma = rho*(sigz/sigx).  default is 0.}

\item{s2}{prior variance for (beta, gamma),where gamma = rho*(sigz/sigx). default is 10}

\item{alamx}{prior shape1 parameter for lam.x, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{blamx}{prior shape2 parameter for lam.x, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{alamz}{prior shape1 parameter for lam.z, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{blamz}{prior shape2 parameter for lam.z, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{asig}{prior shape parameter for sigma2, default is 1. Only used if model is Gaussian}

\item{bsig}{prior scale parameter for sigma2, default is 1. Only used if model is Gaussian}

\item{atau}{prior shape parameter for tau, where tau = sigma2.z*(1-rho^2). default is 1}

\item{btau}{prior scale parameter for tau, where tau = sigma2.z*(1-rho^2).  default is 1}

\item{asigx}{prior shape parameter for sigma2.x, default is 1}

\item{bsigx}{prior scale parameter for sigma2.x, default is 1}

\item{mb0}{prior mean parameter for beta0, default is 0.  Only used if model is not Gaussian}

\item{s2b0}{prior variance parameter for beta0, default is 100. Only used if model is not Gaussian}

\item{me}{prior mean parameter for eta, default is 0. Only used if C is not NULL}

\item{s2e}{prior variance parameter for eta, default is 100. Only used if C is not NULL}

\item{mx}{prior mean parameter for xi, default is 0. Only used for negative binomial model}

\item{s2x}{prior variance parameter for eta, default is 100. Only used for negative binomial model}

\item{tau_cand_sd}{standard deviation for candidate density in Metropolis step for tau.  Default is 1}

\item{sig2_cand_sd}{standard deviation for candidate density in Metropolis step for sig2.  Default is 1. Only used if model is Gaussian}

\item{draws}{number of MCMC iterates to be collected. default is 10000}

\item{burn}{number of MCMC iterates discared as burn-in. default is 5000}

\item{thin}{number by which the MCMC chain is thinned. default is 5}

\item{verbose}{If TRUE, then details associated with data being fit are printed to screen along with MCMC iterate counter}
}


\value{
The function returns an \code{eCAR} object which is a list that contains the following
\item{data_model}{Character indicating which model was fit}
\item{beta_omega}{Matrix that contains respectively, the posterior mean lower and upper quantiles of the (spatial scale)-varying beta at each omega value (for the non Gaussian cases it is the exponentiated beta).}
\item{posterior_draws}{List containing posterior draws of the following parameters

\enumerate{

  \item{beta: vector containing draws from the posterior distribution of exposure coefficient }
  \item{gamma: vector containing draws from the posterior distribution of gamma=rho*(sigz/sigx)}
  \item{tau: vector containing draws from tau=sigma2.z*(1-rho^2)}
  \item{sig2x: vector containing draws from sig2x }
  \item{lamx: vector containing draws from the posterior distribution of lamx }
  \item{lamz: vector containing draws from the posterior distribution of lamz}
  \item{sigma2: vector containing draws from the posterior distribution of sigma2.  Only available if model is Gaussian}
  \item{rho: vector containing draws from the posterior distribution of rho}
  \item{sig2z: vector containing draws from the posterior distribution of sig2z}
  \item{theta: matrix containing draws from the posterior distribution of theta. Only available if model is not Gaussian}
  \item{beta0: vector containing draws from the posterior distribution of beta0. Only available if model is not Gaussian}
  \item{eta: matrix containing draws from the posterior distribution of eta, regression coefficients of additional confouner variables.  Only available if C is not NULL}
  \item{nb_r: matrix containing draws from the posterior distribution of nb_r.  0nly available if model is Negative Binomial}
}
}

\item{DIC}{Not available from parametric model yet}

\item{regrcoef}{Not available from parametric model yet}

}

\references{
Guan, Y; Page, G.L.; Reich, B.J.; Ventrucci, M.; Yang, S; "A spectral adjustment for spatial confounding" <arXiv:2012.11767>
}
\examples{

# Our R-package
library(eCAR)

data(lipcancer)


W <- lipcancer$neighborhood.Matrix
M <- diag(apply(W,1,sum))
R <- M-W
e.dec <- eigen(R)
e.val <- e.dec$values
D.eigval = diag(e.val)


Y <- lipcancer$data$observed
X <- lipcancer$data$pcaff
E <- lipcancer$data$expected

set.seed(101)
fit1 <- par.eCAR.Leroux(y=Y, x=X, W=W, E=E, C=NULL, model="Poisson",
                         draws=10000, burn=5000, thin=1, verbose=FALSE,
                         joint_prior_lamx_lamz=FALSE, m=0, s2=4)


plot(fit1)
}

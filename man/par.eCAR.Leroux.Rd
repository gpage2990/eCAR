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

    * Poisson - E is vector that contains expected counts

    * Binomial - E is vector that contains number of trials

    * Negative Binomial - E is vector that contains an offset.
}

\item{C}{design matrix for the covariates that are included as controls}

\item{model}{Specifies the likelihood or data model.  Options are "Gaussian", "Poisson", "Binomial", "Negative Binomial"}


\item{joint_prior_lamx_lamz}{Logical.  If TRUE, then a uniform prior on space such that lamz > lamx.  If FALSE, independent beta priors are used.  }


\item{m}{prior mean for (beta, gamma), where gamma = rho*(sigz/sigx).  default is 0.}

\item{s2}{prior variance for (beta, gamma),where gamma = rho*(sigz/sigx). default is 10}

\item{alamx}{prior shape1 parameter for lam.x, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{blamx}{prior shape2 parameter for lam.x, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{alamz}{prior shape1 parameter for lam.z, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{blamz}{prior shape2 parameter for lam.z, default is 1. Only used if joint_prior_lamx_lamz = FALSE}

\item{asig}{prior shape parameter for sigma2, default is 1}

\item{bsig}{prior scale parameter for sigma2, default is 1}

\item{atau}{prior shape parameter for tau, where tau = sigma2.z*(1-rho^2). default is 1}

\item{btau}{prior scale parameter for tau, where tau = sigma2.z*(1-rho^2).  default is 1}

\item{asigx}{prior shape parameter for sigma2.x, default is 1}

\item{bsigx}{prior scale parameter for sigma2.x, default is 1}

\item{mb0}{prior mean parameter for beta0, default is 0}

\item{s2b0}{prior variance parameter for beta0, default is 100}

\item{me}{prior mean parameter for eta, default is 0}

\item{s2e}{prior variance parameter for eta, default is 100}

\item{mx}{prior mean parameter for xi, default is 0. Only used for negative binomial model}

\item{s2x}{prior variance parameter for eta, default is 100. Only used for negative binomial model}

\item{tau_cand_sd}{standard deviation for candidate density in Metropolis step for tau.  Default is 1}

\item{sig2_cand_sd}{standard deviation for candidate density in Metropolis step for tau.  Default is 1}

\item{draws}{number of MCMC iterates to be collected. default is 10000}

\item{burn}{number of MCMC iterates discared as burn-in. default is 5000}

\item{thin}{number by which the MCMC chain is thinned. default is 5}

\item{verbose}{If TRUE, then details associated with data being fit are printed to screen along with MCMC iterate counter}
}


\value{
The function returns an \code{eCAR} object which is a list that contains the following
\enumerate{
  \item {data_model: Character indicating which model was fit}
	\item {beta_omega: Matrix comprised of four coulums.  The first contains the estimated scale varying estimate of beta, the second and third contain the 2.5\% and 97.5\% credible band, and the fourth a sequence of values corresponding to eigen-values}
	\item {posterior_draws: a list containing posterior draws for all parameters in the model}
 }
}

\references{
Will put something here.  Maybe arXiv submission?
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

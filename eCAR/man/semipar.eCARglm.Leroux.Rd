\name{semipar.eCARglm.Leroux}
\alias{semipar.eCARglm.Leroux}
\title{R wrapper that accesses INLA to fit semi-parametric Leroux CAR spatial regression model that accounts for spatial confounding.}
\description{
\code{semipar.eCAR.Leroux} is the main function used to fit the semi-parametric Leroux CAR model specified in the spectral domain.
}
\usage{
semipar.eCARglm.Leroux(y, x, W, E,
                       C=NA,
                       model="Poisson",
                       L=20,
                       pcprior.sd=c(0.1,1)
                       s2=2,
                       eval.fineGrid=FALSE,
                       verbose=TRUE)
}

\arguments{
\item{y}{response vector}

\item{x}{covariate vector for which casual effect is desired}

\item{W}{neighborhood matrix comprised of zeros and ones}

\item{E}{Offset value whose specification depends on the data model selected such that for

    * Poisson - E is vector that contains expected counts

    * Binomial - E is vector that contains number of trials

    * Negative Binomial - E is vector that contains an offset.
}

\item{C}{design matrix for the covariates that are included as controls}

\item{model}{Specifies the likelihood or data model.  Options are "Gaussian", "Poisson", "Binomial", "Negative Binomial"}

\item{m}{prior mean for (beta, alpha), default is 0.}

\item{s2}{prior variance for (beta, alpha), default is 10}

\item{alamx}{prior shape1 parameter for lam.x, default is 1}

\item{blamx}{prior shape2 parameter for lam.x, default is 1}

\item{alamz}{prior shape1 parameter for lam.z, default is 1}

\item{blamz}{prior shape2 parameter for lam.z, default is 1}

\item{asig}{prior shape parameter for sigma2, default is 1}

\item{bsig}{prior scale parameter for sigma2, default is 1}

\item{atau}{prior shape parameter for tau, default is 1}

\item{btau}{prior scale parameter for tau, default is 1}

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
Need to describe the function's output
\item{stuff}{If the function outputs a list of that contains stuff, then you can list things here.}
}

\references{
}
\examples{

# Need to put an example

}

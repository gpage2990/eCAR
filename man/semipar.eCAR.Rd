\name{semipar.eCAR}
\alias{semipar.eCAR}
\title{R wrapper that accesses 'INLA' to fit semi-parametric CAR spatial regression model that accounts for spatial confounding.}
\description{
\code{semipar.eCAR} is the main function used to fit the semi-parametric CAR model specified in the spectral domain.  This function calls 'INLA'.
}
\usage{
semipar.eCAR(y, x, W, E,
                       C=NULL,
                       names.covariates=NULL,
                       model="Gaussian",
                       eCAR.model="besag",
                       L=10, rw.ord="rw1",
                       pcprior.sd=c(0.1,1),
                       s2=10,
                       method = "spectral",
                       num.threads.inla = 1,
                       verbose=FALSE, ...)
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

\item{names.covariates}{Specifies the names of the covariates inside C}

\item{model}{Specifies the likelihood or data model.  Options are "Gaussian", "Poisson", "Binomial", "Negative Binomial"}

\item{eCAR.model}{Specifies the assumed gmrf for the spatial residuals.  Options are "besag" and "bym", the latter implementing the BYM model under Riebler/Dean reparametrization (this corresponds to model='bym2' in inla). Default is 'besag'.}

\item{L}{Number of basis functions for the spline model on the (spatial scale)-varying beta. The smoothing method applied here is a Bayesian version of the P-spline approach by Eilers and Marx (1996), assuming a random walk on the spline coefficients and a PC-prior on the precision parameter of the random walk.}

\item{rw.ord}{The order of the random walk prior assumed on the spline coef for the (spatial scale)-varying beta. Default is 'rw1'.}

\item{pcprior.sd}{Vector of length 2 specifying the scaling parameters for the PC-priors assumed on the precision of the (spatial scale)-varying beta and the data y, respectively. Each of the scaling parameters can be interpreted as a guess on the marginal standard deviation (default are 0.1 and 1).}

\item{s2}{Prior variance for the log of the dispersion parameter (only used for model="Negative Binomial", default equal to 10).}

\item{method}{A character defining the type of adjustment; either "spectral" (default choice) which implements the model assuming (spatial scale)-varying beta, or "naive" which implements the standard method with constant beta hence no spectral adjustment.}

\item{verbose}{logical; if TRUE the verbose output from the "inla" call is printed.}

\item{num.threads.inla}{Argument that indicates the number of computing cores that the INLA call will occupy.  For syntax, see ``inla.setOption"}

\item{...}{Arguments to be passed to the "inla" call; for instance control.inla=list(strategy="laplace")}}


\value{
A eCAR object which is a list containing the following
\item{data_model}{Character indicating which model was fit}

\item{beta_omega}{Matrix that contains respectively, the posterior mean lower and upper quantiles of the (spatial scale)-varying beta at each omega value (for the non Gaussian cases it is the exponentiated beta).}

\item{posterior_draws}{List containing posterior draws of the following parameters

\enumerate{
\item{postsample.beta: matrix of dim Lx1000, containing drawns from the posterior distribution of the L spline coefficients}
\item{postsample.prec.beta: vector of length 1000, containing drawns from the posterior distribution of the precision of the random walk on the spline coefficients}
\item{postsample.prec.z: vector of length 1000, containing drawns from the posterior distribution of prec.z}
\item{postsample.lambda.z: vector of length 1000, containing drawns from the posterior distribution of lambda.z}
\item{postsample.c.beta: vector of length 1000, containing drawns from the posterior distribution of the bias corrected beta (i.e. beta evaluated at the largest eigenvalue)}
}
}

\item{DIC}{Deviance information criterion computed by 'INLA'}

\item{regrcoef}{posterior summaries (mean, sd, 0.025quant, 0.975quant) for the regression coefs associated to the covariates inside 'C' (if model is Binomial, Poisson or Negative binomial, the posterior summaries refer to the exponentiated coef)}

}

\references{
Guan, Y; Page, G.L.; Reich, B.J.; Ventrucci, M.; Yang, S; "A spectral adjustment for spatial confounding" <arXiv:2012.11767>
}

\examples{

# Our R-package
library(eCAR)

data(lipcancer)
W <- lipcancer$neighborhood.Matrix
Y <- lipcancer$data$observed
X <- lipcancer$data$pcaff
E <- lipcancer$data$expected

\donttest{
# only run example if INLA is installed.
if (requireNamespace("INLA", quietly = TRUE)) {
  fit1 = semipar.eCAR(y=Y, x=X, W=W, E=E, C=NULL,
                                    pcprior.sd = c(0.1,1),
                                    model="Poisson",
                                    eCAR.model="bym",
                                    L=10, rw.ord="rw1",
                                    num.threads.inla = '1:1',
                                    verbose=FALSE)
  plot(fit1)
}
}
}

\name{semipar.eCARglm.Leroux}
\alias{semipar.eCARglm.Leroux}
\title{R wrapper that accesses INLA to fit semi-parametric Leroux CAR spatial regression model that accounts for spatial confounding.}
\description{
\code{semipar.eCARglm.Leroux} is the main function used to fit the semi-parametric Leroux CAR model specified in the spectral domain.  This function calls INLA.
}
\usage{
semipar.eCARglm.Leroux(y, x, W,
                       E,
                       C=NA,
                       model="Gaussian",
                       L=20,
                       pcprior.sd=c(0.1,1),
                       s2=10,
                       eval.fineGrid=FALSE,
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

\item{model}{Specifies the likelihood or data model.  Options are "Gaussian", "Poisson", "Binomial", "Negative Binomial"}

\item{L}{Number of basis functions for the spline model on the (spatial scale)-varying beta. The smoothing method applied here is a Bayesian version of the P-spline approach by Eilers and Marx (1996), assuming a random walk on the spline coefficients and a PC-prior on the precision parameter of the random walk.}

\item{pcprior.sd}{Vector of length 2 specifying the scaling parameters for the PC-priors assumed on the precision of the (spatial scale)-varying beta and the data y, respectively. Each of the scaling parameters can be interpreted as a guess on the marginal standard deviation (default are 0.1 and 1).}

\item{s2}{Prior variance for the log of the dispersion parameter (only used for model="Negative Binomial", default equal to 10).}

\item{eval.fineGrid}{logical; if TRUE the (spatial scale)-varying beta is computed over a grid of 1000 equally spaced locations on the eigenvalue domain.}

\item{verbose}{logical; if TRUE the verbose output from INLA is printed.}

\item{...}{Arguments to be passed to the inla call; for instance control.inla=list(strategy="laplace")}}


\value{
A list containing

\item{omega}{vector containing eigenvalues (if eval.fineGrid is set to TRUE, omega is a 1000-dimensional vector otherwise is n-dimensional)}

\item{beta.mn,beta.q025,beta.q975}{vectors containing,respectively, the posterior mean lower and upper quantiles of the (spatial scale)-varying beta at each omega value (for the non Gaussian cases it is the exponentiated beta).}

\item{B.pred}{the B-spline prediction matrix, containing L B-spline basis functions evaluated at each omega value (if eval.fineGrid is TRUE, B-splines are evaluated on 1000 equally-spaced locations on the eigenvalue domain)}

\item{postsample.beta}{matrix of dim Lx1000, containing drawns from the posterior distribution of the L spline coefficients}

\item{postsample.prec.beta}{vector of length 1000, containing drawns from the posterior distribution of the precision of the random walk on the spline coefficients}

\item{postsample.prec.z}{vector of length 1000, containing drawns from the posterior distribution of prec.z}

\item{postsample.lambda.z}{vector of length 1000, containing drawns from the posterior distribution of lambda.z}

\item{postsample.c.beta}{vector of length 1000, containing drawns from the posterior distribution of the bias corrected beta (i.e. beta evaluated at the largest eigenvalue)}
}

\references{
Will put something here if it is relevant.
}

\examples{

# Our R-package
library(eCAR)

data(lipcancer)
W <- lipcancer$neighborhood.Matrix
Y <- lipcancer$data$observed
X <- lipcancer$data$pcaff
E <- lipcancer$data$expected

fit1 = semipar.eCARglm.Leroux(y=Y, x=X, W=W, E=E, C=NULL,
                                    pcprior.sd = c(0.1,1),
                                    model="Poisson",
                                    L=20,
                                    eval.fineGrid = FALSE,
                                    verbose=FALSE)
plot.eCAR(fit1, exp.beta = TRUE, ylim=c(0.8,1.3))
}

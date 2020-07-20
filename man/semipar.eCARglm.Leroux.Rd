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
                       pcprior.sd=c(0.1,1),
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

\item{L}{Number of basis.}

\item{pcprior.sd}{To finish}

\item{s2}{To finish}

\item{eval.fineGrid}{To finish}

\item{verbose}{To finish}
}
\value{
Need to describe the function's output
\item{stuff}{If the function outputs a list of that contains stuff, then you can list things here.}
}

\references{
Will put something here if it is relevant.
}
\examples{

# Need to put an example

}

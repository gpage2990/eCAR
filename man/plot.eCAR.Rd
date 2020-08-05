\name{plot.eCAR}
\alias{plot.eCAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the spatial scale varying coefficient
}
\description{This function takes the output obtained from the parametric or semiparametric fit and retunrs a plot of the spatial scale varying coefficient.
}
\usage{
plot.eCAR(obj.eCAR,
          exp.beta = FALSE,
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.eCAR}{
an object of class eCAR (i.e. the output of the par.eCAR.Leroux() or semipar.eCAR.Leroux()
}
  \item{exp.beta}{
logical; if TRUE the exponential of the spatial scale varying coefficient is plotted (this is useful in interpretation of covariate effect in Binomial, Negative Binomial and Poisson models)
}
  \item{\dots}{
  include here other inputs to the plot function
}
}
\value{
This function returns the estimated posterior mean and 95-th credible intervals for the effect of the covariate of interest as a function of eigenvalues.
}
\references{
TO DO...put references to the literature/web site here
}

\examples{
#TO DO...
}

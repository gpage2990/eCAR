\name{plot}
\alias{plot.eCAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the spatial scale varying coefficient
}
\description{This function takes the output obtained from the parametric or semiparametric fit and returns a plot of the spatial scale varying coefficient.
}
\usage{
\method{plot}{eCAR}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class eCAR (i.e. the output of the par.eCAR.Leroux() or semipar.eCAR.Leroux())
}
  \item{\dots}{
  include here other inputs to the plot function
}
}
\value{
This function returns the estimated posterior mean and 95-th credible intervals for the effect of the covariate of interest as a function of eigenvalues.   If model is not Gaussian the exponential of the spatial scale varying coefficient is plotted which is useful in interpretating the covariate effect in the Binomial, Negative Binomial and Poisson models.

}

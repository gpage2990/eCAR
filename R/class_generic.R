#' eCAR class constructor
#'
#' A constructor for the \code{eCAR} class. The class \code{eCAR} is a named list containing
#' the output from the calling the \code{par.eCAR.Leroux} or \code{semipar.eCAR.Leroux} functions.
#'
#' @param data_model a characther indicating what data model was fit;
#' @param beta_omega matrix containing estimated beta as a function of omega with 95\% credible bands, and eigen-values;
#' @param posterior_draws a list containing the posterior draws of all model parameters;
#' @param DIC Deviance information criterion;
#' @param regrcoef posterior summaries (mean, sd, 0.025quant, 0.975quant) for the regression coefs;
#'
#' @export

eCAR.out <- function(
  data_model = NULL,
  beta_omega = NULL,
  posterior_draws = NULL,
  DIC = NULL,
  regrcoef =NULL
){
  value <- list(data_model = data_model,
                beta_omega = beta_omega,
                posterior_draws = posterior_draws,
                DIC=DIC,
                regrcoef = regrcoef)
  attr(value, "class") <- "eCAR"
  value
}


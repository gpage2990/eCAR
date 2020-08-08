plot.eCAR <- function(x,...) {
  model <- x$data_model
  omega <- x$beta_omega[,4]
  beta.mn <- x$beta_omega[,1]
  beta.q025 <- x$beta_omega[,2]
  beta.q975 <- x$beta_omega[,3]
  b0 <- rep(1, length(omega))
  if(model == "Gaussian") b0 <- rep(0, length(omega))
  graphics::matplot(omega, cbind(beta.mn, beta.q025, beta.q975, b0),
                   type="l", lty=c(1,2,2,3), col=c(1,1,1,2),
                   xlab=expression(omega), ylab=expression(beta(omega[k])), ...)

}


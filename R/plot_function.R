

plot.eCAR = function(obj.eCAR, exp.beta=FALSE, ...){
  omega = obj.eCAR$omega
  B = obj.eCAR$B.pred
  if (!exp.beta) {
    beta.mn = apply(B%*%obj.eCAR$postsample.beta, 1, mean)
    beta.q025 = apply(B%*%obj.eCAR$postsample.beta, 1, quantile, 0.025)
    beta.q975 = apply(B%*%obj.eCAR$postsample.beta, 1, quantile, 0.975)
  } else {
      beta.mn = apply(exp(B%*%obj.eCAR$postsample.beta), 1, mean)
      beta.q025 = apply(exp(B%*%obj.eCAR$postsample.beta), 1, quantile, 0.025)
      beta.q975 = apply(exp(B%*%obj.eCAR$postsample.beta), 1, quantile, 0.975)
  }
  matplot(omega, cbind(beta.mn, beta.q025, beta.q975, rep(1, length(omega))),
                   type="l", lty=c(1,2,2,3), col=c(1,1,1,2),
                   xlab=expression(omega), ylab=expression(beta(omega[k])), ...)
  #return(p.plot)
}


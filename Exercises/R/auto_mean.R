#=========================================================================
#
#  Program to compare the efficiency properties of MLE and OLS
#  in the AR(1) regression model in the case where the
#  explanatory variable is a constant. In this case the OLS estimator is
#  the sample mean
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

rho <- 0.6
sig2 <- 10
T    <- c(5, 50, 500)

#   Compute OLS variance
for(k in seq(T)) {
  t <- T[k]
  sum <- 0.0
  for (i in seq(t-1)) {
    for (j in seq(t-i)) {
          sum <- sum + rho^i  
    }    
  }
  var_ols <- sig2*(t + 2*sum)/((1-rho^2)*t^2)

  #  OLS estimator asymptotic variance based on Harvey (1990, p.197)
  sum_h <- (t*(1-rho^2) - 2*rho*(1-rho^t))/((1-rho)^2)
  var_ols_harvey <- sig2*(sum_h)/( (1-rho^2)*t^2 )

  # Compute MLE variance
  var_mle <- sig2/((t-1)*(1-rho)^2)

  #  Print results
  cat('\nSample size              = ', t)
  cat('\nVariance of OLS          = ', var_ols)
  cat('\nVariance of OLS (Harvey) = ', var_ols_harvey)
  cat('\nVariance of MLE          = ', var_mle)
  cat('\nEfficiency (OLS/MLE)     = ', var_ols/var_mle, '\n')  
}
    


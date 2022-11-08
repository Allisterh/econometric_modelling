# ===========================================================================
#
#      Program to generate GMM table and demonstrate empirical moments
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# 
#----------------------Method of Moments Estimation--------------------------
#

gmm_table <- function() {
  t <- 10                        

  yt <- as.matrix(read.table("table.dat"))
  
  m1 <- mean(yt)                 
  
  mt <- cbind(yt, (yt^2), (yt^3), (yt^4), (yt - m1), (yt - m1)^2, log(yt), (1/yt))
  
  m <- colMeans(mt)
  print(unname(rbind(mt, m)))
  
  
  cat('\nNormal distribution (mu)    = ', m[1])
  cat('\nNormal distribution (sig2)  = ', m[6])
  
  cat('\n ' )
  cat('\nStudent t distribution (mu) = ', m[1])
  cat('\nStudent t distribution (nu) = ', 2*m[6]/(m[6]-1))
  
  
  # Fourth central moment 
  # (can be computed directly as c4 = mean((yt-m(1)).^4)
  c4 <- m[4]-4*m[3]*m[1]+6*m[2]*m[1]^2-4*m[1]*m[1]^3+m[1]^4           
  
  cat('\n' )
  cat('\nStudent t distribution (mue) = ', m[1])
  # Solution based on solving a quadratic equation in nu
  cat('\nStudent t distribution (nu)\n' )                            
  cat('\n     first root              = ', (6+sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4)))    
  cat('\n     second root             = ', (6-sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4)))    
  
  cat('\n' )
  cat('\nGamma distribution (alpha)   = ', m[1]/(m[1]-m[8]))
  cat('\nGamma distribution (beta)    = ', (m[1]/(m[1]-m[8]))/m[1])
  
  cat('\n ' )
  cat('\nGamma distribution (alpha)   = ', m[1]^2/(m[2]-m[1]^2))
  cat('\nGamma distribution (beta)    = ', m[1]/(m[2]-m[1]^2))
  
  cat('\n ' )
  cat('\nPareto distribution (alpha)  = ', m[1]/(m[1]-min(yt))) 
}

#=========================================================================
#
#    Generate simulated data by simulating the reduced form.
#    The set of equations is given by the bivariate system
#        y1t = beta1*y2t + alpha1*x1t + u1t
#        y2t = beta2*y1t + alpha2*x2t + u2t
#     where E[ut'ut] = omega
#
#=========================================================================


rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -----------------------------------
#
source("EMTSUtil.R")
library("scatterplot3d")

#
#--------------------------- Simulating a Simultaneous System -----------------------------------
#

linear_simulation <- function() {
  t <- 500 

  beta1  <- 0.6 
  beta2  <- 0.2 
  alpha1 <- 0.4
  alpha2 <- -0.5
  
  omega <- matrix(c(1, 0.5, 0.5, 1.0), nrow = 2, byrow=T)
  
  # Construct population parameter matrices     
  B <- matrix(c(1, -beta2, -beta1, 1.0), nrow = 2, byrow=T) 
  A <- matrix(c(-alpha1, 0, 0, -alpha2), nrow = 2, byrow=T)
  
  # Construct exogenous variables                       
  x <- cbind( 10*rnorm(t), 3*rnorm(t) )
  
  # Construct structural disturbances                   
  u <- matrix(rnorm(t*2), ncol = 2) %*% chol(omega)
  
  # Construct reduced form parameters                   
  invB <- inv(B)
  phi  <- -A %*% invB
  
  cat('\nInverse of B\n')
  print(invB)
  cat('\nReduced form parameter matrix\n')
  print(phi)
  
  # Construct reduced form disturbances                 
  v <- u %*% invB
  
  # Simulate the model by simulating the reduced form   
  y <- array(0, c(t,2))
  
  for (i in seq(t)) {
    y[i,] <- -x[i,] %*% A %*% invB + v[i,]
  }
  
  
  #-----------------------------------------------------------------
  # Generate Graphs
  #-----------------------------------------------------------------
  
  figure()
  t <- seq(t)
  
  par(mfrow=c(2,2))
  
  # Panel (a)
  plot(t,y[,1], type="l",  xaxs="i", yaxs="i",
       main ="(a)",
       xlim = c(0,500),
       ylim = c(-10, 10),
       xlab = expression(t),
       ylab = expression(paste("y"[1]*"",",", ""[t]*"") ))
  
  # Panel (b)
  plot(t,y[,2], type="l",  xaxs="i", yaxs="i",
       main ="(b)",
       xlim = c(0,500),
       ylim = c(-10, 10),
       xlab = expression(t),
       ylab = expression(paste("y"[2]*"",",", ""[t]*"") ))
  
  # Panel (c)
  scatterplot3d(y[,1], x[,1], y[,2], highlight.3d=TRUE,      
        main="(c)", pch=19, 
        xlab = expression(paste("y"[1]*"",",", ""[t]*"")),
        ylab = expression(paste("x"[1]*"",",", ""[t]*"")),
        zlab = expression(paste("y"[2]*"",",", ""[t]*"")),
        xlim = c(-10,10),
        ylim = c(-10, 10),
        zlim = c(-10, 10),
        box = F, grid = F, angle=55, scale.y=0.7)
  
  # Panel (d)
  scatterplot3d(y[,2], x[,2], y[,1], highlight.3d=TRUE,   
        main="(d)", pch=19, 
        xlab = expression(paste("y"[2]*"",",", ""[t]*"")),
        ylab = expression(paste("x"[2]*"",",", ""[t]*"")),
        zlab = expression(paste("y"[1]*"",",", ""[t]*"")),
        xlim = c(-10,10),
        ylim = c(-10, 10),
        zlim = c(-10, 10),
        box = F, grid = F, angle=55, scale.y=0.7)

  
}


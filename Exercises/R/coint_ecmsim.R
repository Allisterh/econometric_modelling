#============================================================================
#
#   Program to generate simulated time series plots for 
#   alternative specifications of the vecm
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - seqa, figure
source("EMTSUtil.R")

# Load require library - Null
library(MASS) 



t <- 200            

experiment <- 4                                                  

if (experiment == 1) {
  beta  <- c(1, -1)                 # Normalized cointegrating vector                              
  beta0 <- 2                      # Long-run constant                                            
  beta1 <- 0                      # Long-run time trend                                          
  
  alpha <- c(-0.1,0.1)             # Error-correction parameters                                  
  
  alpha0 <- Null(alpha %*% t(alpha))
  delta0 <- alpha0*0              # Short-run intercepts 
  delta1 <- alpha0*0              # Short-run time trend using orthogonal complement of alpha                     
} else if (experiment == 2){
  beta <- c(1,-1)                  # Normalized cointegrating vector                              
  beta0 <- 2                      # Long-run constant                                            
  beta1 <- 0                      # Long-run time trend                                          
  
  alpha <- c(-0.0,0.01)           # Error-correction parameters                                  
  
  alpha0 <- Null(alpha %*% t(alpha))
  delta0 <- alpha0*0              # Short-run intercepts using orthogonal complement of alpha    
  delta1 <- alpha0*0              # Short-run time trend using orthogonal complement of alpha    
} else if (experiment == 3) {
  beta <- c(1,-1)                  # Normalized cointegrating vector                             
  beta0 <- 2                      # Long-run constant                                          
  beta1 <- 0                      # Long-run time trend                                        
  
  alpha <- c(-0.1,0.1)             # Error-correction parameters                               
  
  alpha0 <- Null(alpha %*% t(alpha))
  delta0 <- alpha0*2              # Short-run intercepts using orthogonal complement of alpha   
  delta1 <- alpha0*0              # Short-run time trend using orthogonal complement of alpha 
} else {
  beta <- c(1,-1)                  # Normalized cointegrating vector                             
  beta0 <- 2                      # Long-run constant                                           
  beta1 <- 0.1                    # Long-run time trend                                         
  
  alpha <- c(-0.1,0.1)             # Error-correction parameters                                
  
  alpha0 <- Null(alpha %*% t(alpha))
  delta0 <- alpha0*0              # Short-run intercepts using orthogonal complement of alpha  
  delta1 <- alpha0*0              # Short-run time trend using orthogonal complement of alpha  
}
# Simulate the model   

y <- array(0, c(t,2))
v <- array(rnorm(t*2), c(t,2))

for (j in 2:t) {
  u <- beta0 + beta1*j + y[j-1,] %*% beta
  y[j,] <- y[j-1,] + t(delta0 + delta1*j + alpha*u) + v[j,]
}


# Graph simulated data 
figure()

matplot(seq(t),y, type="l",
        xlab = "t",
        ylab = "",
        lty=c(1,4),
        col=c(3,5))
legend("topright",
       legend = c("y1", "y2"),
       lty=c(1,4),
       col=c(3,5))


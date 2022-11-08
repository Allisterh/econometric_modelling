#============================================================================
#
#     Program to compute generalized impulse response function 
#     (Koop, Pesaran and Potter (1996), Journal of Econometrics)
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - recserar
source("EMTSUtil.R")

# load required functions - reshape, zeros
library("matlab") 
#----------------------------------------------------------------------------
#  Model used to compute GIRF
#----------------------------------------------------------------------------
model <- function(v,y0,rho1,rho2) {
  y <- y0 + rep(0, length(v))
  
  for (i in 2:length(y)) {
    # Linear model 
    #         rho = rho1+rho2
    #         y(i) = rho*y(i-1) + v(i)           
    
    # Bilinear model     
    y[i] <- rho1*y[i-1] + rho2*y[i-1]*(y[i-1]>=0) + v[i]
  }
  return(y)
}

#
#--------------------------- GIRF Function ----------------------------------
#
nlm_girf <- function() 
{
  # Parameters
  delta <- -1             # Size of the shock     
  rho1  <- 0.25           # AR(1) parameter                  
  rho2  <- 0.50           # TAR parameter                    
  rho   <- rho1 + rho2
  t     <- 1000           # Choose sample size       
  n     <- 10             # Forecast horizon of impulses     
  
  # Simulate the data   
  v <- rnorm(t)
  y <- rep(0, t)
  
  for (i in 2:t) {
    y[i] <- rho1*y[i-1] + rho2*y[i-1]*(y[i-1]>=0) + v[i]
  }
  
  # Compute generalized impulse response function 
  # average over horizon within loop for a particular history  
  # and then average over histories
  
  # Total number of draws needed for a given initial condition: 
  # t-1 from history and n+1 from number of impulse horizons       
  tn <- (t-1)*(n+1)                                          
  
  impulse <- zeros(t-1,n+1)
  pb <- txtProgressBar(min=0, max=t, style=3)
  # Loop through the data to change the initial condition (history)  
  
  for (i in seq(t-1)) {
    #  Bootstrap residuals 
    ind    <- trunc( runif(tn)*(t-1) + 1 )
    v_boot <- v[ind]              
    v_boot <- reshape(cbind(v_boot),n+1,t-1)
    
    ye0 <- zeros(n+1,t-1)
    ye1 <- zeros(n+1,t-1)
    
    # Loop through horizon of impulse
    for (j in seq(t-1)) {
      # Initial condition based on a boostrap draw              
      ye0[,j] <- model(v_boot[,j],v_boot[1,j],rho1,rho2)        
      
      # Initial condition based on history (i subscript) plus 1  
      ye1[,j] <- model(v_boot[,j],v[i]+delta,rho1,rho2)
    } 
    
    # Average over horizon given an initial condition (history)          
    impulse[i,] <- (rowMeans(ye1) - rowMeans(ye0))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Average each impulse across histories (ie initial conditions)     
  impulse_girf <- cbind(colMeans(impulse))
  
  
  # Linear impulse response    
  impulse_linear <- delta*recserar( cbind(rep(0, 11)) , cbind(1) , cbind(rho) )   
  
  cat('\n         GIRF      Linear ')
  print(cbind(impulse_girf, impulse_linear))
}




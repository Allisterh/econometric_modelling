#============================================================================
#
#   Program to demonstrate the properties of an artificial neural network
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123457, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#
# load required functions - trimr, figure
source("EMTSUtil.R")


#
#----------------------- Artificial Neural Networks -------------------------
#
nlm_annfig <- function () {
  t <- 200      
  
  # Set the parameters of the ANN      
  phi    <- 0.4           
  gam    <- 2.0
  delta0 <- -2.0
  delta1 <- 2.0
  
  
  # Create predictions      				     
  ylag <- seq(-4, 4, 0.1)
  
  f <- 1/( 1+exp(-( delta0 + delta1*ylag)) )
  
  ylin  <- phi*ylag
  ynlin <- gam*f
  
  yhat <- ylin + ynlin  
    
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  
  matplot(ylag,cbind(yhat,ylin, ynlin), type="l",
          main = "",
          ylab = expression(y[t]),
          xlab = expression(y[t-1]),
          bty = "l")  
  
}

#=========================================================================
#
#   Program to estimate the Hamilton and Jorda (2002)
#   ordered probit model of US monetary policy.
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  inv
source("EMTSUtil.R")

#-------------------------------------------------------------------------
#  Unrestricted Probit negative log-likelihood function 
#-------------------------------------------------------------------------
lprobit <- function(b,x,d) {
  # Cut off points
  c <- b[3:6]
  
  # Regression part excluding the intercepts
  xb <- x %*% b[1:2]  
  # Cut off points
  f1 <- pnorm( c[1] - xb,0,1)
  f2 <- pnorm( c[2] - xb,0,1) - f1
  f3 <- pnorm( c[3] - xb,0,1) - f1 - f2
  f4 <- pnorm( c[4] - xb,0,1) - f1 - f2 - f3
  f5 <- 1 - f1 - f2 - f3 - f4
  f  <- cbind(f1,  f2,  f3,  f4,  f5 ) 
  
  # Log-likelihood function
  tp <- d * log(f)
  lf <- -mean( rowSums(tp))    
  
  return(lf) 
}


#
#---------------- Hamilton-Jorda Ordered Probit Model  ----------------------
#
discrete_hamilton_jorda <- function() {
  # Read the data: US weekly yields   
  # 1st week of February 1984 to last week of April 2001
  load("hamilton_jorda.Rdata")
  data <- hamilto_jorda
  
  event   <- data[,1]
  target  <- data[,2]
  change  <- data[,3]
  bin     <- data[,4]
  spread6 <- data[,5]
  
  
  for (i in 2:length(event)) {
    if (event[i] == 1)
      change[i] <- change[i]
    else
      change[i] <- change[i-1]  
  }
  
  # Convert data to event time
  capt      <- length(event)
  eventdata <- array(0, c(capt,2))
  ybin      <- array(0, c(capt,1))
  
  lagchange <- 0
  i <- 1  
  
  for (t in seq(capt)) {    
    if (event[t] == 1) {
      eventdata[i,1] <-  lagchange
      eventdata[i,2] <- spread6[t]
      ybin[i,1]      <- bin[t]
      
      if (t > 1)
        lagchange <- target[t] - target[t-1]
      i <- i + 1
    } 
  }
  
  summ <- sum(event)
  ybin <- ybin[1:summ] 
  x    <- eventdata[1:summ,]
  
  # Create dummy variables for each interest rate change
  d1 <- as.numeric(ybin == -0.50)
  d2 <- as.numeric(ybin == -0.25)
  d3 <- as.numeric(ybin ==  0.00)
  d4 <- as.numeric(ybin ==  0.25)
  d5 <- as.numeric(ybin ==  0.50)
  
  d  <- cbind(d1, d2, d3, d4, d5)
  
  # Choose event days from 1984 to 1997
  x <- x[1:102,]         
  d <- d[1:102,]
  
  # Estimate the ordered probit model
  # Use MATLAB Data
  theta0 <- c(2.5449149,0.54142729, -1.8948826, -0.42001991,-0.0052515480,1.5173916)
  
  estResults <- optim(theta0, lprobit, x=x, d=d, method="BFGS", hessian=T)
  theta1 <- estResults$par
  l1 <- estResults$val
  h <- estResults$hessian
  l1 <- -l1
  
  cat('\nUnrestricted log-likelihood function =     ',l1)
  cat('\nT x unrestricted log-likelihood function = ',t*l1)
  
  cat('\nParameter estimates ', theta1)
  
}

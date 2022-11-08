#============================================================================
#
#   Estimate the Dynamic Conditional Correlation (DCC) model of US yields
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, figure, inv, reshapeg
source("EMTSUtil.R")

#----------------------------------------------------------------------------
# Likelihood function for a GARCH(1,1) model
#----------------------------------------------------------------------------
negloglgarch <- function( b,y ) {
  u <- y       
  h <- recserar(cbind(b[1] + pnorm(b[2])*trimr(c(0.0,u^2),0,1)),cbind(sd(u)^2),cbind(pnorm(b[3])))
  z <- u/sqrt(h)
  f <- - 0.5*log( 2*pi ) - 0.5*log( h ) - 0.5*z^2
  
  lf <- -mean( f )
  return(lf)
}

#----------------------------------------------------------------------------
# Correlation component of DCC
#----------------------------------------------------------------------------
neglogcor <- function(b,y,h1) {
  b    <- pnorm( b ) 
  t <- nrow(y)
  
  f <- rep(0,t)
  u <- y
  z  <- u/sqrt(h1)
  
  qbar <- t(z) %*% z/t                                          
  q    <- qbar
  
  for (i in seq(t)) {
    # Diagonal matrix of conditional standard deviations   
    s  <- diag(sqrt(h1[i,]))   
    
    # Conditional correlation matrix  
    tmp <- inv(diag(sqrt(diag(q))))
    r   <- tmp %*% q %*% tmp
    
    f[i] <- -0.5*log(det(r))-0.5*z[i,] %*% inv(r) %*% cbind(z[i,])+0.5*z[i,] %*% cbind(z[i,])
    
    # Update q
    q  <- abs(1-b[1]-b[2])*qbar + b[1]*cbind(z[i,]) %*% z[i,] + b[2]*q
  }
  lf <- -mean( f )
  return(lf)  
}


#----------------------------------------------------------------------------
# Full log-likelihood for DCC model
#----------------------------------------------------------------------------
neglogf <- function(b,y) {
  ind <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 16, 17)
  b[ind] <- pnorm(b[ind])
  
  t <- nrow(y)
  n <- ncol(y)
  f <- rep(0, t)
  u <- y                                       
  hv <- array(0, c(t,n))
  
  # Construct conditional variances
  hv[,1] <- recserar(cbind(b[1]+b[2]*trimr(c(0.0,u[,1]^2),0,1)),cbind(sd(u[,1]))^2,cbind(b[3]))
  hv[,2] <- recserar(cbind(b[4]+b[5]*trimr(c(0.0,u[,2]^2),0,1)),cbind(sd(u[,2])^2),cbind(b[6]))
  hv[,3] <- recserar(cbind(b[7]+b[8]*trimr(c(0.0,u[,3]^2),0,1)),cbind(sd(u[,3])^2),cbind(b[9]))
  hv[,4] <- recserar(cbind(b[10]+b[11]*trimr(c(0.0,u[,4]^2),0,1)),cbind(sd(u[,4])^2),cbind(b[12]))     
  hv[,5] <- recserar(cbind(b[13]+b[14]*trimr(c(0.0,u[,5]^2),0,1)),cbind(sd(u[,5])^2),cbind(b[15]))    
  
  z    <- u/sqrt(hv)                                                                                           
  qbar <- t(z) %*% z/t                                          
  q    <- qbar
  for (i in seq(t)) {
    # Diagonal matrix of conditional standard deviations   
    s  <- diag(sqrt(hv[i,]))   
    
    # Conditional correlation matrix  
    tmp <- inv(diag(sqrt(diag(q))))
    r   <- tmp %*% q %*% tmp
    
    # Update conditional variance-covariance matrix   
    h  <- s %*% r %*% s                                                                             
    
    f[i] <- -0.5*n*log(2*pi)-0.5*log(det(h))-0.5*u[i,] %*% inv(h) %*% cbind(u[i,])

    # Update q
    q  <- abs(1-b[16]-b[17])*qbar + b[16]*cbind(z[i,]) %*% z[i,] + b[17]*q    
  }
  lf <- -mean( f )
  return(lf)
}



#
# ------------------------ MGARCH DCC Model ----------------------------------
#
mgarch_dcc <- function() {
  # Load data: US daily yields 3-Jan-2000 to 21-Aug-2006
  data <- as.matrix(read.table("daily_finance.dat"))
  
  data <- data[,c(1, 5, 10, 15, 20)]
  
  y    <- 100*(trimr(data,1,0) - trimr(data,0,1))                                   
  y    <- y - mean(y)
  
  t <- nrow(y)
  n <- ncol(y)
  
  # Estimate univariate GARCH models for each variable   
  b1   <- array(0, c( 3,n ))
  st1  <- array(0, c(3,n))
  h1   <- array(0, c( t,n ))
  
  for (k in seq(n)) {
    bstart      <- c(0.1, -2,  2)
    estResults <- optim(bstart, negloglgarch, y=y[,k], method="BFGS")
    tmp       <- estResults$par 
    st1[,k]  <- tmp
    b1[,k]   <- c(tmp[1], pnorm(tmp[2]), pnorm(tmp[3]))
    h1[,k]   <- recserar(cbind(tmp[1] + pnorm(tmp[2])*trimr(c(0.0,y[,k]^2),0,1)),cbind(sd(y[,k])^2),cbind(pnorm(tmp[3]))) 
  }
  
  b1[c(2, 3),] <- pnorm( b1[c(2, 3),] )
  cat('\nGarch parameters\n' )
  print( b1 )   
  
  # Estimate the correlation component
  start <- c(-2, 2)
  estResults <- optim(start, neglogcor, y=y, h1=h1, method="BFGS")
  bc <- estResults$par
  st2   <- bc
  bc    <- pnorm( bc )
  cat('\nCORR parameters\n' )
  print( bc )
  
  # Estimate the correlation component
  start <- c(c(st1), c(st2))
  estResults <- optim(start, neglogf, y=y, method="BFGS")
  bf <- estResults$par
  
  ind <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 16, 17)
  bf[ind] <- pnorm(bf[ind])
  print(cbind(cbind(c(b1), c(bc)),bf))
}

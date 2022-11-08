#============================================================================
#
#   Kalman filter program to estimate a business cycle model for Australia
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()


#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - trimr, inv
source("EMTSUtil.R")

#--------------------------------------------------------------------------
#   Wrapper function to set up and call the Kalman filter
#--------------------------------------------------------------------------
neglog <- function(b,y) {
  Lam <- matrix(c(b[1:6], rep(0, 6)), ncol=2)
  Phi <- matrix(c(b[7],  b[8],
                  1,    0), 2,2, byrow=T)
  
  R   <- diag(b[9:14]^2)
  Q <- diag(c(1,rep(0, (nrow(Phi)-1))))
  
  lf <- -mean( kalman(y,Phi,Lam,R,Q) )
  return(lf)
}

#--------------------------------------------------------------------------
# Kalman filter
#--------------------------------------------------------------------------
kalman <- function(y,Phi,Lam,R,Q) {
  # Allocate arrays  
  t <- nrow(y)
  n <- ncol(y)
  k         <- nrow(Q)
  lnl       <- rep(0, t)
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)*0.1 
 
  mt <- Lam %*% st
  vt <- Lam %*% pt %*% t(Lam) + R
  ut <- y[1,] - mt
  
  lnl[1] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
  Kal <- pt %*% t(Lam) %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt
  
  
  
  # Main loop over observations
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0     
    pt <- Phi %*% p0 %*% cbind(Phi) + Q          
    
    # Observation
    mt <- Lam %*% st
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i,] - mt
    
    # Construct log-likelihood function
    lnl[i] <- - 0.5*n*log(2*pi) - 0.5*log(det(vt)) - 0.5*t(ut) %*% inv(vt) %*% ut
    
    # Update      		
    Kal <- pt %*% t(Lam) %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt   
  }   
  return(lnl)
}


#--------------------------------------------------------------------------
#   Extract smoothed factor
#--------------------------------------------------------------------------
kfsmooth <- function(b,y) {
  # Unpack the parameter vector
  Lam <- matrix(c(b[1:6], rep(0, 6)), ncol=2)
  Phi <- matrix(c(b[7],  b[8],
                  1,    0), 2,2, byrow=T)
  
  R   <- diag(b[9:14]^2)
  Q <- diag(c(1,rep(0, (nrow(Phi)-1))))
 
  # Allocate arrays  
  t <- nrow(y)
  n <- ncol(y) 
  k       <- nrow(Q)
  s10     <- array(0, c(t,k) )    # st|t-1
  s11     <- array(0, c(t,k) )    # st|t  
  p10     <- array(0, c(k,k,t) )
  p11     <- array(0, c(k,k,t) )
  ss      <- rep(0,t)
  
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- diag(k)*0.1  
  s10[1,] <- st 
  p10[,,1] <- pt
  
  mt <- Lam %*% st
  vt <- Lam %*% pt %*% t(Lam) + R
  ut <- y[1,] - mt
  
  Kal <- pt %*% t(Lam) %*% inv(vt)
  s0 <- st + Kal %*% ut
  p0 <- pt - Kal %*% Lam %*% pt  
  
  s11[1,] <- s0
  p11[,,1] <- p0
  
  # Main loop over observations  
  for (i in 2:t) {
    # Prediction 
    st <- Phi %*% s0 
    pt <- Phi %*% p0 %*% cbind(Phi) + Q   
    s10[i,] <- st 
    p10[,,i] <- pt
    
    # Observation
    mt <- Lam %*% st
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i] - mt
    
    # Update          
    Kal <- pt %*% t(Lam) %*% inv(vt)
    s0 <- st + Kal %*% ut
    p0 <- pt - Kal %*% Lam %*% pt
    s11[i,] <- s0
    p11[,,i] <- p0
  }
  
  # Now smooth the factor    
  ss <- s11
  
  for (j in 1:(t-1)) {    
    Jt      <- p11[,,t-j] %*% t(Phi) %*% inv(p10[,,t-j])
    ss[t-j,] <- s11[(t-j),] + Jt %*% t(( t(ss[(t-j+1),]) - t(s10[(t-j+1),]) ) )
  }
  fac <- ss 
  return(fac)
}



#
#--------------------------- Business Cycle ---------------------------------
#
lfac_bcycle <- function( ) {
  # Load the data: 
  # Australian monthly business cycle data Jan 1980 - September 2009
  load('lfac_bcycle.Rdata')
  
  
  coincident   <- data[,1]        # Coincident index                                      
  gdp          <- data[,2]        # GDP (interpolated)               
  unemployment <- data[,3]        # Unemployment rate (p.a.)                             
  employment   <- data[,4]        # Employment                                           
  sales        <- data[,5]        # Retail sales                                          
  income       <- data[,6]        # Household income (interpolated)
  production   <- data[,7]        # Industrial production                                
  
  ind <- cbind(gdp, unemployment,  employment, sales, income, production)
  
  # Transformed indicators 
  y <- cbind(100*(trimr(log(ind[,1]),12,0) - trimr(log(ind[,1]),0,12)),
             trimr(ind[,2],12,0) - trimr(ind[,2],0,12),
             100*(trimr(log(ind[,3:6]),12,0) - trimr(log(ind[,3:6]),0,12)) )
  y <- t (apply(y, 1, '-', colMeans(y)))
  
  t <- nrow(y)
  
  
  # Estimate parameters 
  start  <- c(0.0917,      
              -0.0584,        
              0.0836,         
              0.0644,         
              0.0227,         
              0.0934,         
              1.8952,         
              -0.9125,         
              0.9376,         
              -0.2106,        
              0.7740,         
              2.8910,         
              3.0981,         
              3.8242)
  
  estResults <- optim(start, neglog, y=y, method="BFGS", hessian=T)
  bhat <- estResults$par
  lf <- estResults$val
  hess <- estResults$hess
  
  lf <- -lf
  vc <- (1/t)*inv(hess)
  cat('\nLog-likelihood function    = ',lf)
  cat('\n')
  print(cbind("Params"=bhat, "Std. Erros"=sqrt(diag(vc))) )
  
  # Compute business cycle based on coincident index  
  bc <- 100*(trimr(log(coincident),12,0) - trimr(log(coincident),0,12))      
  bc <- bc - mean(bc)
  
  # Compute and scale the smoothed factor estimate 
  fac <- kfsmooth(bhat,y)
  bm  <- fac[,1]
  bm  <- bm*sd(bc)/sd(bm)
  
  figure()
  matplot(seqa(1981+1/12,1/12,t), cbind(bm, -bc), type='l',
          main = '',
          xlab = '',
          ylab = '',
          bty= 'l')  
}




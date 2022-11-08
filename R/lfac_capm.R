#============================================================================
#
#   Program to estimate a a capital asset pricing model
#   and estimate the excess return on invested wealth
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
neglog <- function(b,y, flag) {
  alpha <- b[1:10]
  Lam <- b[11:20]
  if (flag)
    Phi <- tanh(b[21])
  else
    Phi <- b[21]
  
  R   <- diag(b[22:31]^2)
  Q <- diag(length(Phi))
  
  lf <- -mean( kalman(y,Phi,Lam,R,Q, alpha) )
  return(lf)
}

#--------------------------------------------------------------------------
# Kalman filter
#--------------------------------------------------------------------------
kalman <- function(y,Phi,Lam,R,Q, alpha) {
  # Allocate arrays  
  t <- nrow(y)
  n <- ncol(y)
  k         <- nrow(Q)
  lnl       <- rep(0, t)
  
  # Recursions of the Kalman Filter
  # Initialisation following Harvey and Hamilton  
  st <- rep(0, k)
  pt <- reshape(inv(diag(k^2) - Phi^2)*Q,k,k)  
  
  mt <- cbind(Lam) %*% st
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
    mt <- Lam %*% st + alpha
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
kfsmooth <- function(b,y, flag) {
  # Unpack the parameter vector
  alpha <- b[1:10]
  Lam <- b[11:20]

  if (flag)
    Phi <- tanh(b[21])
  else
    Phi <- b[21]
  
  R   <- diag(b[22:31]^2)
  Q <- diag(length(Phi))
  
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
  
  mt <- cbind(Lam) %*% st + alpha  
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
    mt <- Lam %*% st + alpha
    vt <- Lam %*% pt %*% t(Lam) + R
    ut <- y[i,] - mt
    
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
  return(list(fac=fac, s11=s11))
}


#
#--------------------------- CAPM Model -------------------------------------
#
lfac_capm <- function( ) {
  # Load the data: starting March 1990 and ending March 2010
  # 3 = Microsoft, 
  # 4 = Intel, 
  # 5 = Pfizer, 
  # 6 = Exxon, 
  # 7 = Procter & Gamble, 
  # 8 = Boeing 
  # 9 = AT&T, 
  # 10 = General Electric, 
  # 11 = Chevron, 
  # 12 = Bank of America
  load('capm.Rdata') 
  
  sp500    <- data[,1]               # S&P500 index                                                                                                            **/
  interest <- data[,2]               # US Treasury constant maturity, p.a., 3 months                                              **/
  price    <- data[,(3:ncol(data))]           # stocks
  
  # Compute excess asset excess returns
  tmp <- 4*(trimr(log(price),1,0) - trimr(log(price),0,1)) 
  
  y <- tmp - trimr(interest,1,0)/100
  
  # Compute excess market return on sp500
  market_sp500 <- 4*(trimr(log(sp500),1,0) - trimr(log(sp500),0,1)) - trimr(interest,1,0)/100                                
  
  t <- nrow(y)
  
  # Estimate parameters with restriction imposed
  flag  <- 1
  #start <- c(rep(1, 20), 1, rep(1, 10))
  
  # Use MATLAB results
  start <- c(0.148889067747353,
             0.108546962470647,
             0.0703614202700653,
             0.0532090874308006,
             0.0629648579184943,
             0.0124816228400908,
             -0.00453560600009894,  
             0.0218693408217819,
             0.0433138960520868,
             -0.00411564858472721,
             0.355381245302114,
             0.449074521815885,
             0.235070730893052,
             0.137134603145954,
             0.150166259021341,
             0.307405225575446,
             0.126865715939853,
             0.444433395146189,
             0.148129146775761,
             0.572089722354855,
             0.233608102554924,
             0.475325676639324,
             0.694125943630311,
             0.421990597701908,
             0.271776565325622,
             0.391122592543362,
             0.490866122650479,
             0.390265910866447,
             0.184390699705062,
             0.301287127120619,
             0.769322456168742)
  estResults <- optim(start, neglog, y=y, flag=flag, method="BFGS")
  bhat <- estResults$par
  
  # Restimate without restrictions to get standard errors
  flag  <- 0
  start <- c(bhat[1:20], tanh(bhat[21]),  bhat[22:length(bhat)])
  estResults <- optim(start, neglog, y=y, flag=flag, method="BFGS", hessian=T)
  bhat <- estResults$par
  lf <- estResults$val
  hess <- estResults$hess
  
  cat('\nLog-likelihood function = ',-lf)
  cat('\n')
  
  # Estimate excess returns on wealth - mean 
  s11 <- kfsmooth(bhat,y, flag)$s11
  
  # CAPM regression estimates  
  beta <- lm(y ~ cbind(rep(1, t), market_sp500) - 1)$coef
  
  cat('\nEstimates of Beta \n')
  KF.scaled <- bhat[11:20]*sd(as.numeric(s11))/sd(market_sp500)
  print( cbind("OLS"=beta[2,], "KF"=bhat[11:20], "KF(Scaled)"= KF.scaled))
  
}




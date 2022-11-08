#============================================================================
#
#      Program to estimate a svar with identification based on sign restrictions
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

# Load required functions -  trimr, reshapeg, figure, lag.matrix
source("EMTSUtil.R")

# Load required library - repmat, ones, zeros, eye
library("matlab")

#--------------------------------------------------------------------------
#   Compute impulse responses
#   itype=0 for non-accumulated and itype=1 for accumulated
#--------------------------------------------------------------------------
irf <- function(phi,s,itype,k,q,p) {
  # Construct SI(L) matrices (VMA)     
  si <- cbind(c(eye(k)))
  
  for (i in seq(q)) {
    sum <- 0
    for (j in 1:min(c(p, i))) {
       sum <- sum + reshapeg(phi[,j],k,k) %*% reshapeg(si[,(i-j+1)],k,k)      
    }
    si  <- cbind(si, as.vector(t(sum)))    
  }  
  # Construct orthogonal impulse response functions        
  impulse <- as.vector(t(s))  
  for (i in 2:(q+1)) {
    impulse  <- cbind(impulse, c( t(reshapeg(si[,i],k,k) %*% s )) )
  } 
  impulse <- t(impulse)
  if (itype)
    impulse <- apply(impulse, 2, cumsum)
  return(impulse)
}

#
#------------------------- Sign Restrictions --------------------------------
#   
svar_sign <- function () {
  # GDP and CPI monthly data from 1950Q1 to 2009Q4 for U.S.
  load('sign.Rdata')        

  # Define variables     
  rgdp <- ytdata[,1]
  cpi  <- ytdata[,2]

  # Construct variables for use in VAR: 
  yact <- cbind(log(rgdp),  log(cpi))
  y <- 400*(trimr(yact,1,0) - trimr(yact,0,1))

  # Choose parameters    
  p <- 2          # Order of VAR lags		
  q <- 40         # Order of VMA lags		

  # Estimate the VAR(p) with a constant and a time trend		
  ylag <- lag.matrix(y,1:p)  
  ylag <- ylag[-(1:p),]

  xmat <- cbind(rep(1, nrow(ylag)),  seq(nrow(ylag)),  ylag)
  ymat <- trimr(y,p,0)

  reg <- lm(ymat ~ xmat - 1)
  bar <- reg$coef
  nconst <- 2
  mu <- bar[1:nconst,]

  v   <- reg$residuals
	vc <- t(v) %*% v/nrow(v)

  cat('\nResidual variance-covariance matrix\n')
  print(vc)
  cat('\n')

  # Construct A(L) (VAR) matrices and hence the C(1) matrix 
  bar <- trimr(bar,nconst,0)
	k   <- ncol(v)
  a   <- array(0, c(k^2,p))
  a1 <- eye(k)

  for (i in seq(p)) {
    a[,i] <- c(bar[(1+k*(i-1)):(i*k),])
    a1 <- a1 - reshapeg(a[,i],k,k)    
  }
  # Generate S matrix from Choleski decomposition   
  s <- t(chol(vc))
  
  # Compute impulse responses based on Choleski 
  # Identification is based on short-run restrictions without sign restrictions    
  impulse_nosign <- irf(a,s,TRUE,k,q,p)     

  cat('\nIRF based on short-run restrictions without sign restrictions\n')
  print(impulse_nosign)
  cat('\n')

  # Compute impulse responses rotated by the matrix Q'    
  th   <- 0.2*pi
  qmat <- matrix(c(cos(th), -sin(th),
                   sin(th),  cos(th)), nrow=2, byrow=T)

  impulse_rotate <- irf(a,s %*% t(qmat),TRUE,k,q,p)

  cat('\nIRF based on based on Choleski which is rotated by the matrix Q\n')
  print(impulse_rotate)
  cat('\n')

  # Select IRFs that satisfy sign restrictions      
  nsim <- 10000         # Number of draws to compute IRFs    
  impulse_select <- zeros(q+1,1)

  nsuccess <- 0.0
  pb <- txtProgressBar(min = 0, max=nsim, style=3)
  for (iter in seq(nsim)) {
    th <- pi*runif(1)
    qmat <- matrix(c(cos(th), -sin(th),
                      sin(th),  cos(th)), nrow=2, byrow=T)
    impulse <- irf(a, s %*% t(qmat),TRUE,k,q,p)      
     
    # Choose impulses that satisfy the sign restrictions for all values
    if (min(impulse[,c(1, 2, 4)]) >= 0.0 
          && min(impulse[,3]) < 0.0 ) {
     impulse_select <- cbind(impulse_select, impulse)
     nsuccess <- nsuccess + 1
    }    
    setTxtProgressBar(pb, iter)     
  }
  close(pb)

	cat('\nNumber of simulations          = ', nsim)
  cat('\nNumber of successful draws (#) = ', 100*nsuccess/nsim)

  # Find the model that corresponds to the median impulses         
  impulse_select <- t( trimr(t(impulse_select),1,0) )
    
  # Choose contemporaneous impulses = nsuccess x cols(y)^2
  impulse_contemp <- reshape(rbind(impulse_select[1,]),ncol(impulse_select)/ncol(y)^2,ncol(y)^2)

  sd_impulse_contemp <- apply(impulse_contemp, 2, sd)
  med_impulse_contemp <- apply(impulse_contemp, 2, median)
  zimpulse_contemp <- (impulse_contemp - repmat(med_impulse_contemp,nrow(impulse_contemp),1))/repmat(sd_impulse_contemp,nrow(impulse_contemp),1)        #     Compute deviations  
   
  # Compute total deviations for each model
  total <- t(rowSums(zimpulse_contemp ))       

  # Choose the set of impulses which yields the minimum absolute total deviation for each model  
  kmin <- which.min(abs(total))

  #   Overall median impulse
  impulse_med <- impulse_select[,((kmin-1)*ncol(y)^2+1):(kmin*ncol(y)^2)]   


  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i", mfrow=c(2,2))
    
  # Plot IRFs
  #--------------------------------------------------------#
	plot(seqa(0,1,q+1),impulse_med[1:(q+1),1],type="l",
       main = 'Supply shock',
       ylab = 'Output',
       xlab = 'Quarter',
       bty = "l")


  #--------------------------------------------------------#
  plot(seqa(0,1,q+1),impulse_med[1:(q+1),2],type="l",
       main = 'Demand shock',
       ylab = 'Output',
       xlab = 'Quarter',
       bty = "l")  

	#--------------------------------------------------------#
  plot(seqa(0,1,q+1),impulse_med[1:(q+1),3],type="l",
       main = 'Supply shock',
       ylab = 'Price',
       xlab = 'Quarter',
       bty = "l")
  #--------------------------------------------------------#
  plot(seqa(0,1,q+1),impulse_med[1:(q+1),4],type="l",
       main = 'Demand shock',
       ylab = 'Price',
       xlab = 'Quarter',
       bty = "l")

}


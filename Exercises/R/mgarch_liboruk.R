#============================================================================
#
#   Program to estimate MGARCH model of UK libor rates based on an
#   SVAR model with GARCH conditional variances on the structural
#   disturbances and identification based on short-run restrictions.
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
# Wrapper function for log-liklihood
#----------------------------------------------------------------------------
neglogl <- function( b,v ) {
  
  lf <- - mean( loglt( b,v ) )
  return(lf)
}

#----------------------------------------------------------------------------
# Log-likelihood function for SVAR
#----------------------------------------------------------------------------
loglt <- function( b,v ) {  
  # Unpack parameter vector
  delta <- b[14:19]^2
  alpha <- pnorm( b[20:25] )
  beta  <- pnorm( b[26:31] )
  
  t <- nrow(v)
  n <- ncol(v)
  lf      <- rep(0, t)
  
  
  binv  <- matrix(c(1, 0, 0, 0, 0, 0,
                    b[1],   1,  b[8],  0, 0, 0,
                    b[2], 0,  1,  0,  0, 0,
                    b[3],  0, 0,  1,  0,  0,
                    b[4],  b[6], b[9], b[11], 1, 0,
                    b[5], b[7], b[10], b[12],  b[13],  1), nrow=6, byrow=T)
  
  
  u <- t ( inv(binv) %*% t(v) ) 


  for (i in seq(t)) {
    if (i == 1)
      d <- delta + alpha*t( (apply(u, 2, sd)^2) )+ beta* t( (apply(u, 2, sd)^2) )
    else    
      d <- delta + alpha*t( (u[i-1,]^2) ) + beta*d    
    omega <- binv %*% diag(as.vector(d)) %*% t(binv)   
    lf[i] <- -0.5*n*log(2*pi) - 0.5*log(det(omega)) - 0.5*v[i,] %*% inv(omega) %*% cbind(v[i,])
  }
  return(lf)
}

#--------------------------------------------------------------------------
#   Return matrix to compute time-varying covariances
#--------------------------------------------------------------------------
getSmatrix <- function( b,v ) {
  # Unpack parameter vector
  delta <- b[14:19]^2
  alpha <- pnorm( b[20:25] )
  beta  <- pnorm( b[26:31] )
  
  t <- nrow(v)
  n <- ncol(v)
  S       <- array(0, c( t,n^2 ))
  
  
  binv  <- matrix(c(1, 0, 0, 0, 0, 0,
                    b[1],   1,  b[8],  0, 0, 0,
                    b[2], 0,  1,  0,  0, 0,
                    b[3],  0, 0,  1,  0,  0,
                    b[4],  b[6], b[9], b[11], 1, 0,
                    b[5], b[7], b[10], b[12],  b[13],  1), nrow=6, byrow=T)
  
  u <- t ( inv(binv) %*% t(v) )
  
  for (i in seq(t)) {
    if (i == 1)          
      d <- delta + alpha*t( (apply(u, 2, sd)^2) ) + beta*t( (apply(u, 2, sd)^2) )
    else            
      d <- delta + alpha*t( (u[i-1,]^2) ) + beta*d
    tmp    <- binv %*% diag( as.vector(sqrt(d)) ) 
    
    S[i,] <- as.vector( tmp )
  }       
  
  return(S)
}

#
# ------------------------ MGARCH Libor --------------------------------------
#
mgarch_liboruk <- function( ) {
  # Load data: 1 January 2004 to 28 May 2008
  # Missing observations in credit risk series fixed by interpoloation
  # 11/06/2004, 9/02/2005 & 26/12/2005
  
  # volatility, liquidity, swap, repo, default, libor
  data <- as.matrix(read.table("libor_data.dat"))
  
  y <- data
  y[,5] <- y[,5]*10000  # default*10000
  
  # Dummy variable to allow for structural break in volatility 2/07/2007  
  dum <- rep(0, nrow(y))
  dum[913:length(dum)] <- 1
  
  # Set SVAR parameters
  p <- 2      # Order of VAR lags
  q <- 21     # Order of VMA lags      
  
  # Estimate the VAR with p lags and a constant     
  ylag   <- rep(1, nrow(y))
  nconst <- 1
  
  for (i in seq(p)) {
    ylag <- cbind(trimr(ylag,1,0), trimr(y,0,i))
  }
  
  # OLS Estimates
  bar    <- lm(trimr(y,p,0) ~ ylag - 1)$coef
  v      <- trimr(y,p,0) - ylag %*% bar
  omegav <- t(v) %*% v/nrow(v)
  s      <- t( chol(omegav) )
  
  
  # Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
  bar <- trimr(bar,nconst,0)
  k   <- ncol(v)
  a   <- array(0, c(k^2,p))
  a1  <- diag(k)
  
  for (i in seq(p)) {
    tmp    <- bar[(1+k*(i-1)):(i*k),]
    a[,i] <- c(tmp)
    a1     <- a1 - reshapeg(a[,i],k,k)
  }
  
  # Invert A(1) matrix needed for the MLE
  a1inv <- inv(a1)           
  
  
  bstart <- c(0.001855486206481,
              0.068693371469487,
              0.023880592273368,
              0.001989809184491,
              0.013648672676092,
              0.000281803828058,
              -0.121999669995451,
              -0.002888991832387,
              0.002650363820455,
              0.217906436715095,
              0.001195294344330,
              0.008613191798914,
              0.175400639477836,
              1.103000246244063,
              0.440791629588056,
              0.495605700916872,
              0.365097267442524,
              -0.011572135693424,
              0.650892298805934,
              -1.176100342250688,
              -0.342100219212525,
              -1.059396487450620,
              -0.545902063123713,
              -0.729900334997059,
              -0.402000678739595,
              1.013901399027622,
              -1.601300249119193,
              1.018606365787122,
              0.791195906556810,
              0.818998254746848,
              0.338395871264844)
  estResults <- optim(bstart, neglogl, v=v, method="BFGS")
  b <- estResults$par
  fval <- estResults$val
  
  cat('\n Log-likelihood function ')
  cat('\n', fval, '\n' )
  print(cbind(b))
  
  # Matrix required to compute time-varying variance-covariances              
  s_store <- getSmatrix( b,v )
  
  # Construct C(L) matrices (vector moving average)  
  c <- diag(k)
  c <- cbind(as.vector(c))
  
  for (i in seq(q)) {
    ss <- array(0, c(k,k))
    j   <- 1.0
    
    while (j <= min( c(p, i))) {    
      ss <- ss + reshapeg(a[,j],k,k) %*% reshapeg(c[, (i-j+1)],k,k)
      j   <- j + 1    
    }
    tmp <- t( ss)
    c  <- cbind(c,  as.vector(tmp))
  }
  
  # Construct orthogonal impulse response functions
  tdate   <- 1          #Choose a point in time 
  s       <- t( reshapeg(s_store[tdate,],k,k) )
  impulse <- c(t(s))
  
  for (i in 2:q) {
    tmp      <- t( reshapeg( c[,i],k,k ) %*% s )
    impulse  <- cbind(impulse, c(tmp))
  }
  impulse <- t(impulse)
  
  # Construct variance decompositions 
  tmp0   <- reshapeg( apply(impulse^2, 2, cumsum ),q*k,k )
  tmp1   <- repmat( rowSums(tmp0),1,k )
  decomp <- reshapeg( 100*tmp0/ tmp1 , q , k^2 )
  cat('\nVariance decomposition (#) over time of libor in terms of factors')
  cat('\n*****************************************************************')
  cat('\n ' )
  cat('\nPeriod         y1         y2        y3        y4         y5       y6\n')   
  
  print( cbind(seq(nrow(decomp)), decomp[,(k*k-(k-1)):(k*k)]) )
  
  # Graph the instantaneous effect of shocks on the libor over time  
  rr            <- nrow(v)
  decomp_time1  <- array(0, c( rr,k ))           # Lag 1 (instantaneous)      
  decomp_time5  <- array(0, c( rr,k ))           # Lag 5       
  decomp_time20 <- array(0, c( rr,k ))           # Lag 20
  
  pb <- txtProgressBar(min=0, max=rr, style=3)
  for (j in seq(rr)) {
    s       <- t (reshapeg(s_store[j,],k,k) )
    impulse <- c(t(s))
    
    for (i in 2:q) {
      tmp      <- t( reshapeg( c[,i],k,k ) %*% s )
      impulse  <- cbind(impulse, as.vector(tmp))
    }
    
    impulse <- t(impulse)
    tmp0    <- reshapeg( apply(impulse^2, 2, cumsum),q*k,k )
    tmp1    <- repmat( rowSums(tmp0,2),1,k )
    decomp  <- reshapeg( 100*tmp0/ tmp1 , q , k^2 )
    
    decomp_time1[j,] <- decomp[1,(k*k-(k-1)):(k*k)]
    decomp_time5[j,] <- decomp[5,(k*k-(k-1)):(k*k)] 
    decomp_time20[j,] <- decomp[20,(k*k-(k-1)):(k*k)]
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(mfrow=c(3,3))
  t <- seq(from=2004, by=1/249, to=2008.61)
  
  #--------------------------------------------------------#
  # Panel (a)
  plot(t,decomp_time1[,1],type='l',
       main = 'Global Risk Factor',
       ylab = 'One Period Lag',
       xlab = 't',
       ylim= c(0,100),
       bty = 'l')
  
  #--------------------------------------------------------#
  # Panel (b)
  plot(t,decomp_time1[,3],type='l',
       main = 'Broad Liquidity Factor',
       ylab = "",
       xlab = 't',
       ylim= c(0,100),
       bty = 'l')
  #--------------------------------------------------------#
  # Panel (c)
  plot(t,decomp_time1[,6],type='l',
       main = 'Idiosyncratic Factor',
       xlab = 't',
       ylim= c(0,100),
       ylab = '',
       bty = 'l')
  
  #--------------------------------------------------------#
  # Panel (d)
  plot(t,decomp_time5[,1],type='l',
       main = '',
       ylab = 'Five Period Lag',
       xlab = 't',
       ylim= c(0,100),
       bty = 'l')
  
  #--------------------------------------------------------#
  # Panel (e)
  plot(t,decomp_time5[,3],type='l',
       main = '',
       xlab = 't',
       ylab = '',
       ylim= c(0,100),
       bty = 'l')
  
  #--------------------------------------------------------#
  # Panel (f)
  plot(t,decomp_time5[,6],type='l',
       main = '',
       ylab = '',
       xlab = 't',
       ylim= c(0,100),
       bty = 'l')
  #--------------------------------------------------------#
  # Panel (g)
  plot(t,decomp_time20[,1],type='l',
       ylab = 'Twenty Period Lag',
       main = '',
       xlab = 't',
       ylim= c(0,100),
       bty = 'l')
  
  #--------------------------------------------------------#
  # Panel (h)
  plot(t,decomp_time20[,3],type='l',
       main = '',
       xlab = 't',
       ylab = '',
       ylim= c(0,100),
       bty = 'l')
  
  #--------------------------------------------------------#
  # Panel (i)
  plot(t,decomp_time20[,6],type='l',
       main = '',
       ylab = 'Price',
       xlab = 't',
       ylim= c(0,100),
       bty = 'l')
  
}


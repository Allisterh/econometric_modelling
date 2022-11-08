#============================================================================
#
#   Program to find tau for Union of rejections unit root tests 
#   in the presence of a trend break
#
#============================================================================
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, inv, recserar
source("EMTSUtil.R")


#----------------------------------------------------------------------------
# Detrending function: 
# cbar=-7 constant 
# cbar=-13.5 linear trend 
# cbar=-T for OLS detrending 
#----------------------------------------------------------------------------
glsdetrend <- function(y,x,cbar)  {  
  t <- length(y) 
  yc <- matrix(c(y[1],(trimr(y,1,0)-(1 + cbar/ t) * trimr(y,0,1))), nrow=t)
  xc <- matrix(c(x[1, ],(trimr(x,1,0)-(1 + cbar/ t) * trimr(x,0,1))), nrow=t )
  b <- lm(yc ~ xc - 1)$coef
  u <- y- x %*% b 
  
  return(u)
} 

#---------------------------------------------------------------------------- 
# ADF coefficient and t tests:u must be already de-trended 
#---------------------------------------------------------------------------- 
adftests <- function(u,k)  {
  
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x <- trimr(u,0,1) 
  
  # Set up lags of du 
  if (k > 0)  { 
    ldu <- lag.matrix(du,seq(from = 1,to = k,by = 1)) 
    x <- trimr(cbind(x,ldu),k,0) 
  } 
  
  xxi <- inv( t(x) %*% x ) 
  b <- xxi %*% t(x) %*% trimr(du,k,0) 
  e <- trimr(du,k,0)-x %*% b 
  s2 <- t(e) %*% e/ length(e)
  
  adfstats <- rbind(length(u)*b[1],b[1]/ sqrt(s2*xxi[1,1]))
  
  return(adfstats)
}

#---------------------------------------------------------------------------- 
# M tests 
#---------------------------------------------------------------------------- 
mtests <- function(u,k)  {  
  s2 <- ar1rvar(u,k) 
  n <- length(u) 
  tmp <- sum(u[1:(n-1)]^2)   
  
  u2 <- tmp/n^2 
  mza <- (u[n]^2/n-s2)/ (2*u2)
  msb <- sqrt(u2/ (s2))
  mzt <- msb %*% mza 
  
  tests <- cbind(mza,msb,mzt) 
  
  return(tests)
}
#---------------------------------------------------------------------------- 
# Autoregressive long run variance estimator 
#---------------------------------------------------------------------------- 
ar1rvar <- function(u,k)  {
  
  du <- cbind(trimr(u,1,0)-trimr(u,0,1) )
  x <- trimr(u,0,1) 
  
  if (k > 0)  { 
    x <- cbind(x,lag.matrix(du,seqa(1,1,k))) 
    x <- cbind(x[-(1:k), ])
  }
  x <- cbind(x)  
  b <- lm(trimr(du,k,0) ~ x - 1)$coef
  
  e <- trimr(du,k,0)- x %*% b 
  s2 <- t(e) %*% e/ length(e)  
  if (k > 0)  { 
    s2 <- s2/ (1-sum(trimr(b,1,0)))^2
  }   
  return(s2)
}
          
unit_urbreak <- function() {
 
  t    <- 1000
  reps <- 100000
  
  
  lambdav <- seqa(0.15,0.05,15)
  cbar    <- c(-17.6,-17.8,-18.2,-18.4,-18.6,-18.4,-18.4,-18.2,-18.0, -17.6,-17.4,-17.0,-16.6,-16.0,-15.2)
  dfolscv <- c(-3.58,-3.67,-3.73,-3.77,-3.81,-3.84,-3.86,-3.87,-3.87, -3.88,-3.87,-3.85,-3.83,-3.79,-3.74)
  dfglscv <- c(-3.38,-3.42,-3.43,-3.44,-3.45,-3.45,-3.45,-3.44,-3.43, -3.42,-3.40,-3.37,-3.33,-3.28,-3.23)
  mztols  <- c(-3.42,-3.48,-3.53,-3.56,-3.58,-3.60,-3.61,-3.62,-3.62, -3.61,-3.60,-3.58,-3.55,-3.51,-3.46)
  
  # Initialise arrays
  dfols <- array(0, c(reps,length(lambdav)))
  dfgls <- dfols 
  mztols <- dfols
  gdf <- rep(1, length(lambdav))
  gm <- gdf
  
  
  for (lambdac in 1:length(lambdav)) {
   set.seed(1234)
    
    TB <- floor(lambdav[lambdac]*t)
    x  <- cbind(rep(1,t), seqa(1,1,t), cbind(c(rep(0, TB), seqa(1,1,t-TB))))


   for (rep in seq(reps)) {
     y    <- cumsum(rnorm(t))
     uols <- glsdetrend(y,x,-t)
     ugls <- glsdetrend(y,x,cbar[lambdac])
     
     dfols[rep,lambdac]  <- trimr(adftests(uols,0),1,0)
     dfgls[rep,lambdac]  <- trimr(adftests(ugls,0),1,0)
     mztols[rep,lambdac] <- trimr(t(mtests(uols,0)),2,0)        
   }
   size <- mean((dfols[,lambdac] < (gdf[lambdac]*dfolscv[lambdac])) | (dfgls[,lambdac] < (gdf[lambdac]*dfglscv[lambdac])))
   
   while (size >= 0.05) {
     gdf[lambdac] <- gdf[lambdac] + 0.00001
     
     size <- mean((dfols[,lambdac] < (gdf[lambdac]*dfolscv[lambdac])) | (dfgls[,lambdac] < (gdf[lambdac]*dfglscv[lambdac])))     
   }
   cat('\n', size)
   
   size <- mean((mztols[,lambdac] < (gm[lambdac]*dfolscv[lambdac])) | (dfgls[,lambdac] < (gm[lambdac]*dfglscv[lambdac])))
   while (size >= 0.05) {
     gm[lambdac] <- gm[lambdac] + 0.00001
     
     size <- mean((mztols[,lambdac] < (gm[lambdac]*dfolscv[lambdac])) | (dfgls[,lambdac] < (gm[lambdac]*dfglscv[lambdac])))     
   }
   cat('\n', size)
  }
  
  print(cbind(lambdav, gdf, gm))
    
}


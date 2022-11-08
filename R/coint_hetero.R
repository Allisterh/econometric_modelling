#============================================================================
#
#   Program to investigate the effects of heteroskedasticity 
#   on the size of the trace test of cointegration
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - trimr
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Variance break DGP 
#----------------------------------------------------------------------------
varbreakdgp <- function(e,frac,newvar) {
  t <- nrow(e)
  n <- ncol(e)  
  tb    <- floor(t*frac)
  sig2  <- matrix( c(rep(1, tb*n), newvar*rep(1, (t-tb)*n)), nrow=t, byrow=T )
  y     <- e*sig2  
  return(y)
}

#----------------------------------------------------------------------------
#  GARCH(1,1) DGP 
#----------------------------------------------------------------------------
garchdgp <- function(e,phi1,phi2) {
  t <- nrow(e)
  n <- ncol(e)  
  v2    <-  e
  sig2  <-  matrix(1, t,n)
  for (j in 2:t) {
    sig2[j,] <- (1-phi1-phi2) + phi1*v2[j-1,]^2 + phi2*sig2[j-1,]
    v2[j,]   <- e[j,]*sqrt(sig2[j,])    
  }  
  y <- v2  
  return(y)
}


#----------------------------------------------------------------------------
#  Trace test
#----------------------------------------------------------------------------
tracetest <- function(y,p,model) {
  dy <- trimr(y,1,0)-trimr(y,0,1) 
  z0 <- trimr(dy,p-1,0)
  z1 <- trimr(y,p-1,1)
  z2 <- c()
#   for (j in seq(p-1)) {
#     z2 <- cbind(z2, trimr(dy,p-1-j,j))
#   }
    
  if (model == 2)
    z1 <- cbind(trimr(y,p-1,1), rep(1, length(y)-p,1) )
  
  else if (model == 3 ) 
    z2 <- cbind(rep(1, length(y)-p), z2)  
  else if (model == 4) {
    z1 <- cbind(z1, seqa(1,1,length(y)-p))
    z2 <- cbind(rep(1, length(y)-p), z2)
  } else if (model == 5 )
    z2 <- cbind(rep(1, length(y)-p, seqa(1,1,rows(y)-p), z2))  
  
  if (is.null(ncol(z2))) {
    r0 <- z0 
    r1 <- z1
  }else {
    r0 <- lm(z0 ~ z2 -1)$residuals
    r1 <- lm(z1 ~ z2 - 1)$residuals    
  }   
  C     <- inv(chol(t(r1) %*% r1))
  lambda <- eigen(t(C) %*% t(r1) %*% r0 %*% inv( t(r0) %*% r0) %*% t(r0) %*% r1 %*% C)$values
  tr <- -nrow(r0) %*% sum(log(1-flipud(lambda)))
  return(tr)
}


# 
#------ Heteroskedasticity and the Distribution of the Trace Statistic ------
#            
coint_hetero <- function() {
  reps  <- 10000
  n     <- 4
  model <- 1
  Tv    <- c(100,200,400,800,1600,3200)
  cv1   <- c(4.173,12.285,24.102,39.921,59.829,83.428)
  
  tr <- matrix(0, reps,6)
  
  for (tc in seq(Tv)) {
    t <- Tv[tc]
    pb <- txtProgressBar(min=0, max=reps, style=3)
    for (rep in seq(reps)) {
      e <- matrix(rnorm(t*n), t, n)
      
      # iid 
      tr[rep,1] <- tracetest(apply(e, 2, cumsum),1,model)
      
      # GARCH(1,1)
      tr[rep,2] <- tracetest(apply(garchdgp(e,0.3,0.6), 2, cumsum),1,model)
      tr[rep,3] <- tracetest(apply(garchdgp(e,0.3,0.69), 2, cumsum),1,model)
      
      # Variance break 
      tr[rep,4] <- tracetest(apply(varbreakdgp(e,0.1,2), 2, cumsum),1,model)
      tr[rep,5] <- tracetest(apply(varbreakdgp(e,0.5,2), 2, cumsum),1,model)
      tr[rep,6] <- tracetest(apply(varbreakdgp(e,0.9,2), 2, cumsum),1,model) 
      setTxtProgressBar(pb, rep)
    }  
    close(pb)
    cat('\n       T   iid GARCHa GARCHb Break-a Break-b Break-c\n')    
    print(matrix(c(t, colMeans(tr>cv1[n])), nrow=1))
  }
}



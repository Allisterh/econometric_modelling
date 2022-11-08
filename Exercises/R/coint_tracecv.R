#============================================================================
#
#   Program to approximate asymptotic critical values for the trace test
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

model <- 5                  # Choose model (1 - 5)                  
t     <- 1000                   
nreps <- 100000             
tr    <- rep(0, nreps)

pv <- c(0.90, 0.95,  0.99)

q <- c()


for (nr in seq(6))                #  Up to 6 common trends 
{
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  for (j in seq(nreps)) {
    db <- matrix(rnorm(t*nr), nrow=t, ncol=nr) 
    f  <- as.matrix(apply(db, 2, cumsum))
    if (model == 2)     
      f <- cbind(f, rep(1, t))
    
    else if (model == 3 ) {  
      f[,1] <- seqa(1,1,t)
      f    <- f - colMeans(f)
    } else if (model == 4) {
      f <- cbind(f,  seqa(1,1,t))
      f <- f - colMeans(f)
    }else if (model == 5) {
      f[,1] <- (seqa(1,1,t)^2)      
      x      <- cbind(rep(1, t), seqa(1,1,t))
      f      <- lm(f ~ x - 1)$residuals
    }
    setTxtProgressBar(pb, j)
    #  Note that db is computed as a forward difference relative to f     
    m1 <- t(trimr(f,0,1)) %*% trimr(db,1,0)/t      
    m2 <- t(trimr(f,0,1)) %*% trimr(f,0,1)/(t^2)

    tr[j] <- sum( diag( t(m1) %*% inv(m2) %*% m1 ) )
  }
  close(pb)  
  q <- cbind(q,  quantile(tr,pv))  
}
cat('\nNumber of common trends   = ',model)
cat('\n------------------------------------------------------------------')
cat('\n      n-r       1         2       3         4        5       6\n')   
print( unname(cbind(pv, q) ))
 
       
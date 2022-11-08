#============================================================================
#
#   Program to estimate a bivariate gaussian copula for asset returns
#   Asset price data from 6 August 2010 to 2 January 2001 (note that the
#   data are in reverse order ie from recent to past)
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Copula log-likelihood function
#----------------------------------------------------------------------------
neglog <- function (b,y) {
  m1 <- b[1]                 # Asset 1
  s1 <- abs(b[3])
  z1 <- (y[,1] - m1)/s1
  f1 <- dnorm(z1)/s1
  
  m2 <- b[2]                # Asset 2 
  s2 <- abs(b[4])
  z2 <- (y[,2] - m2)/s2
  f2 <- dnorm(z2)/s2
  
  r  <- b[5]                 # Dependence
  
  lt <- -log(2*pi) - log(s1*s2) - 0.5*log(1 - r^2) - 0.5*(z1^2 - 2*r*z1*z2 + z2^2)/(1 - r^2) + log(f1) + log(f2)
  
  lf <- -mean(lt)
  return(lf)
}

# 
#---------------------------- Copula Model ----------------------------------
#
# Load data
load("diversify.RData")

t <- 2413

# Select appropriate sample  
pt_apple     <- as.matrix(pt_apple)[1:t]
pt_ford      <- as.matrix(pt_ford)[1:t]

# Compute percentage returns  
r_apple <- 100*diff(log(pt_apple)) 
r_ford  <- 100*diff(log(pt_ford))

y <- cbind(r_apple, r_ford)
t <-nrow(y)

# Compute statistics
m <- colMeans(y)
s <- apply(y, 2, sd)
co <- cor(y)
r <- co[1,2]

# Estimate parameters of the copula 
start <- c(m,s,r)
estResults <- optim(start, neglog, y=y, method="BFGS", hessian=T)
bhat <- estResults$par
hess <- estResults$hess

vc <- (1/t)*inv(hess)

# Wald test of independence        
wd <- (bhat[5] - 0)^2/vc[5,5]
cat('\nWald statistic     = ',wd)
cat('\nP-value            = ',1-pchisq(wd,1))




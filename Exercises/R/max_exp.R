#========================================================================
#
#      Program to estimate a exponential model 
#
#========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Wrapper function to calculate inverse of a given matrix
#-------------------------------------------------------------------------
inv <- function (M) {
  return(solve(M))
}

#
#--------------------------- Algorithms -----------------------------------
# 

# Data
y <- c(3.5, 1.0, 1.5)
t <- length(y)  
niter <- 10 

# Newton-Raphson algorithm

cat('\nNewton Raphson algorithm\n\n')

theta <- 1.0
            
rnames <- seq(niter)
cnames <- c("G", "J", "logL", "th(k)")
res <- array( dim = c( length(rnames),length(cnames)), 
              dimnames = list (rnames, cnames))
             
for (k in seq(niter)) {

  g     <- -1.0/theta + mean(y)/theta^2  
  h     <- 1.0/theta^2 - 2*mean(y)/theta^3
  lnl   <- -log(theta) - mean(y)/theta
  
  theta <- theta - inv(h) %*% g
 
  # populate results table with values  
  res[k, "G"] <- g
  res[k, "J"] <- t(g) %*% g
  res[k, "logL"] <- lnl
  res[k, "th(k)"] <- theta    
}
print(round(res, 6))


# Scoring algorithm 
cat('\nScoring algorithm\n\n')
theta <- 1.0
res <- array( dim = c( length(rnames),length(cnames)), 
              dimnames = list (rnames, cnames))
for (k in seq(niter)) { 
  g     <- -1.0/theta + mean(y)/theta^2  
  i     <- 1.0/theta^2
  lnl   <- -log(theta) - mean(y)/theta
  
  theta <- theta + inv(i) %*% g
 
  # populate results table with values  
  res[k, "G"] <- g
  res[k, "J"] <- t(g) %*% g
  res[k, "logL"] <- lnl
  res[k, "th(k)"] <- theta
}
print(round(res, 6))

# BHHH algorithm 
cat('\nBHHH algorithm\n')

theta <- 1.0
for (k in seq(niter)) {
  
  gt    <- -1.0/theta + y/theta^2
  g     <- mean(gt)
  j     <- t(gt) %*% gt/t
  lnl   <- -log(theta) - mean(y)/theta
  theta <- theta + inv(j) %*% g
  
  # populate results table with values  
  res[k, "G"] <- g
  res[k, "J"] <- j
  res[k, "logL"] <- lnl
  res[k, "th(k)"] <- theta  
}
print(round(res, 6))


# BHHH algorithm with squeezing  

cat('\nBHHH algorithm with squeezing\n')

theta <- 1.0  
rnames <- seq(niter)
cnames <- c("G", "J", "logL", "th(k)", "lambda")
res <- array( dim = c( length(rnames),length(cnames)), 
              dimnames = list (rnames, cnames))

for (k in seq(niter)) {
  gt <- cbind(-1.0/theta + y/theta^2)  
  g  <- mean(gt)
  j  <- t(gt) %*% gt/t
  lnl <- -log(theta) - mean(y)/theta
  
  thetaold <- theta
  lnlold   <- lnl
  lambda <- 1

  theta <- theta + inv(j) %*% g
  lnl   <- -log(theta) - mean(y)/theta
  
  if (lnl > lnlold )
     cat('\nFull iteration step successful at iteration k =' , k)
  else {    
     cat('\nFull iteration step not successful at iteration k =' , k)
  
    for (m in 2:10) {
    lambda <- 1/m
    cat('\n', rep('->', m-1), 'Squeezing with step = ', lambda)
    
    theta <- thetaold + lambda*inv(j)*g 
    lnl   <- -log(theta) - mean(y)/theta
    if (lnl > lnlold)  
      break 
    }    
  }
  
  # populate results table with values  
  res[k, "G"] <- g
  res[k, "J"] <- j
  res[k, "logL"] <- lnl
  res[k, "th(k)"] <- theta
  res[k, "lambda"] <- lambda
}
cat('\n\n')
print(round(res, 6))

# BFGS algorithm   
cat('\nBFGS algorithm\n')

theta <- 1.5
h     <- -1
niter <- 9

rnames <- seq(niter)
cnames <- c("G", "H", "logL", "th(k)")
res <- array( dim = c( length(rnames),length(cnames)), 
              dimnames = list (rnames, cnames))

for (k in seq(niter)) {
  g     <- -1.0/theta + mean(y)/theta^2
  lnl   <- -log(theta) - mean(y)/theta
  theta <- theta - inv(h) %*% g
  
  # populate results table with values  
  res[k, "G"] <- g
  res[k, "H"] <- h
  res[k, "logL"] <- lnl
  res[k, "th(k)"] <- theta
  
  dtheta <- -inv(h) %*% g
  gold   <- g
  g      <- -1.0/theta + mean(y)/theta^2
  dg     <- g - gold
  h      <- h - ( h %*% dtheta %*% t(dg) + dg %*% t(dtheta) %*% h)/ (t(dg) %*% dtheta) + ( 1 + (t(dtheta) %*% h %*% dtheta)/( t(dg) %*% dtheta) ) %*% (dg %*% t(dg)) / (t(dg) %*% dtheta)
}
print(round(res, 6))


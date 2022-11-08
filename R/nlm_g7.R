#============================================================================ 
# 
# Program to estimate a bivariate threshold model 
# 
# Anderson,Anthanasopoulos Vahid(2007,22,63–87). 
# 
# Uses quarterly data on percentage growth rate in real gdp and 
# interest rate spread(10yr-3mth). Source:JAE Archive. 
# 
#============================================================================ 
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, figure
source("EMTSUtil.R")

#---------------------------------------------------------------------------- 
# Gets the data for a specfic country 
#---------------------------------------------------------------------------- 

getData <- function(flag)  {
  load("G7Data.RData")
  if (flag == 1 )  {       # Canada
    y <- as.matrix(canada)      
  } else if (flag == 2 ) { # France
    y <- as.matrix(france) 
    
  } else if (flag == 3 ) { # Germany
    y <- as.matrix(germany)         
  } else if (flag == 4) {  # Italy
    y <- as.matrix(italy)  
  } else if (flag == 5 ) { # Japan
    y <- as.matrix(japan) 
  } else if (flag == 6 )  { # UK
    y <- as.matrix(uk)
  } else { # US
    y <- as.matrix(us)
  }
  return(y)
  
} 

#---------------------------------------------------------------------------- 
# Compute the LM statistic to test a bivariate TAR 
# model assuming one lag in the auxiliary model 
#----------------------------------------------------------------------------

btar_test <- function(yvar,dep,trans,itype)  {
  
  # itype=1:regress yt on[constant ylag,ylag*ylag] 
  # itype=2:regress yt on[constant ylag,ylag*ylag,ylag*ylag^2] 
  # itype=3:regress yt on[constant ylag,ylag*ylag,ylag*ylag^2,ylag*ylag^3] 
  # itype=4:regress yt on[constant ylag,ylag*ylag,ylag*ylag^3] 
  
  # First stage regression 
  y <- trimr(yvar[ ,dep],1,0) 
  x <- cbind(matrix(1, nrow = length(y), ncol = 1),trimr(yvar,0,1)) 
  k <- dim(x)[2] 
  u <- lm(y ~ x - 1)$residuals
  
  
  # Second stage regression  
  z <- trimr(yvar[ ,trans],0,1) # Choose the transition variable(lagged one period) 
  
  if (itype == 1)  { 
    tp <- x[, 2:k] * z
    x <- cbind(x, tp)
  } 
  
  if (itype == 2)  { 
    tp <- x[, 2:k] * z
    tt <- x[, 2:k] * z^2
    x <- cbind(x, tp, tt)    
  } 
  
  if (itype == 3)  { 
    tp <- x[, 2:k] * z
    tt <- x[, 2:k] * z^2
    tr <- x[, 2:k] * z^3
    x <- cbind(x, tp, tt, tr)    
  } 
  
  if (itype == 4)  { 
    tp <- x[, 2:k] * z
    tt <- x[, 2:k] * z^3
    x <- cbind(x, tp, tt)
  } 
  
  e <- lm(u ~ x - 1)$residuals
  
  # Compute lm statistic   
  t <- length(y)
  r2 <- 1-sum(e^2)/sum((u-mean(u))^2)
  lm <- t*r2 
  
  dof <- dim(x)[2]-k   
  return(list(lm = lm,dof = dof))
} 


#-------------------------------------------------------------------------
# Unrestricted likelihood function 
#-------------------------------------------------------------------------
neglog <- function(b,nobs,y1,y2,maxlag){
  
  v1 <- rep(0, nobs)
  v2 <- rep(0, nobs)
  lt <- rep(0, nobs)
  
  
  for (t in (maxlag+1):nobs) {
    w1 <- 1/( 1 + exp(-100*(y1[t-2] - b[11])  ) )        
    v1[t] <- y1[t] - ( b[1]  + b[2]*y1[t-1]  + b[3]*y1[t-2]  + b[4]*y2[t-1]  + b[5]*y2[t-2] ) - ( b[6] + b[7]*y1[t-1] + b[8]*y1[t-2] + b[9]*y2[t-1] + b[10]*y2[t-2] )*w1
    
    w2 <- 1/( 1 + exp(-100*(y1[t-1] - b[22]) ) )         
    v2[t] <- y2[t] - ( b[12] + b[13]*y1[t-1] + b[14]*y1[t-2] + b[15]*y2[t-1] + b[16]*y2[t-2] ) - ( b[17] + b[18]*y1[t-1] + b[19]*y1[t-2] + b[20]*y2[t-1] + b[21]*y2[t-2] )*w2
  }
  
  #     Exclude the first maxlag observations  
  v     <- trimr(cbind(v1,v2),maxlag,0) 
  n <- nrow(v)
  k <- ncol(v)
  
  
  vc <- t(v) %*% v/n
  
  for (t in 1:(nobs-maxlag)) {
    lt[t] <- - k*0.5*log(2*pi) - 0.5*log(det(vc)) - 0.5*v[t,] %*% inv(vc) %*% cbind(v[t,])
  }
  
  
  #  Log-likelihood is not defined for the last maxlag observations  
  lf <- - mean( trimr(lt,0,maxlag) )
  return(lf)
}



#---------------------------------------------------------------------------- 
#   Simulate the model 
#----------------------------------------------------------------------------
btarsim <- function(b,ystart,nobs,maxlag) {
  y1 <- rep(0, nobs) + ystart[1]     
  y2 <- rep(0, nobs) + ystart[2]
  
  for (t in (maxlag+1):nobs) {
    w1 <- 1/( 1 + exp(-100*(y1[t-2] - b[11])  ) )        
    y1[t] <-  ( b[1]  + b[2]*y1[t-1]  + b[3]*y1[t-2]  + b[4]*y2[t-1]  + b[5]*y2[t-2] ) - ( b[6] + b[7]*y1[t-1] + b[8]*y1[t-2] + b[9]*y2[t-1] + b[10]*y2[t-2] )*w1
    
    w2 <- 1/( 1 + exp(-100*(y1[t-1] - b[22]) ) )         
    y2[t] <-  ( b[12] + b[13]*y1[t-1] + b[14]*y1[t-2] + b[15]*y2[t-1] + b[16]*y2[t-2] ) - ( b[17] + b[18]*y1[t-1] + b[19]*y1[t-2] + b[20]*y2[t-1] + b[21]*y2[t-2] )*w2
  }
  #     Exclude the first maxlag observations  
  y     <- trimr(cbind(y1, y2),maxlag,0)
  return(y)
}


#
# ------------------ Bivariate LSTAR Model of G7 Countriess -----------------
#
nlm_g7 <- function()  {
# Set country 
country <- 7    
# 1 for Canada 
# 2 for France 
# 3 for Germany 
# 4 for Italy 
# 5 for Japan 
# 6 for UK 
# 7 for US 


# Get the data 
y <- getData(country) 
y1 <- cbind(y[ ,1]) # Growth rate of GDP 
y2 <- cbind(y[ ,2]) # Spread 
nobs <- length(y1) 

# Perform bivariate LM nonlinearity tests on y1 
dep <- 1 
trans <- 2 
b.test <- btar_test(y,dep,trans,1)
lm1 <- b.test$lm
dof1 <- b.test$dof

b.test <- btar_test(y,dep,trans,2)
lm2 <- b.test$lm
dof2 <- b.test$dof

b.test <- btar_test(y,dep,trans,3)
lm3 <- b.test$lm
dof3 <- b.test$dof

b.test <- btar_test(y,dep,trans,4)
lm4 <- b.test$lm
dof4 <- b.test$dof


cat('\nLM tests of nonlinearity in GDP') 
cat('\n **************') 

cat('\n LM statistic[Test 1] = ',lm1) 
cat('\n  p-value = ',1-pchisq(lm1,dof1)) 

cat('\n LM statistic[Test 2] = ',lm2) 
cat('\n  p-value = ',1-pchisq(lm2,dof2)) 

cat('\n LM statistic[Test 3] = ',lm3) 
cat('\n  p-value = ',1-pchisq(lm3,dof3)) 

cat('\n LM statistic[Test 4] = ',lm4) 
cat('\n  p-value = ',1-pchisq(lm4,dof4)) 

cat('\n ') 
cat('\n ') 

# Perform bivariate LM nonlinearity tests on y2 
dep <- 2 
trans <- 2 
b.test <- btar_test(y,dep,trans,1)
lm1 <- b.test$lm
dof1 <- b.test$dof

b.test <- btar_test(y,dep,trans,2)
lm2 <- b.test$lm
dof2 <- b.test$dof

b.test <- btar_test(y,dep,trans,3)
lm3 <- b.test$lm
dof3 <- b.test$dof

b.test <- btar_test(y,dep,trans,4)
lm4 <- b.test$lm
dof4 <- b.test$dof 

cat('\nLM tests of nonlinearity in the Spread') 
cat('\n **************') 

cat('\n LM statistic[Test 1] = ',lm1) 
cat('\n  p-value = ',1-pchisq(lm1,dof1)) 

cat('\n LM statistic[Test 2] = ',lm2) 
cat('\n  p-value = ',1-pchisq(lm2,dof2)) 

cat('\n LM statistic[Test 3] = ',lm3) 
cat('\n  p-value =',1-pchisq(lm3,dof3)) 

cat('\n LM statistic[Test 4] = ',lm4) 
cat('\n  p-value =  ',1-pchisq(lm4,dof4)) 

# Estimate the bivariate model using BFGS  
final  <- matrix(0, 22,7)
maxlag <- 2 # Maximum number of lags in bivariate model 
start  <- 0.1*rep(1,22)   #     Starting parameter estimates        


estResults <- optim(start, neglog, nobs=nobs, y1=y1, y2=y2, maxlag=maxlag, method="BFGS", hessian=T)
bhat <- estResults$par
lf <- estResults$val
hess <- estResults$hess

t <- nobs
lf    <- -lf
omega <- (1/t)*inv(hess)

cat('\n ') 
print( country )
cat('\nLog-likelihood function= ', lf)
cat('\nCovariance matrix')
print(omega)

# For plotting and simulation use GAUSS estimates
final <- as.matrix(read.table('final_estimates.dat'))

# Simulate all models 
ystart <- c(1,1)   # Start values      

# Canada
tmp <- btarsim(final[,1],ystart,10000,maxlag)
y1s_canada <- tmp[,1]
y2s_canada <- tmp[,2]

# France
tmp <- btarsim(final[,2],ystart,10000,maxlag)
y1s_france <- tmp[,1]
y2s_france <- tmp[,2]

# Germany
tmp <- btarsim(final[,3],ystart,10000,maxlag)
y1s_germany <- tmp[,1]
y2s_germany <- tmp[,2]

# Italy
tmp <- btarsim(final[,4],ystart,10000,maxlag)
y1s_italy <- tmp[,1]
y2s_italy <- tmp[,2]

# Japan
tmp <- btarsim(final[,5],ystart,10000,maxlag)
y1s_japan <- tmp[,1]
y2s_japan <- tmp[,2]

# UK
tmp <- btarsim(final[,6],ystart,10000,maxlag)
y1s_uk <- tmp[,1]
y2s_uk <- tmp[,2]

# US
tmp <- btarsim(final[,7],ystart,10000,maxlag)
y1s_us <- tmp[,1]
y2s_us <- tmp[,2]

#**********************************************************************
#***
#***     Generate graph
#***
#**********************************************************************

figure()
par(mfrow = c(4,3))

plot(y1s_canada[1:100],type='l', main = '', xlab = 't', ylab = 'Canada-output')

plot(y2s_canada[1:100], type='l', main = '', xlab = 't', ylab = 'Canada-spread')

plot(y1s_canada,y2s_canada, type='l', main='', xlab='Canada-output', ylab='Canada-spread')

plot(y1s_france[1:100], xlab='t', ylab='France-output ', main='', type='l')

plot(y2s_france[1:100], main='', type='l', xlab='t', ylab='France-spread')


plot(y1s_france,y2s_france, main='', type='l', xlab='France-output', ylab='France-spread')

plot(y1s_germany[1:100], main='', type='l', xlab='t', ylab='Germany-output')

plot(y2s_germany[1:100], main='', type='l', xlab = 't', ylab='Germany-spread')

plot(y1s_germany,y2s_germany, main='', type='l', xlab='Germany-output', ylab='Germany-output')

plot(y1s_italy[1:100], main='', type='l', xlab='t', ylab='Italy-output')

plot(y2s_italy[1:100], main='', type='l', xlab='t', ylab='Italy-spread')

plot(y1s_italy,y2s_italy, main='', type='l', xlab='Italy-output', ylab='Italy-spread')

figure()
par(mfrow=c(3,3))

plot(y1s_japan[1:100], main='', type='l', xlab='t', ylab='Japan-output')

plot(y2s_japan[1:100], main='', type='l', xlab='t', ylab='Japan-spread')

plot(y1s_japan,y2s_japan, main='', type='l', xlab='Japan-output', ylab='Japan-spread')

plot(y1s_uk[1:100], main='', type='l', xlab='t', ylab='UK-output')

plot(y2s_uk[1:100], main='', type='l', xlab='t', ylab='UK-spread')

plot(y1s_uk,y2s_uk, main='', type='l', 
     xlab='UK-output', ylab='UK-spread',
     xlim=c(-8e300, 17e300), ylim=c(0e300, -5e300))

plot(y1s_us[1:100], main='', type='l', xlab='t', ylab='US-output')

plot(y2s_us[1:100], main='',type='l', xlab='t', ylab='US-spread')

plot(y1s_us,y2s_us, main='', type='l', xlab='US-output', ylab='Us-spread')
} 


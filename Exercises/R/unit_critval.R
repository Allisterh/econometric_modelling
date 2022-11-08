#============================================================================ 
# 
# Simulate critical values for DF,M and P tests 
# 
#============================================================================ 
rm(list = ls(all=T))
graphics.off()

#
# ------------------------ Helper Functions ----------------------------------
#
# Load required functions - trimr, inv
source("EMTSUtil.R")

#------------------------------------------------------------------------- --
# Detrending function: 
# cbar=-7 constant; 
# cbar=-13.5 linear trend 
# cbar=-T for OLS detrending 
#----------------------------------------------------------------------------
glsdetrend <- function(y,x,cbar)  {  
  t <- length(y) 
  yc <- matrix(c(y[1],(trimr(y,1,0)-(1 + cbar/ t) * trimr(y,0,1))), nrow=t)
  xc <- matrix(c(x[1, ],(trimr(x,1,0)-(1 + cbar/ t) * trimr(x,0,1))), nrow=t )
  b <- lm(yc ~ xc - 1)$coef
  u <- y- x %*% b 
  
  return(list(u = u,b = b))
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

#---------------------------------------------------------------------------- 
# P-test 
#---------------------------------------------------------------------------- 

Ptest <- function(y,x,cbar,k)  {
  n <- length(y) 
  uc <- glsdetrend(y,x,cbar)$u 
  u0 <- glsdetrend(y,x,0)$u 
  s2 <- ar1rvar(uc,k) 
  
  uc <- matrix(c(uc[1],trimr(uc,1,0)-(1 + cbar/n) * trimr(uc,0,1)), nrow=n)
  u0 <- matrix(c(u0[1],trimr(u0,1,0)-trimr(u0,0,1)), nrow=n)
  
  pt <- (t(uc) %*% uc-(1 + cbar/n) %*% t(u0) %*% u0)/ s2  
  return(pt)
} 


#
# --------------- Critical Values for unit root tests -----------------------
#
unit_critval <- function()  { 
    
    set.seed(1234) 
 
    reps <- 100000
 
    Tv <- cbind(25,50,100,250,500,1000) 
    pv <- cbind(0.01,0.05,0.10) 
 
    adf1 <- matrix(0, nrow = length(Tv), ncol = length(pv)) 
    adf2 <- adf1 
    adf1cols <- adf1 
    adf2cols <- adf1 
    adf1cgls <- adf1 
    adf2cgls <- adf1 
    adf1tols <- adf1 
    adf2tols <- adf1 
    adf1tgls <- adf1 
    adf2tgls <- adf1 
 
    MSB <- adf1 
    MSBcols <- adf1 
    MSBcgls <- adf1 
    MSBtols <- adf1 
    MSBtgls <- adf1 
    MZtcols <- adf1 
    MZtcgls <- adf1 
    MZttols <- adf1 
    MZttgls <- adf1 
 
    Pc <- adf1 
    Pt <- adf1 
 
    adf <- matrix(0, nrow = reps, ncol = 2) 
    M <- matrix(0, nrow = reps, ncol = 3) 
 
    adfcols <- matrix(0, nrow = reps, ncol = 2) 
    Mcols <- matrix(0, nrow = reps, ncol = 3) 
    adfcgls <- matrix(0, nrow = reps, ncol = 2) 
    Mcgls <- matrix(0, nrow = reps, ncol = 3) 
    adftols <- matrix(0, nrow = reps, ncol = 2) 
    Mtols <- matrix(0, nrow = reps, ncol = 3) 
    adftgls <- matrix(0, nrow = reps, ncol = 2) 
    Mtgls <- matrix(0, nrow = reps, ncol = 3) 
    Pcgls <- matrix(0, nrow = reps, ncol = 1) 
    Ptgls <- matrix(0, nrow = reps, ncol = 1) 
 
    for (Tc in (1:(length(Tv)))) { 
 
        t <- Tv[Tc] 
 
        for (rep in (1:reps)) { 
 
            y <- cumsum(rnorm(t)) 
            ycols <- glsdetrend(y,matrix(1, nrow=t, ncol=1),-t)$u 
            xxc <- cbind(rep(1,t), seqa(1,1,t))
            ytols <- glsdetrend(y,xxc,-t)$u 
            ycgls <- glsdetrend(y,matrix(1, nrow=t, ncol=1),-7)$u 
            ytgls <- glsdetrend(y,xxc,-13.5)$u 
            
 
            adf[rep, ] <- t(adftests(y,0)) 
           
            
            M[rep, ] <- t(mtests(y,0)) 
            adfcols[rep, ] <- t(adftests(ycols,0)) 
            Mcols[rep, ] <- t(mtests(ycols,0)) 
            adfcgls[rep, ] <- t(adftests(ycgls,0)) 
            Mcgls[rep, ] <- t(mtests(ycgls,0)) 
            adftols[rep, ] <- t(adftests(ytols,0)) 
            Mtols[rep, ] <- t(mtests(ytols,0)) 
            adftgls[rep, ] <- t(adftests(ytgls,0)) 
            Mtgls[rep, ] <- t(mtests(ytgls,0)) 
            Pcgls[rep] <- Ptest(y,matrix(1, nrow = t, ncol = 1),-7,0) 
            Ptgls[rep] <- Ptest(y,cbind(rep(1,t),seqa(1,1,t)),-13.5,0)  
        } 
        
 
        adf1[Tc, ] <- t(quantile(cbind(adf[ ,1]),pv)) 
        adf2[Tc, ] <- t(quantile(cbind(adf[ ,2]),pv)) 
        adf1cols[Tc, ] <- t(quantile(cbind(adfcols[ ,1]),pv)) 
        adf2cols[Tc, ] <- t(quantile(cbind(adfcols[ ,2]),pv)) 
        adf1cgls[Tc, ] <- t(quantile(cbind(adfcgls[ ,1]),pv)) 
        adf2cgls[Tc, ] <- t(quantile(cbind(adfcgls[ ,2]),pv)) 
        adf1tols[Tc, ] <- t(quantile(cbind(adftols[ ,1]),pv)) 
        adf2tols[Tc, ] <- t(quantile(cbind(adftols[ ,2]),pv)) 
        adf1tgls[Tc, ] <- t(quantile(cbind(adftgls[ ,1]),pv)) 
        adf2tgls[Tc, ] <- t(quantile(cbind(adftgls[ ,2]),pv)) 
        MSB[Tc, ] <- t(quantile(cbind(M[ ,2]),pv)) 
        MSBcols[Tc, ] <- t(quantile(cbind(Mcols[ ,2]),pv)) 
        MSBcgls[Tc, ] <- t(quantile(cbind(Mcgls[ ,2]),pv)) 
        MSBtols[Tc, ] <- t(quantile(cbind(Mtols[ ,2]),pv)) 
        MSBtgls[Tc, ] <- t(quantile(cbind(Mtgls[ ,2]),pv)) 
        MZtcols[Tc, ] <- t(quantile(cbind(Mcols[ ,3]),pv)) 
        MZtcgls[Tc, ] <- t(quantile(cbind(Mcgls[ ,3]),pv)) 
        MZttols[Tc, ] <- t(quantile(cbind(Mtols[ ,3]),pv)) 
        MZttgls[Tc, ] <- t(quantile(cbind(Mtgls[ ,3]),pv)) 
        Pc[Tc, ] <- t(quantile(Pcgls,pv)) 
        Pt[Tc, ] <- t(quantile(Ptgls,pv))  
    } 
 
    cat('\nPercentiles of Distributions of Unit Root Test Statistics') 
    cat('\n---------------------------------------------------------') 
 
    cat('\nDF a test') 
    cat('\n---------------------------------------------------------') 
    cat('\nNo de-trending') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), adf1=adf1))
    cat('\n')
    cat('\nOLS,constant') 
    cat('\n      T         1%        5%         10%\n')
    print(cbind(Tv=t(Tv), adf1cols=adf1cols))    
    cat('\nOLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n')
    print(cbind(Tv=t(Tv), adf1tols=adf1tols))
    cat('\nGLS,constant[ays. equiv to no de-trending]') 
    cat('\n      T         1%        5%         10%\n')
    print(cbind(Tv=t(Tv), adf1cgls=adf1cgls))
    cat('\nGLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), adf1tgls=adf1tgls))
 
    cat('\nDF t test') 
    cat('\n---------------------------------------------------------') 
    cat('\nNo de-trending') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), adf2=adf2))
    cat('\nOLS,constant') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), adf2cols=adf2cols))
    cat('\nOLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), adf2tols=adf2tols))
    cat('\nGLS,constant[ays. equiv to no de-trending]') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), adf2cgls=adf2cgls))
    cat('\nGLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n')
    print(cbind(Tv=t(Tv), adf2tgls))
 
 
 
    cat('\nMZ t test') 
    cat('\n---------------------------------------------------------') 
    cat('\nOLS,constant') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), MZtcols=MZtcols))
    cat('\nOLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), MZttols=MZttols))
 
 
    cat('\nMSB test') 
    cat('\n---------------------------------------------------------') 
    cat('\nNo de-trending') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), MSB=MSB))
    cat('\nOLS,constant') 
    cat('\n      T         1%        5%         10%\n')
    print(cbind(Tv=t(Tv), MSBcols=MSBcols))
    cat('\nOLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), MSBtols=MSBtols))
    cat('\nGLS,constant[ays. equiv to no de-trending]') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), MSBcgls=MSBcgls))
    cat('\nGLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n')
    print(cbind(Tv=t(Tv), MSBtgls=MSBtgls))
 
 
    cat('\nPoint optimal') 
    cat('\n---------------------------------------------------------') 
    cat('\nGLS,constant') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), Pc=Pc))
    cat('\nGLS,constant and linear trend') 
    cat('\n      T         1%        5%         10%\n') 
    print(cbind(Tv=t(Tv), Pt=Pt)) 
} 

#============================================================================
#
#   Program to reproduce Diebold and Yilmaz (2009) spillover model
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  trimr
source("EMTSUtil.R")

#
#------------------------- Diebold-Yilmaz Spillover Index--------------------
#
stsm_spillover <- function() {
  # Read the data
  # stock market returns and volatility (4-Dec-1996 to 23-Nov-2007)
  # Order of countries
  #                         1.  US
  #                         2.  UK
  #                         3.  France
  #                         4.  Germany
  #                         5.  Hong Kong
  #                         6.  Japan
  #                         7.  Australia
  #                         8.  Indonesia
  #                         9.  S. Korea
  #                         10. Malaysia
  #                         11. Philippines
  #                         12. Singapore
  #                         13. Taiwan
  #                         14. Thailand
  #                         15. Argentina
  #                         16. Brazil
  #                         17. Chile
  #                         18. Mexico
  #                         19. Turkey
  
  load('diebold_yilmaz.Rdata')
  
  # Choose the data
  y <- returns        
  stdy <- apply(y, 2, sd)
  z <- t( apply((y-colMeans(y)), 1, "/",stdy))
  
  print(cbind("Mean" = colMeans(y), "Median" = apply(y, 2, median), "Max" =apply(y, 2, max), "Min"= apply(y, 2, min), "Std dev" = apply(y, 2, sd), "Skew" = colMeans(z^3), "Kurt" = colMeans(z^4)))
  
  # Estimate a VAR with lag p=2  
  x <- cbind(rep(0, nrow(y)-2), trimr(y,1,1), trimr(y,0,2))
  var.reg <- lm(trimr(y,2,0) ~ x - 1)
  b <- var.reg$coef
  
  mu  <- t(trimr(b,0,38))  # Vector of intercepts                                              
  phi1 <- t(trimr(b,1,19)) # Lag 1 parameter estimates   
  phi2 <- t(trimr(b,20,0)) # Lag 2 parameter estimates     
    
  # Generate VMA (non-orthogonalized) for horizons 1 to 10 
  si1  <- diag(ncol(y))
  si2  <- phi1
  si3  <- phi1 %*% si2 + phi2
  si4  <- phi1 %*% si3 + phi2 %*% si2
  si5  <- phi1 %*% si4 + phi2 %*% si3
  si6  <- phi1 %*% si5 + phi2 %*% si4
  si7  <- phi1 %*% si6 + phi2 %*% si5
  si8  <- phi1 %*% si7 + phi2 %*% si6
  si9  <- phi1 %*% si8 + phi2 %*% si7
  si10  <- phi1 %*% si9 + phi2 %*% si8
  
  # Generate VMA (orthogonalized) for horizons 1 to 10  
  v  <- var.reg$residuals  # VAR residuals      
  vc <- t(v) %*% v/nrow(v)
  d  <- diag(vc)
  s  <- t(chol(vc))
  
  ir1  <- si1 %*% s
  ir2  <- si2 %*% s
  ir3  <- si3 %*% s
  ir4  <- si4 %*% s
  ir5  <- si5 %*% s
  ir6  <- si6 %*% s
  ir7  <- si7 %*% s
  ir8  <- si8 %*% s
  ir9  <- si9 %*% s
  ir10 <- si10 %*% s
  
  # Compute variance decompositions for horizons 1 to 10 
  vd1  <- ir1^2
  vd1  <- 100* t( apply(vd1, 1, "/", rowSums(vd1)) )
  
  vd2  <- ir1^2 + ir2^2
  vd2  <- 100* t( apply(vd2, 1, "/", rowSums(vd2)) )
  
  vd3  <- ir1^2 + ir2^2 + ir3^2
  vd3  <- 100* t( apply(vd3, 1, "/", rowSums(vd3)) )
  
  vd4  <- ir1^2 + ir2^2 + ir3^2 + ir4^2
  vd4  <- 100*t( apply(vd4, 1, "/", rowSums(vd4)) )
  
  vd5  <- ir1^2 + ir2^2 + ir3^2 + ir4^2 + ir5^2
  vd5  <- 100*t( apply(vd5, 1, "/", rowSums(vd5)) )
  
  vd6  <- ir1^2 + ir2^2 + ir3^2 + ir4^2 + ir5^2 + ir6^2
  vd6  <- 100*t( apply(vd6, 1, "/", rowSums(vd6)) )
  
  vd7  <- ir1^2 + ir2^2 + ir3^2 + ir4^2 + ir5^2 + ir6^2 + ir7^2
  vd7  <- 100*t( apply(vd7, 1, "/", rowSums(vd7)) )
  
  vd8  <- ir1^2 + ir2^2 + ir3^2 + ir4^2 + ir5^2 + ir6^2 + ir7^2 + ir8^2
  vd8  <- 100*t( apply(vd8, 1, "/", rowSums(vd8)) )
  
  vd9  <- ir1^2 + ir2^2 + ir3^2 + ir4^2 + ir5^2 + ir6^2 + ir7^2 + ir8^2 + ir9^2
  vd9  <- 100*t( apply(vd9, 1, "/", rowSums(vd9)) )
  
  vd10 <- ir1^2 + ir2^2 + ir3^2 + ir4^2 + ir5^2 + ir6^2 + ir7^2 + ir8^2 + ir9^2 + ir10^2
  vd10 <- 100*t( apply(vd10, 1, "/", rowSums(vd10)) )
  
  rnames <- c('US     '  ,
           'UK     '  ,    
           'FRA    '  ,
           'GER    '  ,
           'HKG    '  ,
           'JPN    '  ,
           'AUS    '  ,
           'IDN    '  ,    
           'KOR    '  ,    
           'MYS    '  ,    
           'PHL    '  ,    
           'SGP    '  ,    
           'TAI    '  ,    
           'THA    '  ,    
           'ARG    '  ,    
           'BRA    '  ,    
           'CHL    '  ,    
           'MEX    '  ,    
           'TUR    ')
  
  cat('\nVariance decomposition at period 10\n')
  print(vd10) 
  
  tmp <- rowSums(vd10)-diag(vd10)
  cat('\nContribution From Others\n')
  print( matrix(tmp, nrow=19, dimnames=list(rnames, "")))
  
  tmp <- (colSums(vd10)-diag(vd10))
  cat('\nContribution To Others\n')
  print( matrix(tmp, nrow=19, dimnames=list(rnames, "")))
}

#============================================================================
#
#   Estimate a GARCH(1,1) model of daily equity returns with seasonality 
#
#============================================================================
rm(list = ls(all=T))
graphics.off()
#
#--------------------------- Helper Functions -------------------------------
# 
# Load required functions - recserar, trimr
source("EMTSUtil.R")

#-------------------------------------------------------------------------
#  GARCH(1,1) - unrestricted
#-------------------------------------------------------------------------
neglog <- function(b,y,dtue,dwed,dthu,dfri) {
  u <- y  
  d <- b[4]^2*dtue + b[5]^2*dwed + b[6]^2*dthu + b[7]^2*dfri                 
  h <- recserar( cbind(b[1]^2 + b[2]^2*trimr(rbind(0.0,u^2),0,1) + d),cbind(sd(u)^2),cbind(b[3]^2))   
  z <- u/sqrt(h)                                                                       
  f <-  -0.5*log(2*pi) - 0.5*log(h) - 0.5*z^2  
  lf <- -mean( f )
  return(lf)
}


end
#-------------------------------------------------------------------------
#  GARCH(1,1) - restricted
#-------------------------------------------------------------------------
neglog0 <- function(b,y) {
  u <- y                                                                                         
  h <- recserar(cbind(b[1]^2 + b[2]^2*trimr(rbind(0.0,u^2),0,1)),cbind(sd(u)^2),cbind(b[3]^2))                    
  z <- u/sqrt(h)                                                                          
  f <-  - 0.5*log(2*pi) - 0.5*log(h) - 0.5*z^2   
  lf <- -mean( f )
  return(lf)
}


#
#--------------------------- Day of the week effects ------------------------
#

garch_seasonality <- function(  ) {
  # Load daily equity indices 
  # 5 January 1989 to 31 December 2007
  # ftse, dow, nk, dummy vars for day of week (tue, wed, thurs, fri, holiday_ftse, holiday_dow, holiday_nk)
  data <- as.matrix(read.table("equity.dat", 
                               col.names=c('FSTSE', 'DOW', 'NIKKEI', 
                                           'tdum',  'wdum',  'thdum',  'fdum',  'hdum_ftse',  'hdum_dow',	'hdum_nk')))
  equity <- data[,1] # ftse
  
  #     Compute the percentage return 
  y <- 100*(trimr(log(equity),1,0) - trimr(log(equity),0,1))    
  y <- y - mean(y)      
  t <- length(y)
  
  # Trim dummy variables
  dtue <- trimr(data[,4],1,0)                         
  dwed <- trimr(data[,5],1,0) 
  dthu <- trimr(data[,6],1,0) 
  dfri <- trimr(data[,7],1,0)                                      
  
  # Estimate the GARCH(1,1) model with seasonality     
  start   <- c(0.05, 0.1, 0.9, 1, 1, 1, 1)
  estResults <- optim(start, neglog, y=y, dtue=dtue,dwed=dwed,dthu=dthu,dfri=dfri, method="BFGS")
  lf1 <- estResults$val
  
  lf1 <- -lf1
  
  cat('\nLikelihood function (unrestricted)   = ',lf1)
  
  # Estimate the GARCH(1,1) model without seasonality 
  start <- c(0.05, 0.1, 0.9)
  estResults <- optim(start, neglog0, y=y, method="BFGS")
  lf0 <- estResults$val
  
  lf0 <- -lf0
  
  cat('\nLikelihood function (restricted)    = ',lf0)
  
  # LR test
  lr  <- -2*t*(lf0 - lf1)
  pv <- 1-pchisq(lr,5)
  cat('\n')
  cat('\nLR statistic           = ',lr)
  cat('\np-value                = ',pv)
}


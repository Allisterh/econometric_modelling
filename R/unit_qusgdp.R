#============================================================================
#
#   Program to illustrate basic concepts using U.S. GDP
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()


#
#------------------------- Helper Functions----------------------------------
#

#load required functions -  figure
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  Detrending function: 
#       cbar = -7 constant 
#       cbar = -13.5 linear trend  
#       cbar = -T for OLS detrending
#----------------------------------------------------------------------------

glsdetrend <- function( y,x,cbar ){
  t <- length(y)
  yc <- cbind(c(y[1], (trimr(y,1,0)-(1+cbar/t)*trimr(y,0,1))) )
  xc <- rbind(x[1,], cbind(trimr(x,1,0)-(1+cbar/t)*trimr(x,0,1)))  
  reg  <- lm(yc ~ xc - 1)
  return(reg)
}

#
#-------------------------------- U.S. GDP ----------------------------------
#
unit_qusgdp <- function() {
  # Load the data
  gdp <- as.matrix(read.table("usgdp.dat"))
  t  <- length(gdp)
  y <- log(gdp)
    
  # Detrending
  x <- cbind(rep(1, t), seq(from=1,by=1,length.out=t))
  lm <- glsdetrend(y,x,-t)  
  bols <- lm$coef
  uols <- y - x %*% bols

  lm <- glsdetrend(y,x,0)  
  bdif <- lm$coef
  udif <- y - x %*% bdif
 
  lm <- glsdetrend(y,x,-13.5)  
  bgls <- lm$coef
  ugls <- y - x %*% bgls

  # OLS with trend break 1973Q1

  dvec <- seq(from=1947, by=1/4, length.out=t)
  dt   <- cumsum( dvec > 1973 )

  xd   <- cbind(x,  dt) 
  bd   <- lm(y ~ xd - 1)$coef
  yhat <- xd %*% bd
  ud   <- y - yhat

  # Linear trend  no break at 1973Q1
  yhat1 <- xd %*% cbind(c(bd[1:2], 0))

  #**********************************************************************
  #***
  #***     Generate graphs
  #***
  #**********************************************************************

  figure()
  par(xaxs="i", yaxs="i")

  #--------------------------------------------------------#
  plot(dvec,y, type="l",
       xlab = "Year",
       ylab = "Log of Real GDP",
       xlim = c(1940, 2010),
       ylim = c(7.0, 10.0),
       bty = "l")

  figure()
  par(xaxs="i", yaxs="i")
  
  #--------------------------------------------------------#
  matplot(dvec,cbind(uols, udif, ugls), type="l",
          main = "Plots or Residuals",
          lty = c(1,2,3), col=3:5,
          xlab = "Year",
          ylab = expression(hat(u[t])),
          xlim = c(1940, 2010),
          ylim = c(-0.15, 0.15),
          bty = "l")
  legend("topright", 
           legend=c("OLS", "1st DIF", "GLS"),
           lty=c(1,2,3), col=3:5)
  

  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
  #--------------------------------------------------------#
  matplot(dvec,cbind(y,yhat,yhat1), type="l",
          main = "(a)",
          xlab = "Year",
          ylab = "Time Trends, Real GDP",
          xlim = c(1940, 2010),
          ylim = c(7.0, 10.0),
          bty = "l")
  #--------------------------------------------------------#  
  matplot(dvec,cbind(uols,ud), type="l",
          main = "(b)",
          xlab = "Year",
          ylab = expression(hat(u[t])),
          xlim = c(1940, 2010),
          ylim = c(-0.15, 0.15),
          bty = "l")
}

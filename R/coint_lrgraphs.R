#============================================================================
#
#   Program to plot US macroeconomic data and 
#   to identify long-run properties of economic models
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - seqa, figure
source("EMTSUtil.R")

# Load required library - scatterplot3d
library("scatterplot3d")

# 
#-------------- Graphical Analysis of Long-Run Economic Theories ------------
#

coint_lrgraphs <- function() 
{  
  # Load perminent income data (total period is 1947q1 to 2010q2)
  data <- as.matrix(read.table("permincome.dat"))
  rcpc <- data[,1]
  rypc <- data[,2]
  
  # Select desired sample (1984Q1 to 2005Q4)
  c    <- rcpc[150:237]
  y    <- rypc[150:237]
  
  cat('\nOLS estimates of the permanent income equation (1984Q1 to 2005Q4)\n')
  x <- cbind(rep(1, length(y)), log(y))
  bhat <- lm(log(c) ~ x - 1)$coef
  cat(bhat)
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
      
  dvec <- seqa( 1984,1/4,length(c) )
  #--------------------------------------------------------#
  # Panel (a)
  matplot( dvec,cbind(log(c),log(y)), type="l",
           main = '(a)',
           ylab = expression(lrc[t], lry[t]),
           xlab = "Year", 
           ylim = c(9.8, 10.4),
           bty = "l")
  
  
  #--------------------------------------------------------#
  # Panel (b)
  plot(log(c),log(y), pch=19,
           main = '(b)',
           ylab = expression(lrc[t]),
           xlab = expression(lry[t]),
           xlim = c(9.8, 10.4),
           ylim = c(9.8, 10.4),
           bty = "l")
  
  # Load money demand data (1959q1 to 2005q4)  
  data <- as.matrix(read.table("moneydemand.dat"))
  
  cpi <- data[,1]
  fedfunds <- data[,2]
  gdp <- data[,3]
  m2 <- data[,4]
  tbill <- data[,5]
  
  
  lrm    <- log(m2/cpi)
  lry    <- log(gdp/cpi)
  spread <- tbill/100 - fedfunds/100
  
  x    <- cbind(rep(1,length(lry)), lry, spread)
  bhat <- lm(lrm ~ x - 1)$coef
  
  cat('\n\nOLS estimates of the money demand (1959Q1 to 2005Q4)\n')
  cat(bhat) 
  
  u <- lrm - x %*% bhat
  
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
      
  #--------------------------------------------------------#
  # Panel (a)
  scatterplot3d(lry,spread,lrm,pch=19,
                main = "(a)",
                xlab = expression(lry[t]),
                ylab = "Spread",
                zlab = expression(lrm[t]))
  
  #--------------------------------------------------------#
  # Panel (b)
  plot(seqa(1959,1/4,188),u, type="l",
       main = '(b)',
       ylab = "Residuals",
       xlab = "Year",
       xlim = c(1959, 2005),
       ylim = c(-0.2, 0.2),
       bty = "l")  
  
  
  # Load interest rate data (March 1962 to September 2010)
  data <- as.matrix(read.table("usmacro.dat"))
  r10yr <- data[,1]  
  r1yr <- data[,2]
  r5yr <- data[,3]
  
  #**********************************************************************
  #***
  #***     Generate graph
  #***
  #**********************************************************************
  
  figure()
  par(xaxs="i", yaxs="i", mfrow=c(1,2))
      
  dvec <- seqa( 1962,1/4,length(r1yr) )
  #--------------------------------------------------------#
  # Panel (a)
  matplot(dvec,cbind(r1yr,r5yr,r10yr), type="l",
          main = '(a)',
          ylab = "Yield %",
          xlab = "Year",
          bty = "l")
  
  #--------------------------------------------------------#
  # Panel (b)
  scatterplot3d(r10yr,r5yr,r1yr, pch=19,
                main = '(b)',
                xlab = '10 year',
                ylab = '5 year',
                zlab = '1 year',
                bty = "l")
}

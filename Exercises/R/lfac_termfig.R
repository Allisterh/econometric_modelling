#============================================================================
#
#   Program to plot term structure data 
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#

# Load required functions - figure
source("EMTSUtil.R")

# Read the data:
# US daily zero coupon yields starting 10/04/1988 and ending 12/28/2001
#   The variables are:
#
#		1.	tcm3m
#		2.	tcm1y
#		3.	tcm3y
#		4.	tcm5y
#		5.	tcm7y
#		6.	tcm10y

#
#------------------- Term Structure of Interest Rates plots ------------------
#
data <- as.matrix(read.table("lfac_usdata.dat"))

# Set up dates
t <- seq(1988+10/12, 2002.057333333333, 1/250)
t <- t[1:nrow(data)]


#**********************************************************************
#***
#***     Generate graphs
#***
#**********************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(3,2))

#--------------------------------------------------------#
# Panel (a)
plot(t,data[,1],type="l",
     main = '(a) 3 Month Yield',
     xlab = 't',
     ylab = '%',     
     bty = "l")

#--------------------------------------------------------#
# Panel (b)
plot(t,data[,2],type="l",
     main = '(b) 1 Year Yield',
     ylab = '%',
     xlab = 't',
     bty ="l")

#--------------------------------------------------------#
# Panel (c)
plot(t,data[,3],type="l",
     main = '(c) 3 Year Yield',
     ylab = '%',
     xlab = 't',
     bty = "l")

#--------------------------------------------------------#
# Panel (d)
plot(t,data[,4],type="l",
     main = '(d) 5 Year Yield',
     ylab = '%',
     xlab = 't',
     bty = "l")

#--------------------------------------------------------#
# Panel (e)
plot(t,data[,5],type="l",
     main = '(e) 7 Year Yield',
     ylab = '%',
     xlab = 't',
     bty = "l")
#--------------------------------------------------------#
# Panel (f)
plot(t,data[,6],type="l",
     main = '(f) 10 Year Yield',
     ylab = '%', 
     xlab = 't',     
     bty = "l")
      

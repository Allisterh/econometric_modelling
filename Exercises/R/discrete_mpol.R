#=========================================================================
#
#   Plot US Federal Funds target rate
#
#=========================================================================
rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  figure
source("EMTSUtil.R")

# Read the data: loads a time series object "target"
# US weekly yields starting 1st week of Feb 1984 -- first week of June 1997

data <- read.table("FedFunds.dat")
target <- data[,2]
dates <- as.Date(data[,1], format="%d-%b-%Y")

#**********************************************************************
#***
#***     Generate graph
#***
#**********************************************************************
  
figure()

plot(dates,target,type="l",
     main="",
     ylab = 'Federal Funds Rate',
     xlab = 'Time',
     bty = "l")
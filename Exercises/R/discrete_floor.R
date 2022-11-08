#============================================================================
#
#   Money demand equation with a floor interest rate
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(42, kind="Mersenne-Twister")

#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  figure
source("EMTSUtil.R")

t  <- 20
x  <- seqa(0,1,t)
b0 <- 10
b1 <- -0.5

m <- b0 + b1*x
s <- 1
c <- 4

y  <- m + s*rnorm(t,1)
f1 <- pnorm( (c - m)/s )
f2 <- dnorm( (c - m)/s )
pred <- m + s*f2/(1 - f1)


#****************************************************************************
#***
#***     Generate graph
#***
#****************************************************************************

figure()

matplot(x,cbind(m,c*rep(1,t), pred), type="l",
        ylab = 'Interest rate (#)',
        xlab = 'Money',
        bty = "l")

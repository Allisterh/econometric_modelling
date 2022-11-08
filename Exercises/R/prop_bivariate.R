#=========================================================================
#
#   Bivariate normal distribution properties
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

y1 <- seq(-4, 4, 0.2)
y2 <- seq(-4, 4, 0.2)

#********************************************************************
#***
#***     Generate graph
#***
#********************************************************************

mu  <- c(0, 0)
Sig <- rbind(c(1, 0.6), c(0, 0.6))

f <- array(0, c(length(y1),length(y2)) )
for (i in seq(y1)) {
  
    for (j in seq(y2)) {
        y1y2 <- c(y1[i], y2[j])        
        
        lhs <- ( 1/(2*pi*det(Sig)^1/2) )
        rhs <- exp(-0.5 * (y1y2 - mu) %*% solve(Sig) %*% (cbind(y1y2) - mu) )
        
        f[i,j] <- lhs * rhs
    }
}

figure()
par(mfcol=c(2,2))

# Create colored persp plots
fill <- matrix("", nrow = nrow(f)-1, ncol = ncol(f)-1)
fcol <- fill
fcol[] <- terrain.colors(nrow(fcol))
sh <- 0.15
persp(y2, y1, f, theta = 30, phi = 30, 
      ticktype="detailed", nticks = 5, 
      col = fcol, shade = sh,
      main = "rho = 0.6",
      ylab = "\ny1",
      xlab = "\ny2",
      zlab = "\n\n\nf(y1, y2)")

contour(y1, y2, f,
        xlab = "y1",
        ylab = "y2")

mu  <- c(0, 0)
Sig <- rbind(c(1, 0), c(0, 1))

f <- array(0, c(length(y1),length(y2)) )
for (i in seq(y1)) {
  
    for (j in seq(y2)) {
        y1y2 <- c(y1[i], y2[j])        
        
        lhs <- ( 1/(2*pi*det(Sig)^1/2) )
        rhs <- exp(-0.5 * (y1y2 - mu) %*% solve(Sig) %*% (cbind(y1y2) - mu) )
        
        f[i,j] <- lhs * rhs
    }
}

persp(y2, y1, f, theta = 30, phi = 30, 
      ticktype="detailed", nticks = 5,
      col= fcol, shade = sh, 
      main = "rho = 0",
      ylab = "\ny1",
      xlab = "\ny2",
      zlab = "\n\n\nf(y1, y2)",
      zlim = c(0,0.4))

contour(y1, y2, f,
        xlab = "y1",
        ylab = "y2")

#============================================================================
#
#   Program to generate the sampling distribution of the eigenvalues of
#   a bivariate vecm with rank r=1
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(1234)

# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - figure, trimr
source("EMTSUtil.R")
# Load required library - flipud, zeros, kde
library(matlab)
library(ks)

Tv   <- c(100,200,400,800)
reps <- 10000
lam  <- zeros(reps,2)
lr   <- zeros(reps,length(Tv))


for (tc in seq(length(Tv))) {
  t <- Tv[tc]
  
  for (rep in seq(reps)) {
    # Simulate dgp
    y2 <- cumsum(rnorm(t))
    y1 <- y2 + rnorm(t)
    y  <- cbind(y1, y2)
    
    r0 <- trimr(y,1,0)-trimr(y,0,1)     
    r1 <- trimr(y,0,1)
    
    # Construct sums of squares matrices     
    tmp <- nrow(r0)
    s00 <- t(r0) %*% r0/tmp                  
    s11 <- t(r1) %*% r1/tmp
    s01 <- t(r0) %*% r1/tmp
    s10 <- t(s01)

    # Choleski decomposition of s11         
    l <- t(chol(s11))
            
    # Compute eigenvalues  
    
    lam[rep,] <- eigen( inv(l) %*% s10 %*% inv(s00) %*% s01 %*% inv(t(l)) )$values
    # Compute trace statistic (smallest eigenvalue is zero)  
    lr[rep,tc] <- -t*log(1-lam[rep,2])       
  }
  cat('\n     T       Mean 1    Mean 2       StDev1      StDev2\n')
  print(rbind(c(t, colMeans(lam), apply(lam, 2, sd))))
  cat( '\n')                       
}

#  Compute and plot the kernel density estimator of the lr statistic 
#  and compare with chisq dof=1   

minx <- 0.5 
maxx <- 12
xi <- seqa(minx,(maxx-minx)/200,201)

fchi <- dchisq(xi,1)

# Pick last column of lr corresponding to t=800 to compute kernel estimate    
fhat <- kde(lr[,tc],h = 0.07, eval.points=xi)

#**********************************************************************
#***
#***     Generate graph
#***
#**********************************************************************
figure()

plot(fhat,
        ylab="f(LR)",
        xlab = "LR",
        bty = "l")
lines(xi, fchi, lty=3)

#=========================================================================
#
#   Program to estimate an Poisson model and plot the 
#   log-likelihood function.
#
#=========================================================================

# clear all
rm(list = ls(all = TRUE))
graphics.off()

# load any compatibility requirements
source("EMTSUtil.R")

# Data
y <- c(8,3,4)    # Data used in Poisson example    
#y <- c(6, 2, 3, 1)     

t <- length(y)
theta <- mean(y)
lnl_t <- y*log(theta) - theta - log(factorial(y))
g_t <- y/theta - 1
h_t <- -y/theta^2

cat('\n')
cat('Sum of y =', sum(y),'\n')
cat('Mean of y =', theta, '\n')
cat('Log-likelihood function =', mean(lnl_t),'\n')

rnames <- seq(t)
cnames <- c("yt", "lnlt", "gt", "ht")
print(matrix(c(y, lnl_t, g_t, h_t), 
             dimnames=list(rnames, cnames), nrow=t))


# ***    Generate graph   ***

theta <- seq(0.01, 15, 0.01)
lnl <- array(0, c(length(theta), 1) )

for(i in seq(theta))
{
  lnl[i] <- mean( y * log(theta[i]) - theta[i] - log( factorial(y) ))
}

# create subplots

figure()
layout(matrix(c(1,2), byrow=TRUE, ncol=2))

# generate line graph
plot(theta, lnl, type="l",
     col="darkblue",
     main = "(a) Log-likelihood function",
        ylab = expression(paste("ln L"[T]*"","(",theta, ")")),
        xlab = expression(theta),
        bty="l")

# generate custom bar graph
plot(y, abs(lnl_t), type="h", lwd=20,                  
                  main="(b) Log-density function",
                  col="darkblue", yaxt ="n",     
                  ylab = expression(paste("ln f(y"[t]*";",hat(theta)," = 3)")),
                  xlab = expression("y"[t]*""),
                  bty="l")
yaxis.vector <- seq(0, max(abs(lnl_t))+1, 0.2)
axis(2, at = yaxis.vector, labels = -yaxis.vector, xaxs="i", yaxs="i")


#============================================================================
#
#     Program to demonstrate the distribution of the t-statistic
#     from a Student t distribution with degrees of freedom of nu={1,2,3}.
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()

set.seed(12, kind="Mersenne-Twister")


# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - figure
source("EMTSUtil.R")

# 
#--------------------- Central Limit Theorem Student t-----------------------
#
t <- 500                              
r <- 5000                            

#  Generate Student t (nu<-1) ie Cauchy distribution     
nu     <- 1
chi1   <- matrix(rnorm(t*r)^2, t, r)
rstud1 <- matrix(rnorm(t*r)/sqrt( chi1/nu ), t, r)

z1 <- sqrt(t)*colMeans(rstud1)/apply(rstud1, 2, sd)


figure()
x <- seq(min(z1), max(z1), by=0.01)
hist(z1, breaks=41, prob=T, col='blue',
     main='Student t with nu=1',
     xlab = 'Midpoint')
curve(dnorm(x, mean=mean(z1),sd= sd(z1)), add=T)

#  Generate Student t (nu=2) 
nu     <- 2
chi2   <- matrix(rnorm(t*r)^2 + rnorm(t*r)^2, t, r)
rstud2 <- matrix(rnorm(t*r)/sqrt( chi2/nu ), t, r)

z2 <- sqrt(t)*colMeans(rstud2)/apply(rstud2, 2, sd)

figure()

hist(z2, breaks=41, prob=T, col='blue',
     main= 'Student t with nu=2', 
     xlab = 'Midpoint')
x <- seq(min(z2), max(z2), by=0.01)
curve(dnorm(x, mean=mean(z2),sd= sd(z2)), add=T)

#  Generate Student t (nu=3)  
nu     <- 3
chi3   <- matrix(rnorm(t*r)^2 + rnorm(t*r)^2, t, r)
rstud3 <- matrix(rnorm(t*r)/sqrt( chi2/nu ), t, r)

z3 <- sqrt(t)*colMeans(rstud3)/apply(rstud3, 2, sd)

figure()
hist(z3, breaks=41, prob=T, col="blue",
     main = 'Student t with nu=3',
     xlab = 'Midpoint')
x <- seq(min(z3), max(z3), by=0.01)
curve(dnorm(x, mean=mean(z3),sd= sd(z3)), add=T)

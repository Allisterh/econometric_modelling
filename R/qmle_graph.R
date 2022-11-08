#==============================================================================
#
#   Comparing the true likelihood (lnl0) and the quasi-likelihood (lnl) 
#
#==============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------------
#

#load required functions - figure
source("EMTSUtil.R")

#
#--------------------  QMLE Graphical Analysis --------------------------------
#

t <- 10

# Exponential versus normal
mu <- 1
y  <- -mu*log(1 - runif(t))                 
y <- round(y*1000)/1000            # Round y to three decimal places

# True likelihood (exponential)
theta0 <- seq(from=0.01,by=0.1,length.out=521)
lnl0   <- -log(theta0) - mean(y)/theta0         

# Quasi-likelihood (normal with unit variance)
theta <- seq(from=-5,by=0.1,length.out=201)
sig2  <- 1
lnl   <- rep(0, length(theta))

for (i in seq(theta)) {
   lnl[i]   <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*mean( (y - theta[i] )^2 )/sig2  
}    

#********************************************************************
#***
#***     Generate graph
#***
#********************************************************************

figure()
par(xaxs="i", yaxs="i", mfrow=c(1,2))

plot(theta0,lnl0, type="l",
     main = '(a) Exponential - Normal',
     xlab = expression(theta),
     ylab = expression(paste( "ln L"[T]*"", (theta) )),
     xlim = c(-1,3),
     ylim = c(-2, -1))
lines(theta,lnl,lty=2)

   
# Negative Binomial versus Poisson
mu <- 5
p  <- 0.5
y <- rnbinom(t, prob=p, size=mu)                                         


# True likelihood (negative binomial)
theta0 <- seq(from=0.01,by=0.1,length.out=201)
lnl0   <- rep(0, length(theta0))

for (i in seq(theta0)) {
  lnl0[i] <- mean(log(gamma(y + theta0[i]))) - mean(log(gamma(y + 1)))  - log(gamma(theta0[i])) + theta0[i]*log(1-p) + mean(y)*log(p)  
}

# Quasi-likelihood (Poisson)
theta <- theta0
lnl <- mean(y)*log(theta) - theta - mean(log(factorial(y)))  

plot(theta0,lnl0, type="l",
     main = "(b) Negative Binomial - Poisson",
     xlab = expression(theta),
     ylab = expression(paste( "ln L"[T]*"", (theta) )),
     xlim = c(0,10),
     ylim = c(-6, -2))
lines(theta,lnl,lty=2)


# Student t versus Normal
mu  <- 5
sig <- 1
v   <- 10

# Standardized Student t random numbers
y <- mu + sig*sqrt((v-2)/v)*pt(runif(t),v)                    


# True likelihood (Student t)
theta0 <- seq(from=-5, by=0.1,length.out=201)
z      <- array(0, c(t,length(theta0)) )
for (i in seq(length(theta0))) {
  z[,i] <- y-theta0[i]
}
const  <- gamma( (v+1)/2 ) / ( sqrt(pi*(v-2)) * gamma( v/2 ) )               

lnl0 <- log(const) - 0.5*log(sig^2) - 0.5*(v+1)*colMeans( log( 1 + (z^2)/(v-2) ) )            

# Quasi-likelihood (normal)
theta <- theta0
lnl <- -0.5*log(2*pi) - 0.5*log(sig^2) - 0.5*colMeans( z^2 )         

figure()
par(xaxs="i", yaxs="i")
plot(theta,lnl, type="l", lty=2,
     main = "t - Normal",
     xlab = expression(theta),
     ylab = expression(paste( "ln L"[T]*"", (theta) )) )
lines(theta0,lnl0)

 
# Poisson versus normal
mu <- 5

# Poisson random numbers 
y <- rpois(t, mu)

# True likelihood (Poisson)
theta0 <- seq(from=0.01,by=0.1,length.out=201)
lnl0 <- mean(y)*log(theta0) - theta0 - mean(log(factorial(y)))  

# Quasi-likelihood (normal with variance <- mu)
theta <- seq(from=-5,by=0.1,length.out=201)
sig2  <- mu 
for (i in seq(theta0)) {
  z[,i] <- y-theta0[i]/sig2
}

lnl <- -0.5*log(2*pi) - 0.5*log(sig2) - 0.5*colMeans( z )

figure()
par(xaxs="i", yaxs="i")
plot(theta0,lnl0, type="l",
     main = "Poisson - Normal",
     xlab = expression(theta),
     ylab = expression(paste( "ln L"[T]*"", (theta) )),
     xlim = c(-5, 25))
    
lines(theta,lnl, lty=2)

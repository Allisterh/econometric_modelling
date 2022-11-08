#=========================================================================
#
#    Program to demonstrate the Lindberg-Feller central limit theorem
#    using a regression model with gamma distributed errors.
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123456, kind="Mersenne-Twister")

t      <- 500
NDraws <- 5000                       # Number of draws                                    


b0    <- cbind(c(1,2))                   # Population parameters   
rho   <- 0.25                        # Parameters of the gamma distribution     
alpha <- 0.1


trigs <- matrix(rgamma(t*NDraws, rho), nrow = t)
x0    <- cbind(rep(1, t), rnorm(t))

# For the sample size of T = 100  
t <- 100
x  <- x0[1:t,]

zeros <- rep(0, NDraws)
z1 <- zeros
z2 <- zeros

pb <- txtProgressBar(min= 0, max = NDraws, style = 3)
for (i in seq(NDraws) ) {

    u    <- alpha * trigs[1:t,i] - rho * alpha  
    
    y    <- x %*% b0 + u
    b    <- as.matrix( lm(y~x[,2])$coef )
    e    <- y - x %*% b
    s2   <- as.numeric(t(e) %*% e/t)    
    vcov <- s2 * solve(t(x) %*% x)
    z1[i]<- (b[1] - b0[1]) / sqrt(vcov[1,1])
    z2[i]<- (b[2] - b0[2]) / sqrt(vcov[2,2])   
    setTxtProgressBar(pb, i)
}
close(pb)


# For the sample size of T = 500
t <- 500
x  <- x0[1:t,]

zeros <- rep(0, NDraws)
z3 <- zeros
z4 <- zeros

pb <- txtProgressBar(min= 0, max = NDraws, style = 3)
for (i in seq(NDraws) ) {

    u    <- alpha * trigs[1:t,i] - rho * alpha  
    
    y    <- x %*% b0 + u
    b    <- as.matrix( lm(y~x[,2])$coef )
    e    <- y - x %*% b
    s2   <- as.numeric(t(e) %*% e/t)    
    vcov <- s2 * solve(t(x) %*% x)
    z3[i]<- (b[1] - b0[1]) / sqrt(vcov[1,1])
    z4[i]<- (b[2] - b0[2]) / sqrt(vcov[2,2])   
    setTxtProgressBar(pb, i)
}
close(pb)

#**************************************************************************
#**     
#**     Generate graph       
#**                                         
#**************************************************************************

figure()

par(mfrow=c(2,2), xaxs="i", yaxs="i")
hist(z1, breaks = 21,     
     main = expression(
       paste("(a) Distribution of z"[hat(beta)[0]]*" (T = 100)")),
     ylab = expression(f(z[hat(beta)[0]])),
     xlab = expression(z[hat(beta)[0]]) )
     
hist(z2, breaks = 21,
     main = expression(
       paste("(b) Distribution of z"[hat(beta)[1]]*" (T = 100)")),
     ylab = expression(f(z[hat(beta)[1]])),
     xlab = expression(z[hat(beta)[1]]) )


hist(z3, breaks = 21,
     main = expression(
       paste("(c) Distribution of z"[hat(beta)[0]]*" (T = 500)")),
     ylab = expression(f(z[hat(beta)[0]])),
      xlab = expression(z[hat(beta)[0]]) )

hist(z4, breaks = 21,
     main = expression(
       paste("(d) Distribution of z"[hat(beta)[1]]*" (T = 500)")),
     ylab = expression(f(z[hat(beta)[1]])),
     xlab = expression(z[hat(beta)[1]]) )

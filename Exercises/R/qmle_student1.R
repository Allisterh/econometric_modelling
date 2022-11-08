#===============================================================================
#   Program to demonstrate the effects of missspecifying the 
#   log-likelihood function on the relationship between the 
#   information matrix and the outer product of gradients.
#
#   True model is student t, misspecified model is normal.
#===============================================================================

 
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")

#
#------------------------- Helper Functions------------------------------------
#

#
#--------------------  QMLE Student t Distribution ----------------------------
#

mu  <- 1     # Population mean      
sig <- 1     # Population standard deviation  
gam <- 5    # Degrees of freedom
t   <- 500000

# Generate data from the true model (Student t)
v <- qt( runif(t), gam )    
y <- mu + sig*sqrt( (gam-2)/gam )*v

# Estimate the misspecified model (normal)    
m  <- mean(y)
s2 <- mean( (y - m)^2 )

# Compute gradients of the misspecified model
g1 <- (y - m)/s2
g2 <- -0.5/s2 + 0.5*(y - m)^2/s2^2
g <- cbind(g1,  g2)
 
# Compute Hessian and information matrix of the misspecified model  
H <- array(0, c(2,2))

H[1,1] <- -1/s2
H[1,2] <- -mean(y - m)/s2^2
H[2,1] <- H[1,2]
H[2,2] <- 0.5/s2^2 - mean((y - m)^2)/s2^3
    
I  <- -H
    
# Compute opg of the misspecified model (normal) 
J <- t(g) %*% g/t
    
cat('\nInformation Matrix')
cat('\n------------------\n')
print(unname(I))
    
  
cat('\nOPG Matrix' )
cat('\n------------------\n')
print(unname(J))
    

cat('\nDifference' )
cat('\n------------------\n')
print(unname(I-J))

 

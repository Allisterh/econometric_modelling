# 

#=======================================================================
#
#    Program to estimate a Normal model with unknown mean 
#   and known variance equal to one
#
#=======================================================================

# clear all
rm(list = ls(all = TRUE))
graphics.off()
 
# Sample data for y  
y <- cbind(1,2,5,1,2)


t <- length( y );     # Define the sample size

# Compute the MLE
theta <- mean(y)

cat('\nMLE\n')
cat('---\n')
cat(theta, '\n')


# Define the log of the likelihood at each observation

lnlt <- -0.5*log(2*pi) - 0.5*(y - theta)^2

# Evaluate log like at each obs
cat('\n')
cat('\nLog like at MLE for each obs\n')
cat('----------------------------\n')
rnames <- seq(t)
cnames <- c("yt", "lnlt")
print(matrix(c(y, lnlt), 
             dimnames=list(rnames, cnames), nrow=t))


# Evaluate log likelihood 
 
cat('\nLog Likelihood at MLE\n')
cat('---------------------\n')
cat(mean(lnlt), '\n')

 
# Evaluate the gradient at the MLE  for each obs
 
g_t <- y - theta
 
cat('\nGradient of the log like at MLE for each obs\n')
cat('--------------------------------------------\n')
cat(g_t, '\n')
 
# Evaluate the gradient at the MLE
g <- mean(g_t)
 
cat('\nGradient of the log likelihood at MLE\n')
cat('-------------------------------------\n')
tol <- 1e-9
cat(ifelse(g < tol,0.0,g), '\n')


# Evaluate the Hessian at the MLE  for each obs
h <- -t

cat('\nHessian of the log likelihood at MLE\n')
cat('------------------------------------\n')
cat(h, '\n')

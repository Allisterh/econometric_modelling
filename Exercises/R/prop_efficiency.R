#========================================================================================
#  Program to demonstrate the efficiency property of MLE for the normal distribution
#
#  Note that the results will not match the numbers reported in the
#  text exactly because of differences in the Gauss and Matlab
#  random number generation.
#=========================================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(123457, kind="Mersenne-Twister")

mue  <- 1       # Population mean         
sig2 <- 2       # Population variance     
r    <- 10000   # Number of replications                            
t    <- 100     # Sample size             


# Generate sample means    

u <- matrix(rnorm(t*r), nrow=t, ncol=r)     # Generate N(0,1) random numbers     
y <- mue + sqrt(sig2) * u                   # Generate N(mue,sig2) random numbers    
meany   <- apply(y, 2, mean)               # Compute the means of the r samples - equivalent to colMeans(y)    
mediany <- apply(y, 2, median)

var_meany   <- mean( (meany   - mue)^2 )
var_mediany <- mean( (mediany - mue)^2 )

cat('\nTheoretical variance of the sample mean   = ', sig2/t, '\n')
cat('Simulated variance of the sample mean     = ', var_meany, '\n')

cat('\nTheoretical variance of the sample median = ', pi*sig2/(2*t), '\n')
cat('Simulated variance of the sample median   = ', var_mediany, '\n')


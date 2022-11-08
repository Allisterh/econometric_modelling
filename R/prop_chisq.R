#==========================================================================
#
#     Monte Carlo program to demonstrate the Central Limit Theorem
#     where 10000 numbers are drawn from a Chi-squared distribution
#     with one degree of freedom.
#
#     For each sample of size 5, the sample mean is computed and the
#     standardized random variable constructed where the population
#     mean and variance are equal to 1 and 2 respectively for the
#     Chi-squared distribution with one degree of freedom
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()

# Load any compatibility requirements
source("EMTSUtil.R")

set.seed(123457, kind="Mersenne-Twister")
        
t <- 5                                    # Sample size                   
r <- 1000                           	    # Number of draws                 
rndNorm <- matrix(rnorm(t*r), nrow = t)  	# Generate N(0,1) random numbers   
rchi1 <- rndNorm^2                   	    # Chi-squared (1) random numbers   
z <- sqrt(t)*(colMeans(rchi1) - 1)/sqrt(2)   # Standardized random variable

figure()
hist(z,breaks = 21, col="darkblue")         # Plot the histogram of z

cat('\nReplications       = ', r)
cat('\nSample size        = ', t)
cat('\nMean               = ', mean(z))
cat('\nVariance           = ', sd(z)^2, '\n')
cat('\nStandard deviation = ', sd(z), '\n')

zsort <- sort(z)                 # Sort the data

cat('\nThe empirical distribution\n')
cat('\n**************************\n')
cat('\n 1.0 per cent      = ', zsort[0.01*r])
cat('\n 2.5 per cent      = ', zsort[0.025*r])
cat('\n 5.0 per cent      = ', zsort[0.05*r])
cat('\n10.0 per cent      = ', zsort[0.10*r])
cat('\n50.0 per cent      = ', zsort[0.50*r])
cat('\n90.0 per cent      = ', zsort[0.90*r])
cat('\n95.0 per cent      = ', zsort[0.95*r])
cat('\n97.5 per cent      = ', zsort[0.975*r])
cat('\n99.0 per cent      = ', zsort[0.99*r], '\n')


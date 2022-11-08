#==========================================================================
#     Program to demonstrate the biasedness property of MLE for the 
#     normal distribution and an adjustment to achieve unbiasedness.
#
#     Note that the results will not match the numbers reported in the
#     text exactly because of differences in the Gauss and Matlab
#     random number generation.
#==========================================================================

rm (list = ls(all=TRUE))
graphics.off()

set.seed(123457, kind="Mersenne-Twister")


mue  <- 1       # Population mean        
sig2 <- 2       # Population variance    
r    <- 20000   # Number of replications                            
t    <- 5       # Sample size            


#     Generate sample means assuming population mean is unknown      

u   <- matrix(rnorm(t*r), nrow = t)                            # Generate N(0,1) random numbers          
y   <- mue + sqrt(sig2) * u                                    # Generate N(mue,sig2) random numbers     

colMeany <- matrix(colMeans(y), nrow = 1)
tmp <- kronecker(matrix(1,nrow(y),1),colMeany)                 # Repeat copies of colMeans(y) by nrow(y) in a matrix


vary <- colSums( (y - tmp)^2 )/t                               # Compute the MLEs of the r samples       
vary_unbiased <- colSums( (y - tmp)^2 )/(t-1)                  # Compute the unbiased estimates of the r samples       

cat('\nResults based on unknown population mean\n')
cat('****************************************\n')
cat('Theoretical value                               = ', sig2, '\n')
cat('Simulated expectation of the MLE                = ', mean(vary), '\n')
cat('Simulated expectation of the unbiased estimator = ', mean(vary_unbiased), '\n\n')


#     Generate sample means assuming population mean is known
u   <- matrix(rnorm(t*r), nrow = t)                            # Generate N(0,1) random numbers          
y   <- mue + sqrt(sig2) * u                                    # Generate N(mue,sig2) random numbers     

vary          <- colSums( (y - mue)^2 )/t     	               # Compute the MLEs of the r samples      
vary_unbiased <- colSums( (y - mue)^2 )/(t-1) 	               # Compute the unbiased estimates of the r samples      


cat('\nResults based on known population mean\n')
cat('****************************************\n')
cat('Theoretical value                               = ', sig2, '\n')
cat('Simulated expectation of the MLE                = ', mean(vary), '\n')
cat('Simulated expectation of the unbiased estimator = ', mean(vary_unbiased), '\n\n')



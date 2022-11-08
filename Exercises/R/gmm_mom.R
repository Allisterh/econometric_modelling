#============================================================================
#
#   Program to compute method of moment estimates of various distributions
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(12, kind="Mersenne-Twister")
# 
#---------------------------------- Helper functions ------------------------
#
# Load require functions - inv, trimr
source("EMTSUtil.R")


t <- 10                            

# Generate data and compute sam
yt <- round(5 + 2*rnorm(t))    
m1 <- mean(yt)  
mt <- cbind(yt , (yt^2) , (yt^3) , (yt^4) , (yt - m1) , (yt - m1)^2 , log(yt) , (1/yt))
m  <- colMeans(mt) 

print(rbind(mt, m))
cat('\n')

# Compute various moments used for method of moments estimation
m2 <- mean(yt^2)
c2 <- mean((yt - m1)^2)
c4 <- m[4]-4*m[3]*m[1]+6*m[2]*m[1]^2-4*m[1]*m[1]^3+m[1]^4      
h1 <- mean(1/yt)

cat('\nNormal distribution (mu)     = ',m1)
cat('\nNormal distribution (sig2)   = ',c2) 
cat('\n')

cat('\nNormal distribution (mu)    = ',m1) 
cat('\nNormal distribution (sig2)  = ',sqrt(c4/3))
cat('\n ' )

cat('\nStudent t distribution (mu) =  ',m1)
cat('\nStudent t distribution (nu) = ',2*c2/(c2-1)) 
cat('\n ' )

cat('\nStudent t distribution (mue) = ',m1) 
cat('\nStudent t distribution (nu) based on solving a quadratic equation in nu')
cat('\n     first root              = ',(6+sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4)))    
cat('\n     second root             = ',(6-sqrt(36 - 4*(1-3/c4)*8))/(2*(1-3/c4)))    
cat('\n ' )

cat('\nGamma distribution (alpha)   = ',m1/(m1-h1))
cat('\nGamma distribution (beta)    = ',(m1/(m1-h1))/m1)
cat('\n ' )

cat('\nGamma distribution (alpha)   = ',m1^2/(m2-m1^2)) 
cat('\nGamma distribution (beta)    = ',m1/(m2-m1^2))
cat('\n ' )

ymin <- min(yt)
cat('\nPareto distribution (alpha)  = ',m1/(m1-ymin))
cat('\n ' )

# Solve a cubic equation to find the roots of alpha
z <- c(1, -4, 5-(ymin^2)/m2, -2)
alpha_roots <- polyroot(z)        

# The first root is real whereas the other two roots are the complex conjugates
cat('\nPareto distribution (alpha)  = ',alpha_roots[1])    



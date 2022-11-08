#=========================================================================
#
#   Program to estimate a Cauchy model (one iteration)
#   using the NEWTON-RAPHSON, SCORING and BHHH algorithms
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

#
#--------------------------- Helper Functions -----------------------------------
# 

#-------------------------------------------------------------------------
# Wrapper function to calculate inverse of a given matrix
#-------------------------------------------------------------------------
inv <- function (M) {
  return(solve(M))
}

#
#--------------------------- Algorithms -----------------------------------
#

# Load data                     
yt <- c(2, 5, -2, 3, 3 )
t  <- length(yt)    

theta0 <- 3.0        

# Newton-Raphson
g <- 2*sum( (yt - theta0)/(1 + (yt - theta0)^2) )
h <- 2*sum( ((yt - theta0)^2 - 1)/(1 + (yt - theta0)^2)^2 )

thetaNR <- theta0 - inv(h) %*% g


# Method of Scoring
g <- 2*sum( (yt - theta0)/(1 + (yt - theta0)^2) )
i <- t/2

thetaSC <- theta0 + inv(i) %*% g


# BHHH
gt <- 2*( (yt - theta0)/(1 + (yt - theta0)^2) )

thetaBH <- theta0 + inv(gt %*% cbind(gt)) * sum(gt)


cat('\nNewton-Raphson     ' , thetaNR )
cat('\nMethod of Scoring  ', thetaSC )
cat('\nBHHH               ', thetaBH )


cat('\n' )
cat('\nGradient                       ', g)
cat('\nStandard error (Hessian)       ' , sqrt(-inv(h)) )
cat('\nStandard error (Information)   ',  sqrt( inv(i)) )
cat('\nStandard error (OPG)           ', sqrt(inv(gt %*% cbind(gt))) )




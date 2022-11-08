##########################################################################
##
##  Compute the steps of simplex and find the step thatis used 
##  	in first iteration.
##
##########################################################################

ch3_7 <- function()
{
  rm (list = ls(all=TRUE))
  graphics.off()
  
  # Read in the data for y
  y <- c(3.5, 1.0, 1.5)
  t <- length( y )         # Define the sample size
  
  th  <- c(1, 3)     # Initialize the estimates
  
  neglog <- lnl_neg(th,y,t)
  
  if(neglog[1] > neglog[2]) {
    th <- th[2:1]
    neglog <- neglog[2:1]    
  }  
  th_avg  <- mean(th)
  
  #   The three steps that may be followed
  
  #   Reflect
  alpha <- 0.5
  th_r  <- th_avg + alpha*(th_avg - th[2]) 
  
  #   Expand
  beta <- 1.1
  th_e <- th_avg + beta*(th_avg - th[2])
  
  #   Contract
  gamma <- 0.5
  th_c <- th_avg + gamma*(th_avg - th[2])
  
  # Decide which step is followed
  
  llen_r <- lnl_neg(th_r,y,t)
  
  if(llen_r < neglog[2]) {
    
    if(llen_r < neglog[1]) {
      cat('\nExpanding the simplex')
      llen_e <- lnl_neg(th_e,y,t)
      
      if (llen_e < llen_r) {
        cat('\nTheta2 is replaced by theta_e \n')
        th[2] <- th_e            
      }
      else {
        cat('\nTheta2 is replaced by theta_r \n')
        th[2] <- th_r
      }      
    }
    else {
      cat('\nReflecting the simplex \n')
      cat('\nTheta2 is replaced by theta_r \n')
      th[2] <- th_r   
    }      
  }       
  else {
    llen_c <- lnl_neg(th_c,y,t)
      if(llen_c < neglog[2]) {
        cat('\nContracting the simplex \n')
        cat('\nTheta2 is replaced by theta_c \n')
        th[2] <- th_c       
      }          
      else {
         fprintf('Shrinking the simplex \n')
          th <- (th + th_avg) / 2        
      } 
  }  
}


###########################################################################
#   SUBROUTINES
###########################################################################

# Negative log likelihood evaluator

lnl_neg <- function(theta,y,t){
  t*log(theta) + sum(y)/theta
}



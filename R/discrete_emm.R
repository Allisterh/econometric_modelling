#============================================================================
#
#   Sampling distribution of EMM estimator for the INARMA(1,1) model
#   derived using Monte Carlo methods.
#
#   THIS TAKES AN AGE TO RUN ... reduce niter, ndraws &/or n to get an idea
#
#============================================================================
rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")
#
#--------------------------- Helper Functions -------------------------------
#
# Load required functions -  inv, trimr, 
source("EMTSUtil.R")

#----------------------------------------------------------------------------
#  The objective function to compute the EMM estimator
#----------------------------------------------------------------------------
q <- function(b,n,lag,model,iinv,t,bhat) {
  set.seed(1, kind="Mersenne-Twister")
  rho_emm  <- b[1]
  lam_emm  <- b[2]
  beta_emm <- b[3]
  
  # Generate simulated data    
  us <- rpois(n, lambda=lam_emm)
  ys <- rpois(n, lambda=lam_emm*(1+beta_emm)/(1-rho_emm))
  
  zs <- ys                                                             
  
  for (i in 2:n) {
    # AR part
    if (zs[i-1] >0) {
      # Generate t_t-1 uniform random numbers    
      vs <- runif(zs[i-1])           
      # Sum bernoulli random draws  
      b_thin <- sum( vs < rho_emm )      
      # Generate realisation of y at time t   
      zs[i] <- b_thin + us[i] 
    } else {
      # There are no survivors
      zs[i] <- us[i]
    }
      
    # MA part
    if (us[i] >0  ) {
      # Generate uniform random numbers    
      vs <- runif(us[i],1)           
      # Sum bernoulli random draws  
      b_thin <- sum( vs < beta_emm )      
      # Generate realisation of y at time t   
      ys[i] <- b_thin + zs[i-1]
    } else {
      # There are no survivors
      ys[i] <- zs[i-1]
    }
  }
  
  if (lag == 1) {
    xs    <- cbind(rep(1, length(ys)-lag),   trimr(ys,0,1))
    ehats <- trimr(ys,1,0) - xs %*% bhat 
  }else if (lag == 2)
  {
    xs    <- cbind(rep(1, length(ys)-lag), trimr(ys,1,1),  trimr(ys,0,2))
    ehats <- trimr(ys,2,0) - xs %*% bhat    
  }else if (lag == 3) {
    xs    <- cbind(rep(1, length(ys)-lag), trimr(ys,2,1), trimr(ys,1,2), trimr(ys,0,3))
    ehats <- trimr(ys,3,0) - xs %*% bhat
  }  
  # Choose auxiliary model type      
  s2   <- mean(ehats^2)
  if (model == 1) {
    gs <- as.numeric(ehats) * xs
  } else if (model == 2) {
    gs <- as.matrix(c(as.numeric(ehats),   (ehats^2 - s2)))
  } 
  meang <- apply(gs, 2, mean)
  f  <- meang %*% iinv %*% cbind(meang)
  
  lf <- -0.5*t*f
  return(lf)
}


#-------------------------------------------------------------------------
#  Random search function
#-------------------------------------------------------------------------
search <- function(b0,n,lag,model,iinv,t,niter,bhat) {
  # Function evaluation at initial parameters
  f0 <- q(b0,n,lag,model,iinv,t,bhat)                           
  b <- rep(0, length(b0))
  
  for (i in seq(niter)) {
    b[1] <- runif(1)
    b[2] <- runif(1)*5
    b[3] <- runif(1)                             
    
    f <- q(b,n,lag,model,iinv,t,bhat)
    
    if (f < f0) {
      f0 <- f 
      b0 <- b
    }
  } 
 
  return(b0)
}


#
#---------------- EMM Estimator of Integer model  ---------------------------
#
discrete_emm <- function( ) {
  t      <- 50                                                     
  rho    <- 0.3                                                                  
  beta   <- 0.7
  lam    <- 3.5        
  theta0 <- c(rho,lam,beta)
  
  # Control length of simulation run
  h      <- 50     
  n      <- t*h             
  
  # Number of replications to generate finite sample properties
  ndraws <- 10
  
  # Maximum number of searches
  niter  <- 5                                                      
  
  # Loop over auxiliary models and lag length
  for (m in 1:2) {
    # Start lag at k<-2 for identification 
    for (k in 2:3) {
      model <- m
      lag   <- k
      
      # Main DO LOOP to generate sampling distribution
      theta_emm  <- array(0, c(ndraws,3))
      
      for (j in seq(ndraws)) {
        
        cat('\n# draws  = ',j ) 
        
        
        # Generate the actual data for binomial thinning model 
        # Generate Poisson random variables
        u <- rpois(t+100, lambda=lam)
        # Initialize y using the unconditional distribution 
        y <- rpois(t+100, lambda=lam*(1+beta)/(1-rho))
        z <- y                                                                                                                     
        for (i in 2:(t+100)) {
          # AR part
          if (z[i-1] > 0)
          {
            # Generate t_t-1 uniform random numbers    
            v <- runif(z[i-1],1)           
            # Sum bernoulli random draws  
            b_thin <- sum( v < rho )      
            # Generate realisation of y at time t   
            z[i] <- b_thin + u[i]
          } else {
            # There are no survivors
            z[i] <- u[i]
          }              
          # MA part
          if (u[i] > 0) {
            # Generate uniform random numbers    
            v <- runif(u[i],1)           
            # Sum bernoulli random draws  
            b_thin <- sum( v < beta )      
            # Generate realisation of y at time t   
            y[i] <- b_thin + z[i-1]
          }else {
            # There are no survivors
            y[i] <- z[i-1]
          }
        }
        # Trim first 100 observations to overcome startup problems
        y <- trimr(y,100,0)                    
        
        # Conditional least squares
        x    <- cbind(rep(1, length(y)-1), trimr(y,0,1))
        reg <- lm(trimr(y,1,0) ~ x - 1)
        bhat <- coefficients(reg)
        ehat <- residuals(reg)
        s2   <- mean(ehat^2)
        
        theta_0 <- c(bhat[2],  bhat[1],  0.1)
        
        # Estimate the auxiliary model using actual data  
        if (lag == 1) {
          x    <- cbind(rep(1, length(y)-lag),trimr(y,0,1))
          reg <- lm(trimr(y,1,0) ~ x - 1)
          bhat <- coefficients(reg)
          ehat <- residuals(reg)
          s2   <- mean(ehat^2)
        } else if (lag == 2){
          x    <- cbind(rep(1, length(y)-lag), trimr(y,1,1), trimr(y,0,2))
          reg <- lm(trimr(y,2,0) ~ x - 1)
          bhat <- coefficients(reg)
          ehat <- residuals(reg)
          s2   <- mean(ehat^2)              
        } else if (lag == 3) {
          x    <- cbind(rep(1, length(y)-lag), trimr(y,2,1),  trimr(y,1,2),  trimr(y,0,3))
          reg <- lm(trimr(y,3,0) ~ x - 1)
          bhat <- coefficients(reg)
          ehat <- residuals(reg)
          s2   <- mean(ehat^2)              
        }
        # Choose auxiliary model type          
        if (model == 1) {
          g <- ehat * x
          rowsg <- nrow(g)
        } else if (model == 2) {
          g <- as.matrix(c(ehat * x, (as.numeric(ehat)^2 - s2)))
        }
        
        # Compute the optimal weighting matrix      
        i <- t(g) %*% g
        p <- 0.0        
     
        for (l in seq(p )) {
          gam <- t( g[(l+1):rowsg,] ) %*% g[1:(rowsg-l),] 
          i   <- i + (1.0 - l/(p+1))*(gam + t(gam))        
        }
        i <- i/rowsg
        iinv <- inv(i)		
        
        theta_emm[j,] <- search(abs(theta_0),n,lag,model,iinv,t,niter,bhat)
      }       
      
      # Sampling distribution of the emm estimator   
      mm    <- colMeans(theta_emm)
      stdev <- sqrt((theta_emm - mm)^2)
      rmse  <- sqrt((theta_emm - theta0)^2)     
      
      cat('\n')
      cat('\nEMM estimator results') 
      cat('\n------------------------------------------------')
      cat('\nNumber of replications              = ',ndraws)
      cat('\nSample size                         = ',t) 
      cat('\nAuxiliarly model                    = ',model)
      cat('\nLags                                = ',lag) 
      cat('\nNumber of searches                  = ',niter)
      cat('\nNumber of simulation paths          = ',h, '\n')      
      print(list(True=theta0, Mean=t(mm), Std=colMeans(stdev), RMSE=colMeans(rmse)))
      cat('\n')
    }
  }
}



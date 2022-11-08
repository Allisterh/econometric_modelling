#=========================================================================
#
#   Program to compare the edgeworth expansion of the finite sample distribution
#   and the asymptotic distribution (exponential distribution) 
#
#=========================================================================

rm (list = ls(all=TRUE))
graphics.off()

t <- 5           
s <- seq(-2, 2, 1)   

# Finite sample value (based on a gamma distribution)
f_finite <- 1 - pgamma( t/(1 + s/sqrt(t)),t,1 )             

# Asymptotic value (based on the normal distribution)
f_asy <- pnorm( s )                

# Compute Edgeworth expansion values 
h2 <- s^2 - 1               
h3 <- s^3 - 3*s
h5 <- s^5 - 10*s^3 + 15*s

f_ew <- f_asy - dnorm(s) * ( (1 + (2/3)*h2)/sqrt(t) + ( (5/2)*s + (11/12)*h3 + (2/9)*h5 )/t )   

# Print results
cat('\nSample size    = ' ,t ,'\n')

cat('\n' )
rnames <- seq(t)
cnames <- c("s", "Finite", "Edgeworth", "Asymptotic")
print(matrix(c(s, f_finite, f_ew, f_asy), 
             dimnames=list(rnames, cnames), nrow=t))


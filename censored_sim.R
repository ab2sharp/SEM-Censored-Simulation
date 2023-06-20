library(ggplot2)
library(ggthemes)
library(ExtDist)


#Simulation parameters
n = 1000  
t = 2.4     

#Generate data
set.seed( 90053 )
y = rexp(n, rate=0.5)
x = y[ y < t]
r = sum( y >= t)



# Likelihood function -----------------------------------------------------

create.loglik <- function(x=NULL, r=NULL, n=NULL, t=NULL ) {
  ## for every value of theta 
  ## it evaluates the log-liklihood
  function(theta) {
    sumx = sum(x)
    pi0 = 1 - exp(-t/theta) 
    val = -length(x)*log(theta) - sumx/theta  + r*log(1-pi0)
    return(val)
  }
}

lik =  create.loglik(x=x, r=r, n=n, t=t)


create.qfn = function(x=NULL, r=NULL, t=NULL)
{
  meanx = mean(x)
  n = length(x) + r
  
  function(theta, th0){
    val = -n*log(theta) - (length(x)*meanx)/theta - ( r*(th0+t) )/theta
    return(val)
  }
}

qfn = create.qfn(x=x, r=r, t=t)


# Simulation --------------------------------------------------------------


#max vs mean with SEM

mm = 1e4; #Number of simulations
m =1e3; #Number of iterations of the SEM algorithm

sim.theta = matrix(0, nrow=mm, ncol=m)
sim.obs.loglik = matrix(0, nrow=mm, ncol=m)
theta0 = mean(y)

for (k in 1:mm) {
  sim.theta[k,1] = theta0
  sim.obs.loglik[k,1] = lik( sim.theta[k,1] )
  for (i in 2:m) {
    ### E-step
    y1 = t + rexp( r, rate=1/sim.theta[k,i-1] )
    ### M-step
    sim.theta[k,i] = ( sum(x) + sum(y1) )/ n
    sim.obs.loglik[k,i] = lik( sim.theta[k,i] )
  }
}

theta.avg = Rfast::rowmeans(sim.theta)
theta.max = sapply(1:k, function(k)
{ 
  sim.theta[k,which.max(sim.obs.loglik[k,]) ] 
})

theta.est = cbind(theta.avg, theta.max)
apply(theta.est, 2, mean)
apply(theta.est, 2, sd)


# Plot distributions ------------------------------------------------------

med = median(theta.max)
am = mean(abs(theta.max-med))


par(mfrow=c(1,3))
hist(theta.avg, breaks="FD", main="Tail Average", 
     xlim=c(1.978, 1.988), xlab=expression(theta), cex.lab=1.5)
hist(theta.max, breaks="FD", main="Minimum Likelihood Ratio",
     xlim=c(1.978, 1.988), xlab=expression(theta), cex.lab=1.5,
     ylab="")
hist(theta.max, breaks="FD", main="Minimum Likelihood Ratio, Rescaled", 
     freq=FALSE, xlab=expression(theta), cex.lab=1.5, ylab="")
curve(dLaplace(x,mu=med,b=am), add=TRUE, col="red")

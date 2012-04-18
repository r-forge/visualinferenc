pdf(file="Figure7.pdf",width=10)

# DEMO OF AN INCONSISTENT MAXIMUM LIKELIHOOD ESTIMATOR WITH CONTINUOUS PARAMETER
#
# Data is i.i.d. real, with distribution (1/2)N(0,1) + (1/2)N(t,exp(-1/t^2)^2),
# where t is a positive real parameter.


# PRODUCE DEMO PLOTS.  Leaves the generated data in the global variable x for
# possible later examination.

demo <- function()
{ 
  par(mfrow=c(2,3))

  plot.density(2.5)
  plot.density(0.6)
  plot.density(0.2)

  set.seed(1)
  x <<- gen(0.6,100)
  plot.lik(x[1:10])
  plot.lik(x[1:30])
  plot.lik(x)
}


# GENERATE DATA FROM THE MODEL WITH A GIVEN PARAMETER VALUE.  Arguments are
# the parameter value, t, and the number of data points to generate, m.

gen <- function (t,n)
{
  m <- rep(0,n)
  s <- rep(1,n)
  
  w <- runif(n) < 1/2
  m[w] <- t
  s[w] <- exp(-1/t^2)

  rnorm(n,m,s)
}


# COMPUTE LOG DENSITIES OF DATA VALUES FOR A GIVEN PARAMETER VALUE.  Arguments
# are a vector of data values, x, and the parameter value, t.  
#
# Special care is taken to compute the log of extreme density values without 
# overflow.  Data points exactly equal to the parameter value are treated
# specially, since for them the log of the exponential part of the density 
# vanishes, producing results that don't overflow even when the standard
# deviation, exp(-1/t^2), overflows.

log.density <- function (x,t)
{
  ll1 <- dnorm(x,0,1,log=TRUE) + log(0.5)

  ll2 <- dnorm(x,t,exp(-1/t^2),log=TRUE) + log(0.5)
  ll2[x==t] <- -0.5*log(2*pi) + 1/t^2 + log(0.5)

  ll <- rep(NA,length(x))
  for (i in 1:length(x))
  { ll[i] <- add.logs(ll1[i],ll2[i])
  }

  ll
}


# COMPUTE THE LOG LIKELIHOOD GIVEN A DATA VECTOR AND PARAMETER VALUE.  Arguments
# are the vector of data values, x, and the parameter value, t.

log.lik <- function (x,t)
{
  sum(log.density(x,t))
}


# PLOT THE DENSITY FUNCTION FOR A GIVEN PARAMETER VALUE.  Arguments are the
# parameter value, t, and the grid of data values at which to evaluate the
# density.  The grid defaults to 200 points from -3 to 5 (or the parameter
# value, if greater), plus the value of the parameter, where there may be a
# narrow peak in the density that would be missed by the grid.  The vertical
# scale is scaled by a power of ten close to the maximum density, since the
# density values might otherwise overflow.

plot.density <- function (t, grid=sort(c(t,seq(-3,max(5,t),length=200))))
{
  ld <- log.density(grid,t)
  max.ld <- round(max(ld)/log(10))
  ld <- ld - max.ld*log(10)

  plot(grid,exp(ld),type="l",
    xlab="data value",
    ylab=paste("density (x 10^",max.ld,")",sep=""))

  title(paste("Density with parameter",t))
}


# PLOT THE LIKELIHOOD FUNCTION FOR GIVEN DATA, AND FIND MLE.  The argument 
# is a vector of data values, x.  Additional arguments may be given, and are
# passed on the the plot funcition.  The likelihood is plotted for a grid of
# parameter values from 0.01 to 3, in steps of 0.01, plus all positive
# data values less than 3.  An approximation to the MLE is found by taking 
# the maximum over values on this grid.  Data values in (-3,3) are shown at 
# the bottom of the plot. The vertical scale is scaled by a power of ten close 
# to the maximum likelihood, since likelihood values might otherwise overflow.

plot.lik <- function (x,title="Log Likelihood",grid=NULL,...)
{
  if(is.null(grid)) grid <- sort(c(x[x>0 & x<3],seq(0.01,3,by=0.001)))

  ll <- rep(NA,length(grid))
  for (i in 1:length(grid))
  { ll[i] <- log.lik(x,grid[i])
  }

  mlv <- max(ll)
  mle <- grid[ll==mlv][1]
  max.ll <- round(mlv/log(10))

  #ll <- ll - max.ll*log(10)
  #lik <- exp(ll)
mlvi=(1:length(grid))[ll==mlv]
ll <- ll - ll[ll==mlv] + 2
lik <- ll

  plot(range(grid),c(-6,max(lik) + 4),type="n",xlim=c(-1,2),
    xlab="Data / Parameter",
    ylab=title,
    ...)
  points(x,rep(0,length(x)),lty=19)
  lines(grid,lik,lty=2,col=1)

  title(paste(title, "-",length(x),"data points, \n",
              "MLE approximately",round(mle,4)))
  iol=NULL
  for(k in 1:length(x)){
  for(ki in 1:length(grid))
  iol[ki]=log.lik(x[k],grid[ki]) 
  #iol=iol - iol[ll==mlv] + 2/length(x)
  lines(grid,iol - iol[mlvi] + 2/length(x),col=4)}

}


# ADD NUMBERS REPRESENTED BY LOGS.  Computes log(exp(a)+exp(b)) in a way
# that avoids overflow/underflow problems.

add.logs <- function (a,b)
{ if (a>b)
  { a + log(1+exp(b-a))
  }
  else
  { b + log(1+exp(a-b))
  }
}




demoKOR <- function(n,title="Log Likelihood",grid=NULL)
{ 
  set.seed(1)
  x <<- gen(0.6,n)
  plot.lik(x,title,grid=grid)
}

par(mfrow=c(1,2))

demoKOR(20,title="Log Continuous Likelihoods")


#demoKOR(15,sort(c(x[x>0.1 & x<.15],seq(0.13,.14,by=0.001))))

#change to interval
kdnorm=function(x,m,s,log)
log(pnorm(x + eps,m,s,log=FALSE) - pnorm(x - eps,m,s,log=FALSE))

eps=.01


log.density <- function (x,t)
{
  ll1 <- kdnorm(x,0,1,log=TRUE) + log(0.5)

  ll2 <- kdnorm(x,t,exp(-1/t^2),log=TRUE) + log(0.5)
  #ll2[x==t] <- -0.5*log(2*pi) + 1/t^2 + log(0.5)

  ll <- rep(NA,length(x))
  for (i in 1:length(x))
  { ll[i] <- add.logs(ll1[i],ll2[i])
  }

  ll
}

demoKOR(20,title="Log Observed Likelihoods")
dev.off()


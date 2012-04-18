pdf(file="Figure2.pdf",width=10)
par(mfrow=c(1,2))

#Profile likelihood
# simple regression example recast as a sum of individual log likelihoods

#usual regresion without an intercept
x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 17, 17, 17)
y = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 25, 25)


#write out the sum of individual log likleihoods
pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"], b * x[",i,"],s,log=T)",sep="")
}
pp<-paste("fn<-function(pars){  b=pars[1]; s=pars[2]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

#list log likelihood function
fn

#find the MLE and its SE by optimization
p0=c(1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)

#range of parameter values for b
ii=seq(0,3,(5 - -2)/1000)


#now reparameterize sd to ease optimization
pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],b * x[",i,"],exp(s),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ b=pars[1]; s=pars[2]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

p0=c(1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)

#recall from theory (parameter is exp(s))
log(sqrt(summary(glm(y~-1 + x))$disp))
log(sqrt(summary(glm(y~-1 + x))$disp) * sqrt((13 - 1)/13))

#record the joint profile path
proflik.path=function(b){

pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"], b * x[",i,"],exp(s),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ b=b; s=pars[1]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))


oo=optim(.9 * dmle[-2], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))
return(c(-oo$value,oo$par))
}

#check that the mle is correct
proflik.path(dmle[1])
oo$par
oo$value

#calculate and plot the combined log likelihood
cc=matrix(NA,ncol=2,nrow=length(ii))
for(j in 1:length(ii))
cc[j,]=proflik.path(ii[j])
mle=min(abs(ii-dmle[1])) == abs(ii-dmle[1])
plot(ii,cc[,1] - cc[mle,1] + 0 * 2 ,type="n",ylim=c(-4,0 * 4),ylab="Log Likelihood",xlab="Beta")
title(main="Profile Marginalization to Beta Parameter")
lines(ii,cc[,1] - cc[mle,1] + 0 * 2,col=4,lty=2,lwd=2)

#mark the mle for slope
abline(v=dmle[1])

#individual likelihoods

#matrix to hold results
il=matrix(NA,nrow=length(y),ncol=length(ii))

#write out, evaluate and plot the individual log likelihoods
for(j in 1:length(y)){
pp<-paste("dnorm(y[",j,"], b * x[",j,"],exp(s),log=T)",sep="")
pp<-paste("fn<-function(pars){ b=pars[1]; s=pars[2]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
ic=NULL
for(k in 1:length(ii))
ic[k]=-fn(c(ii[k],cc[k,-1]))
il[j,]=ic - max(ic)
lines(ii,il[j,] + 0 * rnorm(1,0,.05),col=x[j],lwd=2)
#text(ii,il[j,],j)
}


#emphasize the "outliers" and combined log likelihoods
lines(ii,il[j,] + 0 * rnorm(1,0,.05),col=x[j],lwd=4)
lines(ii,cc[,1] - cc[mle,1] + 0 * 2,col=4,lty=2,lwd=3)




#Integrtaed log likelihood
library(pscl)

#Write out single observation likelihoods (not on log scale)
pp<-""
for(i in 1:length(y)){
pp<-paste("dnorm(y[",i,"],b1 * x[",i,"],sqrt(sig.sq),log=F)",sep="") 
pp<-paste("datamodel",i,"<-function(pars){ b1=pars[1]; sig.sq=pars[2]; return((",pp,"))}",sep="")
eval(parse(text=pp))
}
datamodel1
datamodel13

#Write out prior
pp<-""
pp<-paste(pp," 
dnorm(b1,u.0,tau.0) * densigamma(sig.sq,v.0/2,(v.0 * sig.sq.0)/2) ",sep=" ")
pp<-paste("prior<-function(pars){u=pars[1];b1=pars[2];sig.sq=pars[3];return((",pp,"))}",sep="")
eval(parse(text=pp))

prior 

#set prior parameters
u.0=0
tau.0=1/2
v.0=4
sig.sq.0=4


#Combine idividual observation likelihoods by multiplication
pp<-""
for(i in 1:(length(y) - 1)) {
pp<-paste(pp," 
datamodel",i,"(pars)"," * ",sep="")
}
i=length(y)
pp<-paste(pp," 
datamodel",i,"(pars)",sep="")
pp<-paste("datamodel<-function(pars){ b1=pars[1]; sig.sq=pars[2]; return((",pp,"))}",sep="")
eval(parse(text=pp)) 

datamodel 

#Focus on b1 parameter
b1=seq(-0,3,by=1/80) 

b1.prior=dnorm(b1,mean=u.0,tau.0,log=T) 
plot(b1,b1.prior - max(b1.prior),ylim=c(-4,0),type="n",ylab="Log Prior Likelihood Posterior",xlab="Beta") 
title(main="Integrated Marginalization to Beta Parameter")
lines(b1,b1.prior - max(b1.prior),col=2,lwd=2,lty=3) 


#re-write out the sum of individual likelihoods - note b1 will be a vector rather than a function parameter
pp<-""
for(i in 1:length(y)){
pp<-paste("dnorm(y[",i,"],b1 * x[",i,"],sqrt(sig.sq),log=F)",sep="") 
pp<-paste("datamodel",i,"<-function(pars){sig.sq=pars[1]; return((",pp,"))}",sep="")
eval(parse(text=pp))
} 
pp<-""
for(i in 1:(length(y) - 1)) {
pp<-paste(pp," 
datamodel",i,"(pars)"," * ",sep="")
}
i=length(y)
pp<-paste(pp," 
datamodel",i,"(pars)",sep="")
pp<-paste("fn<-function(pars){sig.sq=pars[1]; return((",pp,"))}",sep="")
eval(parse(text=pp)) 

fn
datamodel



#draw simple random sample of sig.sq from its prior (independent of b1 by assumptions), average and plot log of it
foos=function(n1,n2){
ff=fn(1) - fn(1)
for(k in 1:n1){
n=n2
rs=rigamma(n,v.0/2,(v.0 * sig.sq.0)/2)
rp=cbind(rs,rs)
ffi=apply(apply(rp,1,fn),1,sum)
ff=ffi + ff 
}
lines(b1,log(ff/max(ff,na.rm=T)),col=4,lty=2,lwd=2) 
return(ff)
}

nloop=20;ndraw=1000
for(j in 1:1){
ff=foos(nloop,ndraw)
}


#Verify the posterior using MCMC
library(MCMCpack)
line <- list(X1 = x,Y = y) 
for(I in 1:1) {
posterior <- MCMCregress(Y~-1 + X1, data=line, seed=trunc(runif(1) *10000),mcmc=10^6,c0=v.0,d0= v.0 * sig.sq.0,B0=1/tau.0^2) 
beta1=posterior[,1]
dd=density(beta1,from=min(b1),to=max(b1),n=length(b1))
lines(b1,log(dd$y) - max(log(dd$y)),col=6,lty=4,lwd=2)

lines(b1,log(dd$y) - b1.prior - max(log(dd$y)-b1.prior),col=4,lty=2,lwd=2)
}


#Obtain the profile likelihhods to be deformed to add to combined integrated likelihood
ii=b1
dmle=glm(y~-1 + x)$coef


#record the joint profile path
proflik.path=function(b){

pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],b * x[",i,"],exp(s),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ b=b; s=pars[1]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))


oo=optim(c(1,1), fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))
return(c(-oo$value,oo$par))
}


#check that the mle is correct
proflik.path(dmle[1])
oo$par
oo$value


#calculate the combined profile log likelihood
cc=matrix(NA,ncol=3,nrow=length(ii))
for(j in 1:length(ii))
cc[j,]=proflik.path(ii[j])
mle=min(abs(ii-dmle[1])) == abs(ii-dmle[1])
lines(ii,cc[,1] - cc[mle,1],col=3,lty=2,lwd=2)


#ratio of integrated to profile
intlRpl=log(ff)/cc[,1]

#individual likelihoods over the profile path
ill=matrix(NA,ncol=length(b1),nrow=length(y))

for(j in 1:length(y)){
pp<-paste("dnorm(y[",j,"], b * x[",j,"],sqrt(sig.sqv),log=T)",sep="")
pp<-paste("fni<-function(pars){ b=pars[1]; sig.sqv=pars[2]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
ic=NULL
for(k in 1:length(b1))
ic[k]=-fni(c(b1[k],exp(cc[k,2])^2))
offset=runif(1,-.1,.1)
lines(b1,intlRpl * (ic - max(ic)) + 0 * offset,col=x[j],lwd=2)
#text(b1[round(b1/.25,2) - round(b1/.25,0) == 0],(intlRpl * (ic - max(ic)) + offset)[round(b1/.25,2) - round(b1/.25,0) == 0],letters[j],cex=.75)
#lines(b1,(ic - max(ic)) + runif(1,0,.1) ,col=4)
ill[j,]=intlRpl * ic
}

#emphasize the "outliers"
lines(b1,intlRpl * (ic - max(ic)) + 0 * offset,col=x[j],lwd=4)

#check that these sum correctly to combined log likleihood - not run
#cpl=apply(ill,2,sum)
#text(b1,cpl - max(cpl),col=4,"s")

#Plot combined log posterior, integrated log likelihood and profile log likelihood 
lines(b1,log(dd$y) - max(log(dd$y)),col=6,lty=4,lwd=3)
lines(b1,log(dd$y) - b1.prior - max(log(dd$y)-b1.prior),col=4,lty=2,lwd=3)
lines(ii,cc[,1] - cc[mle,1],col=3,lty=2,lwd=3)
dev.off()
pdf(file="Figure1.pdf",width=10)
#set for two plots
par(mfrow=c(1,2))

#Estimated log likelihood plot

# simple regression example recast as a sum of individual log likelihoods

#usual regresion without an intercept
x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 17, 17, 17)
y = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 25, 25)

#linear model without intercept
glm1=glm(y~-1 + x)

#write out the sum of individual log likleihoods
pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"], b * x[",i,"],s,log=T)",sep="")
}
pp<-paste("fn<-function(pars){  b=pars[1];  return(-(",pp,"))}",sep="")
eval(parse(text=pp))



#list the log likelihood function
fn


#set s equal to the usual REML estimate and treat as known
s=sqrt(summary(glm1)$disp)

#set mle for slope
mle.b=summary(glm1)$coef[1,1]
#range of parameter values for b
ii=seq(0,3,(5 - -2)/1000)

#calculate and plot the combined log likelihood
cc=NULL
for(j in 1:length(ii))
cc[j]=-fn(ii[j])
mle=min(abs(ii-mle.b)) == abs(ii-mle.b)
plot(ii,cc - cc[mle] + 0 * 2 ,type="n",ylim=c(-4,0 * 4),ylab="Log Likelihood",xlab="Beta")
title(main="Estimated Marginalization to Beta Parameter")
lines(ii,cc - cc[mle] + 0 * 2,col=4,lty=2,lwd=2)

#mark the mle  
abline(v=mle.b)

#write out, evaluate and plot the individual log likelihoods
for(j in 1:length(y)){
pp<-paste("dnorm(y[",j,"], b * x[",j,"],s,log=T)",sep="")
pp<-paste("fn<-function(pars){ b=pars[1]; s=s; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

ic=NULL
for(k in 1:length(ii))
ic[k]=-fn(ii[k])

i1=ic - max(ic)
lines(ii,i1 + 0 * rnorm(1,0,.05),col=x[j],ylim=c(-4,4),lwd=2)

}

#emphasize the "outliers" and combined log likelihoods
lines(ii,i1 + 0 * rnorm(1,0,.05),col=x[j],lwd=4)
lines(ii,cc - cc[mle] + 0 * 2,col=4,lty=2,lwd=3)


#now add an intercept (another parameter) and note not all log likelihoods are unimodal 
glm2=glm(y~x)

#write out the log likelihood
pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"], u + b * x[",i,"],s,log=T)",sep="")
}
pp<-paste("fn<-function(pars){  b=pars[1]; u=u; s=s; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

#list the log likelihood function
fn

#set u and s equal to the usual REML estimates and treat as known
u=summary(glm2)$coef[1,1]
s=sqrt(summary(glm2)$disp)
mle.b=summary(glm2)$coef[2,1]


#calculate and plot the combined log likelihood
cc=NULL
for(j in 1:length(ii))
cc[j]=-fn(ii[j])
mle=min(abs(ii-mle.b)) == abs(ii-mle.b)
plot(ii,cc - cc[mle] + 0 * 2 ,type="n",ylim=c(-4,0 * 4),xlim=c(0,3),ylab="Log Likelihood",xlab="Beta")
title(main="Estimated Marginalization to Beta Parameter")
lines(ii,cc - cc[mle] + 0 * 2,col=4,lty=2)

#mark mle for slope  
abline(v=mle.b)

#write out, evaluate and plot the individual log likelihoods
for(j in 1:length(y)){
pp<-paste("dnorm(y[",j,"], u + b * x[",j,"],s,log=T)",sep="")
pp<-paste("fn<-function(pars){ b=pars[1]; u=u; s=s; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
ic=NULL
for(k in 1:length(ii))
ic[k]=-fn(ii[k])

i1=ic - max(ic)
lines(ii,i1 + 0 * rnorm(1,0,.05),col=x[j],ylim=c(-4,4),lwd=2)

}


#emphasize the "outliers" and combined log likelihoods
lines(ii,i1 + 0 * rnorm(1,0,.05),col=x[j],ylim=c(-4,4),lwd=4)
lines(ii,cc - cc[mle] + 0 * 2,col=4,lty=2,lwd=3)
dev.off()


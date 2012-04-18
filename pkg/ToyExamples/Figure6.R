pdf(file="Figure6.pdf",width=10)
par(mfrow=c(1,2))
#need this library
library(MASS)

#Make it reproducible
set.seed(1)

#In simulation use published bivariate group means and thier variances and covariances
d<- c( 0.47 , -.32 , 0.2 , -0.6 , 0.4 , -.12 , 0.26 , -.31 , .56 , -.39 )

s1<-c(0.0075466 ,0.0030361 , 0.0030361 , 0.0077251)

s2<-c(0.005712 , 0.0009184 , 0.0009184 ,0.0008452)

s3<-c(0.0020799, 0.0007368 , 0.0007368, 0.0014451 )

s4<-c(0.0029269, 0.0009398 , 0.0009398 , 0.0014867)

s5<-c(0.0148493, 0.0072249 , 0.0072249 ,0.0304102)

sl<-list(s1,s2,s3,s4,s5)

s<-matrix(0,ncol=10,nrow=10)st

for(j in 1:5){

i<- 2 * j

s[(i-1):i,(i-1):i]<-matrix(sl[[j]],ncol=2)

}
 

y<-matrix(d,ncol=2,byrow=T)

s<-list()

for(i in 1:5)
s[[i]]<-matrix(sl[[i]],ncol=2)

nump<-c(14,15,78,89,16)



study<-list()

for(i in 1:2)
study[[i]]<-cbind(mvrnorm(n = nump[i], y[i,], nump[i] * s[[i]], tol = 1e-6, empirical = T),rnorm(nump[i]),rnorm(nump[i]))


#A standard parametrization
#N(c|u.cp, s.cp) * N(e|µ.ep+ßp(c – µ.cp), s.p(e|c)) { control } 
#N(c|u.ct, s.ct) * N(e|µ.et+ßt(c – µ.ct), s.t(e|c)) { treatment }

#With covariates x.e and x.c
#N(c|u.c + ß.c X.c,s.c) * 
#N(e|µ.e+ß(c – (µ.c + ß.c X.c)) + ß.e X.e,s.(e|c)) [Nixon]
#re-parameterization µ.et -> NMB/K + µ.ep + u.ct/K - u.cp/K 
# 
 



parlabs=c("u.cp","u.ep","ßp","s.cp","s.p(e|c)","u.ct","u.et","ßt","s.ct","s.t(e|c)","ßcr","ßxr")
 

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11] * study[[",i,"]][,3])) + pars[12] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))",

" + sum(dnorm(study[[",i + 1,"]][,1],pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3],pars[",(i-1)*3," + 4 + 5], log=T) + 

dnorm(study[[",i + 1,"]][,2],pars[2 + 5] + pars[",(i-1)*3," + 3 + 5] * (study[[",i + 1,"]][,1] - (pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3])) + pars[12] * study[[",i + 1,"]][,4],pars[",(i-1)*3," + 5 + 5] ,log=T))",sep="")

}

pp<-paste("fn<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))
fn


#get starting value

p0=c(apply(study[[1]],2,mean)[1:2],glm(study[[1]][,2] ~ study[[1]][,1])$coef[2],apply(study[[1]],2,sd)[1:2],apply(study[[2]],2,mean)[1:2],

glm(study[[2]][,2] ~ study[[2]][,1])$coef[2],apply(study[[2]],2,sd)[1:2],c(0,0))

#get mles

foo<-function(d){

abs(fn(p0) -fn(p0 + d))

}

sc<-NULL

for(i in 1:length(p0))

sc[i]<-foo( (1:length(p0) == i)/10)

sc

oo<-optim(p0,fn,method = "BFGS",hessian=T,control = list(parscale=1/sc))

rbind(parlabs,round(oo$par,3))


#re-parameterization µ.et -> NMB/K + µ.ep + u.ct/K - u.cp/K 

parlabs=c("u.cp","u.ep","ßp","s.cp","s.p(e|c)","u.ct","NMB","ßt","s.ct","s.t(e|c)","ßcr","ßxr")

 

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11] * study[[",i,"]][,3])) + pars[12] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))",

" + sum(dnorm(study[[",i + 1,"]][,1],pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3],pars[",(i-1)*3," + 4 + 5], log=T) + 

dnorm(study[[",i + 1,"]][,2],pars[2 + 5]/k + pars[2] + pars[6]/k - pars[1]/k + pars[",(i-1)*3," + 3 + 5] * (study[[",i + 1,"]][,1] - (pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3])) + pars[12] * study[[",i + 1,"]][,4],pars[",(i-1)*3," + 5 + 5] ,log=T))",sep="")

}

pp<-paste("fn<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

#Reparameterized log likelihood
fn

#get transformed starting value
p0[7]=(2 * (p0[7] - p0[2]) - (p0[6] - p0[1]))

#set value of K (trade off for NMB)
k=2

#get mles

foo<-function(d){

abs(fn(p0) -fn(p0 + d))

}

sc<-NULL

for(i in 1:length(p0))

sc[i]<-foo( (1:length(p0) == i)/10)

sc

oo<-optim(p0,fn,method = "BFGS",hessian=T,control = list(parscale=1/sc))

rbind(parlabs,round(oo$par,3))

#get mle for NMB

dmle<-oo$par[7]

dmle.se<-sqrt(diag(solve(oo$hessian)))[7]

 

#Choose range for NMB
ii<-seq(-1.5,.5,by=.1)
#

#with NMB parameterization optim over other parameters 

proflik<-function(NMB){

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11 -1] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11 -1] * study[[",i,"]][,3])) + pars[12 - 1] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))",

" + sum(dnorm(study[[",i + 1,"]][,1],pars[1 + 5] + pars[11 -1] * study[[",i + 1,"]][,3],pars[",(i-1)*3," + 4 + 5 -1], log=T) + 

dnorm(study[[",i + 1,"]][,2],NMB/k + pars[2] + pars[6]/k - pars[1]/k + pars[",(i-1)*3," + 3 + 5 -1] * (study[[",i + 1,"]][,1] - (pars[1 + 5] + pars[11 -1] * study[[",i + 1,"]][,3])) + pars[12 -1 ] * study[[",i + 1,"]][,4],pars[",(i-1)*3," + 5 + 5 -1] ,log=T))",sep="")

}

pp<-paste("fn<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

oo<-optim(oo$par[-7], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))

return(c(-oo$value,oo$par))

}

proflik(dmle)

cc<-matrix(NA,nrow=length(ii),ncol=12)

for(j in 1:length(ii))

cc[j,]<-proflik(ii[j])

mle=min(abs(ii-dmle[1])) == abs(ii-dmle[1])

plot(ii,cc[,1] - cc[mle,1] + 2 ,type="n",ylim=c(-2.5,3.5),xlab="Net Monetary Benefit",
ylab="Log Likelihood")
title("No Outliers")

lines(ii,cc[,1] - cc[mle,1] + 2,col=4,lty=2,lwd=2)

abline(v=dmle[1])

 

 

#individual full likelihood

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11] * study[[",i,"]][,3])) + pars[12] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))/nrow(study[[",i,"+1]])",

" + sum(dnorm(study[[",i + 1,"]][j,1],pars[1 + 5] + pars[11] * study[[",i + 1,"]][j,3],pars[",(i-1)*3," + 4 + 5], log=T) + 

dnorm(study[[",i + 1,"]][j,2],pars[2 + 5]/k + pars[2] + pars[6]/k - pars[1]/k + pars[",(i-1)*3," + 3 + 5] * (study[[",i + 1,"]][j,1] - (pars[1 + 5] + pars[11] * study[[",i + 1,"]][j,3])) + pars[12] * study[[",i + 1,"]][j,4],pars[",(i-1)*3," + 5 + 5] ,log=T))",sep="")

}

pp<-paste("fn1<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

ccpar=cbind(cc[,2:7],ii,cc[,8:12])

 

ic=array(NA,c(length(ii),nrow(study[[2]])))

icnorm=array(NA,c(length(ii),nrow(study[[2]])))

for(j in 1:nrow(study[[2]])){

for(l in 1:length(ii))

ic[l,j]=-fn1(ccpar[l,])

i1=ic[,j] - ic[mle,j] + 2/nrow(study[[2]])

lines(ii,i1,col=j + 1,lwd=2)

icnorm[,j]=i1

}

 
#Check that sum of individual profile log likelihoods equals the combined profile log likelihood
sumic=apply(icnorm,1,sum)

sumic

lines(ii,sumic,col=4,lty=2,lwd=2)


#make outlier and repeat above
study[[2]][1]=2.5


#standard parametrization

#N(c|u.cp, s.cp) * N(e|µ.ep+ßp(c – µ.cp), s.p(e|c)) { control } 

#N(c|u.ct, s.ct) * N(e|µ.et+ßt(c – µ.ct), s.t(e|c)) { treatment }

#par[1] = u.cp

parlabs=c("u.cp","u.ep","ßp","s.cp","s.p(e|c)","u.ct","u.et","ßt","s.ct","s.t(e|c)","ßcr","ßxr")

 

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11] * study[[",i,"]][,3])) + pars[12] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))",

" + sum(dnorm(study[[",i + 1,"]][,1],pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3],pars[",(i-1)*3," + 4 + 5], log=T) + 

dnorm(study[[",i + 1,"]][,2],pars[2 + 5] + pars[",(i-1)*3," + 3 + 5] * (study[[",i + 1,"]][,1] - (pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3])) + pars[12] * study[[",i + 1,"]][,4],pars[",(i-1)*3," + 5 + 5] ,log=T))",sep="")

}

pp<-paste("fn<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

#get starting value


p0=c(apply(study[[1]],2,mean)[1:2],glm(study[[1]][,2] ~ study[[1]][,1])$coef[2],apply(study[[1]],2,sd)[1:2],apply(study[[2]],2,mean)[1:2],

glm(study[[2]][,2] ~ study[[2]][,1])$coef[2],apply(study[[2]],2,sd)[1:2],c(0,0))

#get mles

foo<-function(d){

abs(fn(p0) -fn(p0 + d))

}

sc<-NULL

for(i in 1:length(p0))

sc[i]<-foo( (1:length(p0) == i)/10)

sc

oo<-optim(p0,fn,method = "BFGS",hessian=T,control = list(parscale=1/sc))

rbind(parlabs,round(oo$par,3))

 

#re-parameterization µ.et -> NMB/K + µ.ep + u.ct/K - u.cp/K 

parlabs=c("u.cp","u.ep","ßp","s.cp","s.p(e|c)","u.ct","NMB","ßt","s.ct","s.t(e|c)","ßcr","ßxr")

 

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11] * study[[",i,"]][,3])) + pars[12] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))",

" + sum(dnorm(study[[",i + 1,"]][,1],pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3],pars[",(i-1)*3," + 4 + 5], log=T) + 

dnorm(study[[",i + 1,"]][,2],pars[2 + 5]/k + pars[2] + pars[6]/k - pars[1]/k + pars[",(i-1)*3," + 3 + 5] * (study[[",i + 1,"]][,1] - (pars[1 + 5] + pars[11] * study[[",i + 1,"]][,3])) + pars[12] * study[[",i + 1,"]][,4],pars[",(i-1)*3," + 5 + 5] ,log=T))",sep="")

}

pp<-paste("fn<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

p0[7]=(2 * (p0[7] - p0[2]) - (p0[6] - p0[1]))

#set value of K

k=2

#get mles

foo<-function(d){

abs(fn(p0) -fn(p0 + d))

}

sc<-NULL

for(i in 1:length(p0))

sc[i]<-foo( (1:length(p0) == i)/10)

sc

oo<-optim(p0,fn,method = "BFGS",hessian=T,control = list(parscale=1/sc))

rbind(parlabs,round(oo$par,3))

#get mle for NMB

dmle<-oo$par[7]

dmle.se<-sqrt(diag(solve(oo$hessian)))[7]

 

#with NMB parameterization optim over other parameters 

proflik<-function(NMB){

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11 -1] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11 -1] * study[[",i,"]][,3])) + pars[12 - 1] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))",

" + sum(dnorm(study[[",i + 1,"]][,1],pars[1 + 5] + pars[11 -1] * study[[",i + 1,"]][,3],pars[",(i-1)*3," + 4 + 5 -1], log=T) + 

dnorm(study[[",i + 1,"]][,2],NMB/k + pars[2] + pars[6]/k - pars[1]/k + pars[",(i-1)*3," + 3 + 5 -1] * (study[[",i + 1,"]][,1] - (pars[1 + 5] + pars[11 -1] * study[[",i + 1,"]][,3])) + pars[12 -1 ] * study[[",i + 1,"]][,4],pars[",(i-1)*3," + 5 + 5 -1] ,log=T))",sep="")

}

pp<-paste("fn<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

oo<-optim(oo$par[-7], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))

return(c(-oo$value,oo$par))

}

proflik(dmle)

cc<-matrix(NA,nrow=length(ii),ncol=12)

for(j in 1:length(ii))

cc[j,]<-proflik(ii[j])

mle=min(abs(ii-dmle[1])) == abs(ii-dmle[1])

plot(ii,cc[,1] - cc[mle,1] + 2 ,type="n",ylim=c(-2.5,3.5),xlab="Net Monetary Benefit",
ylab="Log Likelihood")
title("Just One Outlier")

lines(ii,cc[,1] - cc[mle,1] + 2,col=4,lty=2,lwd=2)

abline(v=dmle[1])

 

 

#individual full likelihood

pp<-""

for(i in 1:1){

pp<-paste(pp," + 

sum(dnorm(study[[",i,"]][,1],pars[1] + pars[11] * study[[",i,"]][,3],pars[",(i-1)*3," + 4], log=T) + 

dnorm(study[[",i,"]][,2],pars[2] + pars[",(i-1)*3," + 3] * (study[[",i,"]][,1] - (pars[1] + pars[11] * study[[",i,"]][,3])) + pars[12] * study[[",i,"]][,4],pars[",(i-1)*3," + 5] ,log=T))/nrow(study[[",i,"+1]])",

" + sum(dnorm(study[[",i + 1,"]][j,1],pars[1 + 5] + pars[11] * study[[",i + 1,"]][j,3],pars[",(i-1)*3," + 4 + 5], log=T) + 

dnorm(study[[",i + 1,"]][j,2],pars[2 + 5]/k + pars[2] + pars[6]/k - pars[1]/k + pars[",(i-1)*3," + 3 + 5] * (study[[",i + 1,"]][j,1] - (pars[1 + 5] + pars[11] * study[[",i + 1,"]][j,3])) + pars[12] * study[[",i + 1,"]][j,4],pars[",(i-1)*3," + 5 + 5] ,log=T))",sep="")

}

pp<-paste("fn1<-function(pars){ return(-(",pp,"))}",sep="")

eval(parse(text=pp))

ccpar=cbind(cc[,2:7],ii,cc[,8:12])

 

ic=array(NA,c(length(ii),nrow(study[[2]])))

icnorm=array(NA,c(length(ii),nrow(study[[2]])))

for(j in 1:nrow(study[[2]])){

for(l in 1:length(ii))

ic[l,j]=-fn1(ccpar[l,])

i1=ic[,j] - ic[mle,j] + 2/nrow(study[[2]])

lines(ii,i1,col=j + 1,lwd=2)

icnorm[,j]=i1

}

 

sumic=apply(icnorm,1,sum)

sumic

lines(ii,sumic,col=4,lty=2,lwd=2)
dev.off()


pdf(file="Figure4.pdf",width=10)
par(mfrow=c(1,2))
#colours
library(RColorBrewer)
mypalette.M<-brewer.pal(9,"Blues")[-1]
mypalette.F<-brewer.pal(9,"RdPu")[c(-1,-9)]
colourl=list(mypalette.M,mypalette.F)

myletters=c(letters,"@","#","$","%")
#Generate data according to the model with an intercation being present
set.seed(1)
n=30
x1=rnorm(n)
x2=rbinom(n,1,.5)
y=rnorm(n,mean=x1 - .5 * x1 * x2)

#Fit and display results from linear models
coefc=glm(y~x1 * x2)$coef
coef1=glm(y[x2==1] ~ x1[x2==1])$coef
coef0=glm(y[!x2] ~ x1[!x2])$coef

coefc
coef0
c(coefc[1] + coefc[3],coefc[2] + coefc[4])
coef1

ssc=summary(glm(y~x1 * x2))
ss1=summary(glm(y[x2==1] ~ x1[x2==1]))
ss0=summary(glm(y[x2==0] ~ x1[x2==0]))
ssc$disp
ss1$disp; ss0$disp

ssc$disp * ssc$df[2]
ss1$disp * ss1$df[2] + ss0$disp * ss0$df[2]

#Write out combined likelihood
pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] + b2 * x2[",i,"]+ b3 * x1[",i,"] * x2[",i,"],exp(s),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ u=pars[1]; b1=pars[2]; b2=pars[3]; b3=pars[4];s=pars[5]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
fn

#Obtain joint mle
p0=c(1,1,1,1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)

#recall from theory
sqrt(ssc$disp) * sqrt((30 - 4)/30)
exp(dmle[5])

#Write out log combined log likelihood from just females (x2==1)
pp<-""
for(i in 1:length(y)){
if(x2[i]==1){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] ,exp(s),log=T)",sep="")
}
}
pp<-paste("fn<-function(pars){ u=pars[1]; b1=pars[2]; s=pars[3]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
fn

p0=c(1,1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)

#recall from theory
k=sum(x2==1)
sqrt(ss1$disp) * sqrt((k - 2)/k)
exp(dmle[3])

#Now for just males (x2=0)
pp<-""
for(i in 1:length(y)){
if(x2[i]==0){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] ,exp(s),log=T)",sep="")
}
}
pp<-paste("fn<-function(pars){ u=pars[1]; b1=pars[2]; s=pars[3]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

p0=c(1,1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)

#recall from theory
k=sum(x2==0)
sqrt(ss0$disp) * sqrt((k - 2)/k)
exp(dmle[3])

#back to full model

pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] + b2 * x2[",i,"]+ b3 * x1[",i,"] * x2[",i,"],exp(s),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ u=pars[1]; b1=pars[2]; b2=pars[3]; b3=pars[4];s=pars[5]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
fn

p0=c(1,1,1,1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)



#record the joint profile path for b3
proflik.path=function(b1){


pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] + b2 * x2[",i,"]+ b3 * x1[",i,"] * x2[",i,"],exp(s),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ b1=b1; u=pars[1]; b2=pars[2]; b3=pars[3]; s=pars[4]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))


oo=optim(.9 * dmle[-4], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))
return(c(-oo$value,oo$par))
}

#Check mle
proflik.path(dmle[2])
oo$par
oo$value

#Calculate and plot combined log likelihood
ii=seq(-1,3,.1 * 1)
cc=matrix(NA,ncol=5,nrow=length(ii))
for(j in 1:length(ii))
cc[j,]=proflik.path(ii[j])
mle=min(abs(ii-dmle[2])) == abs(ii-dmle[2])
plot(ii,cc[,1] - cc[mle,1] + 2 ,type="n",ylim=c(-3,3),ylab="Log Likelihood",xlab="Male Slope Coefficient")
lines(ii,cc[,1] - cc[mle,1] + 2,lty=2,lwd=2,col=6)

abline(h=0,col=8)  
abline(v=dmle[2],col=8)

#individual likelihoods
il=matrix(NA,nrow=length(y),ncol=length(ii))

for(j in 1:length(y)){
pp<-""
pp<-paste(pp," 
+ dnorm(y[",j,"],u + b1 * x1[",j,"] + b2 * x2[",j,"]+ b3 * x1[",j,"] * x2[",j,"],exp(s),log=T)",sep="")

pp<-paste("fn<-function(pars){ b1=pars[1]; u=pars[2]; b2=pars[3]; b3=pars[4]; s=pars[5]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

ic=NULL
for(k in 1:length(ii))
ic[k]=-fn(c(ii[k],cc[k,-1]))

il[j,]=ic - ic[mle] + 2/length(y)
lines(ii,il[j,] + rnorm(1,0,.025),col=colourl[[(1 + x2[j])]][rep(1:7,5)[j]],lty=c(1,3)[x2[j] + 1],lwd=2)

#text(ii,il[j,],myletters[j],col=colourl[[(1 + x2[j])]][rep(1:7,5)[j]])
}

#Plot combined seperately by male and female log likelihoods
lines(ii,apply(il[x2==0,],2,sum),col=4,lty=2,lwd=2)
lines(ii,apply(il[x2==1,],2,sum),col=2,lty=2,lwd=2)

title("When forced to be as variable as men \n - women have varied ideas about the male slope")



pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] + b2 * x2[",i,"]+ b3 * x1[",i,"] * x2[",i,"],
exp(s1 * (x2[",i,"] == 0) + s2 * (x2[",i,"] == 1)),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ u=pars[1]; b1=pars[2]; b2=pars[3]; b3=pars[4];s1=pars[5];s2=pars[6]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))
fn

p0=c(1,1,1,1,1,1)
oo<-optim(p0,fn,method = "BFGS",hessian=T)
dmle<-oo$par
dmle.se<-sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)

#record the joint profile path for b3
proflik.path=function(b1){


pp<-""
for(i in 1:length(y)){
pp<-paste(pp," 
+ dnorm(y[",i,"],u + b1 * x1[",i,"] + b2 * x2[",i,"]+ b3 * x1[",i,"] * x2[",i,"],
exp(s1 * (x2[",i,"] == 0) + s2 * (x2[",i,"] == 1)),log=T)",sep="")
}
pp<-paste("fn<-function(pars){ b1=b1; u=pars[1]; b2=pars[2]; b3=pars[3];s1=pars[4];s2=pars[5]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))


oo=optim(.9 * dmle[-4], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))
return(c(-oo$value,oo$par))
}

proflik.path(dmle[2])


cc=matrix(NA,ncol=6,nrow=length(ii))
for(j in 1:length(ii))
cc[j,]=proflik.path(ii[j])
mle=min(abs(ii-dmle[2])) == abs(ii-dmle[2])
plot(ii,cc[,1] - cc[mle,1] + 2 ,type="n",ylim=c(-3,3),ylab="Log Likelihood",xlab="Male Slope Coefficient")
lines(ii,cc[,1] - cc[mle,1] + 2,lty=2,lwd=2,col=6)

abline(h=0,col=8)  
abline(v=dmle[2],col=8)

il=matrix(NA,nrow=length(y),ncol=length(ii))
#individual likelihoods
for(j in 1:length(y)){
pp<-""
pp<-paste(pp," 
+ dnorm(y[",j,"],u + b1 * x1[",j,"] + b2 * x2[",j,"]+ b3 * x1[",j,"] * x2[",j,"],exp(dmle[5]),log=T)",sep="")

pp<-paste("fn<-function(pars){ b1=pars[1]; u=pars[2]; b2=pars[3]; b3=pars[4]; return(-(",pp,"))}",sep="")
eval(parse(text=pp))

ic=NULL
for(k in 1:length(ii))
ic[k]=-fn(c(ii[k],cc[k,-1]))

il[j,]=ic - ic[mle] + 2/length(y)
lines(ii,il[j,] + rnorm(1,0,.025),col=colourl[[(1 + x2[j])]][rep(1:7,5)[j]],lty=c(1,3)[x2[j] + 1],lwd=2)


#text(ii,il[j,],myletters[j],col=colourl[[(1 + x2[j])]][rep(1:7,5)[j]])
}

lines(ii,apply(il[x2==0,],2,sum),col=4,lty=2,lwd=2)
lines(ii,apply(il[x2==1,],2,sum),col=2,lty=2,lwd=2)

title("When allowed to vary differently from men \n - women have no idea at all about the male slope")
dev.off()

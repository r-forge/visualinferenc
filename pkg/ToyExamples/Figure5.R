pdf(file="Figure5.pdf",width=10)

par(mfrow=c(1,2))


#prior for log odds
pr=function(lodds) lodds/((2 * pi^2) * sinh(lodds/2))


#numerically integrate likelihood 

f=function(p,lodds) 
(p^(nc) * (1 - p)^n_c) / (1 - p * ( 1 - exp(lodds)))^(ne + 1)

f(.5,3)

intll=function(lodds) lodds^-1 * (exp(lodds) -1 ) * exp(nce * lodds) * integrate(f,0,1,lodds=lodds)$value

intll(3)


#Meta_Analysis - 2 studies, no successes
study =c(rep(0, 2), rep(1, 2))
treatgroup =rep(1:2, 2)
success = c(0,0,0,0)
total = c(35,5,10,5)
failure = total - success

#Write out the likelihood
p=NULL
p[c(1,3)] = paste("lpc",1:2,sep="")
p[c(2,4)] = paste("lor + ",p[c(1,3)],sep="")
pp=""
for(i in 1:4){
pp=paste(pp," 
+ dbinom(success[",i,"],total[",i,"],","exp(",p[i],")/(1 + exp(",p[i],")),log=T)",sep="")
}


pp=paste("fn=function(pars){lor=pars[1];lpc1=pars[2];lpc2=pars[3];
return(-(",pp,"))}",sep="")
eval(parse(text=pp))

#List the log likelihood
fn

#Obtain mles (no very sensible here)
p0=c(0,0,0)
fn(p0)
oo= optim(p0,fn,method = "BFGS",hessian=T)
dmle=oo$par
dmle.se=sqrt(diag(solve(oo$hessian)))
rbind(dmle,dmle.se)


#Profile likelihood forlor
proflik=function(lor){

p=NULL
p[c(1,3)] = paste("lpc",1:2,sep="")
p[c(2,4)] = paste("lor + ",p[c(1,3)],sep="")
pp=""
for(i in 1:4){
pp=paste(pp," 
+ dbinom(success[",i,"],total[",i,"],","exp(",p[i],")/(1 + exp(",p[i],")),log=T)",sep="")
}


pp=paste("fn=function(pars){lor=lor;lpc1=pars[1];lpc2=pars[2];
return(-(",pp,"))}",sep="")
eval(parse(text=pp))

oo=optim(.9 * dmle[-1], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))
return(-oo$value)
}

#Calculate combined profile log likelihood
cc=NULL
ii=seq(-6,5.5,by=.2)

for(j in 1:length(ii))
cc[j]=proflik(ii[j])



#individual likelihoods
for(j in 1:2){
pp=""
for(i in (j*2 - 1):(j*2)){
pp=paste(pp," 
+ dbinom(success[",i,"],total[",i,"],","exp(",p[i],")/(1 + exp(",p[i],")),log=T)",sep="")
}


pp=paste("fn",j,"=function(pars){lor=pars[1];lpc1=pars[2];lpc2=pars[3];
return(-(",pp,"))}",sep="")
eval(parse(text=pp))
}

#List the study log likelihoods
fn1
fn2

#record the joint profile path
proflik.path=function(lor){

p=NULL
p[c(1,3)] = paste("lpc",1:2,sep="")
p[c(2,4)] = paste("lor + ",p[c(1,3)],sep="")
pp=""
for(i in 1:4){
pp=paste(pp," 
+ dbinom(success[",i,"],total[",i,"],","exp(",p[i],")/(1 + exp(",p[i],")),log=T)",sep="")
}


pp=paste("fn=function(pars){lor=lor;lpc1=pars[1];lpc2=pars[2];
return(-(",pp,"))}",sep="")
eval(parse(text=pp))

oo=optim(.9 * dmle[-1], fn,method="BFGS",hessian=T,control=list(reltol=.Machine$double.eps^.75))
return(c(-oo$value,oo$par))
}

#Plot the profile combined log likelihood
cc=matrix(NA,ncol=2+1,nrow=length(ii))
for(j in 1:length(ii))
cc[j,]=proflik.path(ii[j])
mle=min(abs(ii-dmle[1])) == abs(ii-dmle[1])
plot(ii,cc[,1] - cc[mle,1] + 1.36 ,type="n",ylim=c(-2,2),xlim=c(-4,8),ylab="Log Likelihood/Prior/Posterior",xlab="Log Odds Ratio")
lines(ii,cc[,1] - cc[mle,1] + 1.36,col=8,lty=2,lwd=2)

#Plot the study profile log likelihoos
ic=NULL
i=1
for(j in 1:length(ii))
ic[j]=-fn1(c(ii[j],cc[j,-1]))
i1=ic - ic[mle] + 1.36/2
lines(ii,i1,col=3,lty=1,lwd=2)

ic=NULL
i=1
for(j in 1:length(ii))
ic[j]=-fn2(c(ii[j],cc[j,-1]))
i2=ic - ic[mle] + 1.36/2
lines(ii,i2,col=3,lty=1,lwd=2)

#Plot the sum of study log likelihoods
lines(ii,i1 + i2,col=3,lty=2,lwd=2)

#Plot zero reference line
abline(h=0,col=8)

#try bayes using numerical integration

#recode outcomes
xc=success[1];nc=total[1];xt=success[2];nt=total[2]


nce=xt
n_ce=nt-xt
nc_e=xc
n_c_e=nc-xc
nc=nce + nc_e
n_c=n_ce + n_c_e
ne=nce + n_ce

#Calculate study1 integrated log likelihood
il=rep(NA,length(ii))
for(j in 1:length(ii))
il[j]=log(intll(ii[j] + .01))

xc=success[3];nc=total[3];xt=success[4];nt=total[4]


nce=xt
n_ce=nt-xt
nc_e=xc
n_c_e=nc-xc
nc=nce + nc_e
n_c=n_ce + n_c_e
ne=nce + n_ce

#Calculate study1 integrated log likelihood
i2=rep(NA,length(ii))
for(j in 1:length(ii))
i2[j]=log(intll(ii[j] + .01))

#Add the study log integrated likelihoods
ic= il + i2

#Find the integrated mle
milei=(1:length(ii))[ic == max(ic)]
mile=ii[milei]

#Calculate log prior and posterior
cc=NULL
for(j in 1:length(ii))
cc[j]=log(pr(ii[j] + .01))
#add logposterior
logpo=cc + il + i2

#Find the maximum posterior
mpei=(1:length(ii))[logpo == max(logpo)]
mpe=ii[mpei]

#Plot study log integrated likelihoods, combined log likelihood, prior and posterior
lines(ii,il - il[mpei] + 1.36/2,col=4,lty=1,lwd=2)

lines(ii,i2 - i2[mpei] + 1.36/2,col=4,lty=1,lwd=2)

lines(ii,il - il[mpei] + 1.36/2 + i2 - i2[mpei] + 1.36/2,col=4,lty=2,lwd=2)

lines(ii,cc - cc[mpei] + 1.36 ,col=2,lwd=2,lty=3)

lines(ii,cc - cc[mpei] + il - il[mpei] + 1.36/2 + i2 - i2[mpei] + 1.36/2,col=6,lty=4,lwd=2)


abline(v=c(0,mpe,5.5),col=8)
title("Individual Study and Pooled Log Likelihoods \n Using Numerical Integration")


#Now do Bayes by direct simulation of joint probability model
ii=seq(-6,8,by=.2)
n=28000000

#draw individual parameter from prior
pc1=rbeta(n,.5,.5)
pc2=rbeta(n,.5,.5)

#calculate parameter of focus log odds ratio
rlpc1=log(pc1/(1-pc1))
rlpc2=log(pc2/(1-pc2))
prior.or=rlpc2 - rlpc1

#obtain a density estimate for the prior distribution of this parameter
prior.dd=density(prior.or,from=min(ii),to=max(ii),n=length(ii))

#for study 1 obtain the posterior
#generate potential data from sampled prior
x1=rbinom(n,size=35, p=pc1)
x2=rbinom(n,size=5, p=pc2)

#condition on actual data to obtain posterior and get density estimate
is=(x1==success[1] & x2==success[2])
dd1=density(prior.or[is],from=min(ii),to=max(ii),n=length(ii))

#repeat for study 2
x1=rbinom(n,size=10, p=pc1)
x2=rbinom(n,size=5, p=pc2)
is=(x1==success[3] & x2==success[4])
dd2=density(prior.or[is],from=min(ii),to=max(ii),n=length(ii))

#subtract log prior to obtain the log likelihoods
ll1=(log(dd1$y) - log(prior.dd$y))
ll2=(log(dd2$y) - log(prior.dd$y))

#add to obtain the combined posterior
post=log(prior.dd$y) + ll1 + ll2

#obtain approximate max posterior
mpei=(1:length(ii))[post == max(post)]
mpe=ii[mpei]

#plot them
plot(ii,log(prior.dd$y) - log(prior.dd$y)[mpei] + 1.36,type="n",ylim=c(-2,2),xlim=c(-4,8),ylab="Log Likelihood/Prior/Posterior",xlab="Log Odds Ratio")
lines(ii,log(prior.dd$y) - log(prior.dd$y)[mpei] + 1.36,col=2,lty=3,lwd=2)

lines(ii,ll1 - ll1[mpei] + 1.36/2,col=4,lty=1,lwd=2)
lines(ii,ll2 - ll2[mpei] + 1.36/2,col=4,lty=1,lwd=2)
lines(ii,(ll1 + ll2) - (ll1 + ll2)[mpei] + 1.36,col=4,lty=2,lwd=2)
lines(ii,(log(prior.dd$y) + ll1 + ll2) - (log(prior.dd$y) + ll1 + ll2)[mpei] + 1.36,col=6,lty=4,lwd=2)


lines(ii,(log(dd1$y) - log(prior.dd$y)) +(log(dd2$y) - log(prior.dd$y)) - ((log(dd1$y) - log(prior.dd$y)) +(log(dd2$y) - log(prior.dd$y)))[mpei] + 1.36,col=4,lty=2,lwd=2) 
lines(ii,(log(dd1$y) - log(prior.dd$y)) +(log(dd2$y) ) - ((log(dd1$y) - log(prior.dd$y)) +(log(dd2$y) ))[mpei] + 1.36,col=6,lty=4,lwd=2) 
abline(v=c(0,mpe,5.5),col=8)

abline(h=1.36,col=3,lty=2,lwd=2)
abline(h=1.36/2,col=3,lty=1,lwd=2)

#Plot zero reference line
abline(h=0,col=8)

title("Individual Study and Pooled Log Likelihoods \n Using Direct Simulation")
dev.off()



pdf(file="Figure3.pdf",width=10)
par(mfrow=c(1,2))

#data from published meta-analysis
x=c("Allen",26,1,
"Alvi",43,3,
"Buller",116,24,
"Codegoni",31,8,
"Dellas",66,20,
"Fujita",47,8,
"Geisler",107,21,
"Gras",42,2,
"Han",19,1,
"Iwabuchi",95,6,
"King",41,7,
"Kobayashi",68,2,
"Krajinovic",12,2,
"Osborne",25,2,
"Shih",31,0,
"Sood",68,25,
"Sood",109,13,
"Tangir",31,0)

X=matrix(x,ncol=3,byrow=T)
pos=as.numeric(X[,3])
neg=as.numeric(X[,2])
num=pos + neg
pos/num

#Glm models
summary(glm(cbind(pos,neg) ~ 1,family=binomial))
summary(glm(cbind(pos,neg) ~ 1,family=quasibinomial))

#Write out the log likelihood
loglik=function(x,n,p)
x * log(p) + (n - x) * log(1 - p)

#set range for unknown p
p=(1:1499)/1500

#calculate combined log likelihood
sumll=loglik(sum(pos),sum(num),p) - max(loglik(sum(pos),sum(num),p))

#Calculate log prior and posterior
prior=dbeta(p,.5,.5,log=T) - max(dbeta(p,.5,.5,log=T)) + 2
post=sumll + prior - max(sumll + prior) + 2

#plot combined log likelihood
plot(p,sumll + 2,type="n",ylim=c(-4,8),xlim=c(0,1),
xlab="Proportion",ylab="Log Prior/Likelihood/Posterior")
lines(p,sumll + 2,lty=1,col=4,lwd=3)
title(main="Individual Log Likelihoods for Proportion \n Full View")

#plot individual log likelihoods
j=(1:length(sumll))[sumll==max(sumll)]
for(i in 1:length(num))
lines(p,loglik(pos[i],num[i],p)  - loglik(pos[i],num[i],p[j])  + 2/length(num),lty=1,col=8)

#Highlight within drop of 2 from maximum 
for(i in 1:length(num)){
il=loglik(pos[i],num[i],p)  - loglik(pos[i],num[i],p[j])  + 2/length(num)
maxil=max(il)
peak=(1:length(il))[ il - maxil > -2]
lines(p[peak],il[peak],col=4,lwd=2)
}


mpi=(1:length(p))[post == max(post)]
p[mpi]

#plot log prior and psoterior
lines(p,prior -prior[mpi] + 2,col=2,lwd=2,lty=3)
lines(p,post,col=6,lty=4,lwd=3)
abline(h=0,col=8)

#Replot over smaller interval for unknown p
plot(p,sumll + 2,type="n",ylim=c(-4,8),xlim=c(0,.25),
xlab="Proportion",ylab="Log Prior/Likelihood/Posterior")
lines(p,sumll + 2,lty=1,col=4,lwd=3)
title(main="Individual Log Likelihoods for Proportion \n Focused View")
j=(1:length(sumll))[sumll==max(sumll)]
for(i in 1:length(num))
lines(p,loglik(pos[i],num[i],p)  - loglik(pos[i],num[i],p[j])  + 2/length(num),lty=1,col=8)
for(i in 1:length(num)){
il=loglik(pos[i],num[i],p)  - loglik(pos[i],num[i],p[j])  + 2/length(num)
maxil=max(il)
peak=(1:length(il))[ il - maxil > -2]
lines(p[peak],il[peak],col=4,lwd=2)
}


mpi=(1:length(p))[post == max(post)]
p[mpi]

lines(p,prior -prior[mpi] + 2,col=2,lwd=2,lty=3)
lines(p,post,col=6,lty=4,lwd=3)
abline(h=0,col=8)
dev.off()


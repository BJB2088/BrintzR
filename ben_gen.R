head(birds)
n.it=t(birds)[8:57,]
X.const=rep(1,50)
Z.const=rep(1,300)
dim(birds)
data=n.it
R=50
T=6
#Tests
p=.5
lam=5
nmixbayes=function(params,data)
  R=dim(data)[1]
T=dim(data)[2]
N = as.vector(apply(n.it, 1, max))
Nrep=as.vector(rep(apply(n.it, 1, max),6))
vec=as.vector(data)
Bin.mat = matrix(dbinom(vec,Nrep,p),nrow=R,ncol=T)
N.pois.prior = dpois(N,lam)
lam.gam.prior = dgamma(lam,1000,rate=.001))


optimize(fn,interval=c(0,1),maximum = TRUE)





qplot(rgamma(10000,1000,rate=.0001),xlim=c(0,100),binwidth=.5)



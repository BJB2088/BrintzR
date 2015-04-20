birds=read.csv("birds.csv",header=TRUE)
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
R=dim(data)[1] # sites
T=dim(data)[2] # visits
nmixlik=function(p,N.it,lam){
#################### Create my likelihoods ####################
  Bin.mat = matrix(dbinom(vec,rep(N.it,T),p),nrow=R,ncol=T)
  N.pois.prior = dpois(N.it,lam)
  lam.gam.prior = dgamma(lam,shape=1.5,scale=100)
  lik=prod(as.vector(Bin.mat),N.pois.prior,lam.gam.prior)

}

# Gibbs
nmixbayes= function(p.in,lam.in){
  nparams = 52
  for (i in 1:nparams){
  # of params = 52 Ni's and p and lam
  N.in = as.vector(apply(n.it, 1, max)) # Initial values for N.i 
  Nrep=as.vector(rep(apply(n.it, 1, max),6)) # Put in binom matrix format
  vec=as.vector(data)
  upd=matrix(c(N.in,p.in,lam.in,rep(0,52*100)),nrow=52,ncol=101) 
}
}


fn=function(x) return(exp(-1/2*x^2))
1/(sum(fn(seq(-100,100))))

1/sqrt(2*pi)






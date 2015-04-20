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
nmixlik=function(params){
  #################### Create my likelihoods ####################
  Bin.mat = matrix(dbinom(vec,rep(params[1:50],T),params[51],log=TRUE),nrow=R,ncol=T)
  p.prior=dunif(params[51],0,1,log=TRUE)
  N.pois.prior = dpois(params[1:50],params[52],log=TRUE)
  lam.gam.prior = dgamma(params[52],shape=1.5,scale=100,log=TRUE)
  lik=sum(as.vector(Bin.mat),p.prior,N.pois.prior,lam.gam.prior)
  return(lik)
}
# Gibbs
nmixbayes= function(p.in,lam.in,n.reps){
  nparams = 52
  for(j in 1:n.reps){
  for (i in 1:nparams){
    # of params = 52 Ni's and p and lam
    N.in = as.vector(apply(n.it, 1, max)) # Initial values for N.i
    Nrep=as.vector(rep(apply(n.it, 1, max),6)) # Put in binom matrix format
    vec=as.vector(data)
    p.in=runif(1,0,1)
    lam.in=pgamma(1,shape=1,scale=100)
    upd=matrix(c(N.in,p.in,lam.in,rep(0,52*n.reps)),nrow=52,ncol=n.reps+1)
    
      upd[i,j]=
      params=upd[,j]
    nmixlik(upd[,j])
  }
}}



fn=function(x,y) return(exp(-1/2*x^2)*y)



1/(sum(fn(seq(-100,100))))
1/sqrt(2*pi)


hist(rgamma(10000,shape=2,scale=100))
rgamma(1,shape=1,scale=100)





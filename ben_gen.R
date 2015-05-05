x=rpois(50,3)
sum(x)
mat=matrix(rep(x,6),nrow=50,ncol=6)
n=matrix(rbinom(300,x,.5),nrow=50,ncol=6)


birds=read.csv("birds.csv",header=TRUE)
head(birds)
#n.it=t(birds)[8:57,]
n.it=n
start=c(as.vector(apply(n.it,1,max)),.2,2)
posterior <- function(param){
  N = rep(param[1:50],6)
  p = param[51]
  lam = param[52]
  singlelikelihoods = dbinom(as.vector(n.it), N, p, log = T)
  Nprior=dpois(param[1:50],param[52],log=TRUE)
  pprior=dunif(param[51],0,1,log=TRUE)
  lamprior=dgamma(param[52],shape=1.5,scale=100,log=TRUE)
  sumll = sum(singlelikelihoods+Nprior+pprior+lamprior)
  return(sumll)
}

proposalfunction <- function(param){
  mx=apply(n.it,1,max)
  dif=rpois(50,param[1:50])
  p1=apply(cbind(mx,dif),1,max)
  p2=runif(1,0,2*param[51])
  p3=rnorm(1,param[52],1)
  return(c(p1,p2,p3))
}
run <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,52))
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    proposal = proposalfunction(chain[i,])
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

out=run(start,500)
out[500,]

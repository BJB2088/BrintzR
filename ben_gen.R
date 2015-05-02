birds=read.csv("birds.csv",header=TRUE)
head(birds)
n.it=t(birds)[8:57,]
start=c(as.vector(apply(n.it,1,max)),.5,3)
posterior <- function(param){
  N = rep(param[1:50],6)
  p = param[51]
  lam = param[52]  
  singlelikelihoods = dbinom(as.vector(n.it), N, p, log = T)
  
  Nprior=dpois(params[1:50],params[52],log=TRUE)
  pprior=dunif(params[51],0,1,log=TRUE)
  lamprior=dgamma(params[52],shape=1.5,scale=100,log=TRUE)
  
  sumll = sum(singlelikelihoods+Nprior+pprior+lamprior)
  return(sumll)   
}
mx=apply(n.it,1,max)
proposalfunction <- function(param){
  dif=rpois(50,param[1:50])
  p1=mx+ifelse(mx<dif,rpois(50,dif-mx),0)
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






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
  N = rep(param[1:50],each=6)
  p = param[51]
  lam = param[52]
  singlelikelihoods = dbinom(as.vector(t(n.it)), N, p,log=T)
  Nprior=dpois(param[1:50],param[52],log=T)
  pprior=dunif(param[51],0,1,log=T)
  lamprior=dgamma(param[52],shape=1.5,scale=100,log=T)
  sumll = sum(singlelikelihoods,Nprior,pprior,lamprior)
  return(sumll)
}

proposalfunction <- function(param){
  mx=apply(n.it,1,max)
  dif=round(rnorm(50,param[1:50]))
  p1=apply(cbind(mx,dif),1,max)
  pval=rnorm(1,param[51],.2)
  
  p2= if(pval>1){1} else if(pval<0){0} else{pval}
  lval=rnorm(1,param[52],1)
  p3=ifelse(lval<0,0,lval)
  
  return(c(p1,p2,p3))
}

run <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,52))
  chain[1,] = startvalue
  count=0
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    
    if (runif(1) < probab){
      count=count+1
      print(count)
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
  
}

out=run(start,500)
out[500,]

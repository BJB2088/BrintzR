#Simulation
x=1:30
n=length(x)
y=7+3*x+rnorm(n,0,1)


colik.a <- function(a,b){
  pred = a + b*x
  singlelikelihoods = dnorm(y, mean = pred, sd=1,log=T)
  aprior = dnorm(a, 0, 1000,log=T)
  bprior = dnorm(b, 0, 1000,log=T)
  Sumll = sum(singlelikelihoods,aprior,bprior)
  return(prodll)   
}

# C=1/(integrate(Vectorize(coblikelihood),lower=-Inf,upper=Inf,b=2)$value)
proposalfunction <- function(prop){
  return(rnorm(2,mean = prop, sd= 0.1))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(colik.a(proposal[1],proposal[2])-colik.a(chain[i,1],chain[i,2]))
    
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

gibbsfn = function(a.start,b.start,its){ #a.start only used in first M-H
  lista=rep(0,its+1)
  listb=rep(0,its+1)
  lista[1]=a.start
  listb[1]=b.start
  #lista[1]=run_metropolis_MCMC(c(a.start,b.start),5000)[5000,1]
  #listb[1]=run_metropolis_MCMC(c(lista[1],b.start),5000)[5000,2]
  for (i in 1:its){
    lista[i+1]=run_metropolis_MCMC(c(lista[i],listb[i]),5000)[5000,1]
    listb[i+1]=run_metropolis_MCMC(c(lista[i+1],listb[i]),5000)[5000,2]
  }
  return(cbind(lista,listb))
}
res=gibbsfn(5,5,1000)


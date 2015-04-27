library(coda)
library(lattice)
library(rjags)
library(R2jags)
#sim test
x=rpois(50,3)
sum(x)
mat=matrix(rep(x,6),nrow=50,ncol=6)
n=matrix(rbinom(300,x,.5),nrow=50,ncol=6)

# bird test
n.it=t(birds)[8:57,]
Date=matrix(rep(2003:2008,50),nrow=50,ncol=6,byrow=TRUE)
# Royle Intercept Model
model.3 <- function(){
  for(i in 1:R){
    for (m in 1:V){
      y[i,m] ~ dbin(alpha.p, alpha.N[i])
    }
    alpha.N[i] ~ dpois(lambda)
  }
  
  alpha.p ~ dunif(0,1)

  lambda ~ dgamma(1.0E-6,1.0E-6)
  }

# --------------------
# 3. fit the model
# --------------------


# Initial values

#zst <- array(1,dim=c(Y,R))
#ust <- array(1,dim=c(S,Y,R))

n.it=n
inits <- function(){list(alpha.p = runif(1,0,1),
                         alpha.N = rep(8,nrow(n.it)),
                         lambda=3)}
params <- c("alpha.p","lambda","alpha.N")
test.data <- list(y = n.it, R = nrow(n.it), V = ncol(n.it))
system.time(out3 <- jags(test.data, inits, params, model.3, n.chains = 3,
                         n.thin = 20, n.iter = 20000, n.burnin = 1000))
print(out3, 2)
traceplot(out3)


sum(out3$BUGSoutput[10]$summary[1:50,1])

mean((as.vector(out3$BUGSoutput[10]$summary[1:50,1])-x)>1)


# Royle linear term model
Datevec=Date
model.3 <- function(){
  for(i in 1:R){
    for (m in 1:V){
      y[i,m] ~ dbin(alpha.p[i,m], alpha.N[i])
      alpha.p[i,m] <- p.int + p.cov.2004*ifelse(Date[m]==2004,1,0) + 
        p.cov.2005*ifelse(Date[m]==2005,1,0)+
        p.cov.2006*ifelse(Date[m]==2006,1,0) +
        p.cov.2007*ifelse(Date[m]==2007,1,0) +
        p.cov.2008*ifelse(Date[m]==2008,1,0)
    }
    alpha.N[i] ~ dpois(lambda)
  }
  
  p.int ~ dunif(0,1)
  p.cov.2004 ~ dnorm(0,1)
  p.cov.2005 ~ dnorm(0,1)
  p.cov.2006 ~ dnorm(0,1)
  p.cov.2007 ~ dnorm(0,1)
  p.cov.2008 ~ dnorm(0,1)
  
  lambda ~ dgamma(1.0E-6,1.0E-6)
}

inits <- function(){list(
                         alpha.N = rep(8,nrow(n.it)),
                         lambda=3,
                         p.int=.5, 
                         p.cov.2004=.1,
                         p.cov.2005=.1,
                         p.cov.2006=.1,
                         p.cov.2007=.1,
                         p.cov.2008=.1
                         )}
params <- c("p.int","p.cov.2004","p.cov.2005","p.cov.2006",
            "p.cov.2007","p.cov.2008","lambda","alpha.N")
test.data <- list(y = n.it, R = nrow(n.it), V = ncol(n.it),Date=Datevec[1,])
system.time(out3 <- jags(test.data, inits, params, model.3, n.chains = 3,
                         n.thin = 50, n.iter = 30000, n.burnin = 1000))
print(out3, 2)
traceplot(out3)


sum(out3$BUGSoutput[10]$summary[1:50,1])

# Dave's Model
model.3 <- function(){
  for(i in 1:R){
    y[i,1] ~ dbin(alpha.p, alpha.N[i,1])
    alpha.N[i,1] ~ dpois(lambda)
    for (m in 2:V){
      y[i,m] ~ dbin(alpha.p, alpha.N[i,m])
      alpha.N[i,m] ~ dpois(lambda)
    }
  }
  
  alpha.p ~ dunif(0,1)
  
  lambda ~ dgamma(1.0E-6,1.0E-6)
}

# --------------------
# 3. fit the model
# --------------------


# Initial values

#zst <- array(1,dim=c(Y,R))
#ust <- array(1,dim=c(S,Y,R))

n.it=n
inits <- function(){list(alpha.p = runif(1,0,1),
                         alpha.N = rep(8,nrow(n.it)),
                         lambda=3)}
params <- c("alpha.p","lambda","alpha.N")
test.data <- list(y = n.it, R = nrow(n.it), V = ncol(n.it))
system.time(out3 <- jags(test.data, inits, params, model.3, n.chains = 3,
                         n.thin = 20, n.iter = 20000, n.burnin = 1000))
print(out3, 2)

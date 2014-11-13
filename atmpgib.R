model.3 <- function(){
  for (i in 1:R){ # sites
    for (m in 1:Y){ # visits
      y[i,m] ~ dbin(p.val[i], N.est[i])
      N.est[i] ~ dpois(lambda[i])
      logit(lambda[i]) <- alpha.zero + alpha.one*location[i]
        }
}

lambda ~ dgamma(1000,1000)

}
zst <- apply(y.array, c(3,4), max)
ust <- apply(y.array, c(2,3,4), max)
inits <- function(){list(z=zst, u=ust)}

params <- c("alpha.zero","alpha.one","alpha.two","alpha.three")
test.data <- list(y = y.array, R = R, Y = Y, S = S, V = V, trt=trt,year=year)

# Call JAGS
system.time(out3 <- jags(test.data, inits, params, model.3, n.chains = 3,
                         n.thin = 50, n.iter = 50000, n.burnin = 7500))
# print(out3, 2)

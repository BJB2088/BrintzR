## INFO at http://wwwn.cdc.gov/nchs/nhanes/2007-2008/MCQ_E.htm#MCQ220 and 
## http://wwwn.cdc.gov/nchs/nhanes/2007-2008/DEMO_E.htm
# read in data files
MCQ <- read.xport("MCQ.xpt")
MCQweights <- read.xport("MCQ_weight.xpt")
head(MCQ)
head(MCQweights)
#combine weights with outomes # could also use WTMEC2YR
MCQ$weight <- MCQweights$WTINT2YR[ match(MCQ$SEQN,MCQweights$SEQN)]

# weights not even close to matching what XIA cited

# take out important data
Canc = MCQ %.% filter(MCQ220 == 1 | MCQ220 == 2) %.% 
  filter(!is.na(weight)) %.%
  select(SEQN,MCQ220,weight) %.%
  mutate(y=ifelse(MCQ220==2,1,0))
head(Canc,20)
table(Canc$MCQ220) #This matches the document! 

# Stab at MCMC
model <- function(){
for (i in 1:N {
y[i] ~ dbern(psi[i])
logit(psi[i]) <- nu + b1*weight[i]

}
# Priors
nu ~ dnorm(0,sqrt(3))
b1 ~ dnorm(0,sqrt(3))

# hyper paramter
p ~ dunif(0,1)
}
inits <- function(){list(nu=rnorm(1,0,sqrt(3)),b1=rnorm(1,0,sqrt(3)),w=.5)}
params <- c("nu","b1","p","w")
test.data <- list(y = Canc$y, weight = Canc$weight,N = length(Canc$y))

# Call JAGS
system.time(out <- jags(test.data, inits, params, model, n.chains = 3,
                         n.thin = 50, n.iter = 50000, n.burnin = 7500))
# print(out3, 2)

bird=t(birds)[8:57,]
colnames(bird) = c("2003","2004","2005","2006","2007","2008")

Tutorial on fitting the generalized N-mixture models in R

David Dail, April 2010

This tutorial demonstrates the closure test and population dynamics and abundance
estimation using the generalized N-mixture models as described by Dail and Madsen,
Biometrics, 2010. It largely follows the form of the tutorial provided by J.A. Royle
for fitting the original N-mixture model, given in the supplemental appendix of 
Kery, Royle, Schmid, Modeling avian abundance from replicated counts using binomial
mixture models, Ecological Applications, 2005.

This tutorial assumes the reader has some familiarity with R.

The "nmix.mig" takes a while to run, so the time required for finding MLEs will be 
reduced by working on a 64 bit computer (though this is not required).

A command intended to be entered as R code into the R console is depicted by the 
" ", with commentary between the commands describing the closure test or model 
selection using the generalized model.

First, load the R file gen.nmix.R into an active R session, using the command

  source("<pathway /gen.nmix.R")

where <pathway  is the complete path name of the directory containing the file
gen.nmix.R. There are 3 files contained in this archive, verified using the command:
  
    ls()

The data file "mallard.data" is the data file provided by Kery, Royle and Schmid 
(2005). It consists of Mallard point counts at 239 sites, each of which is sampled
up to 3 times, along with an assortment of possible covariates.

The other two files are R functions. The function "nmix.mig" fits the generalized 
model to data, finding the MLEs of the parameters. The function "ests" provides 
the estimated total abundance (with asymptotic standard error) during each primary
period in the study.

These functions require the logit and reverse-logit transform, which are built by 
the following commands:
  
expit = function(x) { return((exp(x))/(1+exp(x))) }
logit = function(x) { return(log(x/(1-x))) }


A couple preliminaries are required to fit the nmix.mig model to these data. First,
the covariate matrices for lambda and p be built (X and Z respectively), according
to the generalized linear models log(lambda) = XB and logit(p) = ZV, where B and V
are vectors of parameters.
t(birds) 

 data=bird
 #elev<-data[,"elev"]
 #length<-data[,"length"]
 #forest<-data[,"forest"]
 n.it<-data[,1:6]
 R<-nrow(n.it)
 T<-ncol(n.it)
 #DATE<-data[,c("date1","date2","date3")]
 #IVEL<-data[,c("ivel1","ivel2","ivel3")]
 #DATE2<-DATE^2


To fit the best model found in Kery, Royle, and Schmid (2005), we need the X and Z 
to have the form of X.lam and Z.p below:
  
   X.const = rep(1,R)
 #X.lam = cbind(X.const,elev,forest)

 p.date.lin = as.vector(t(DATE))
 Z.const = rep(1,R*T)
 Z.p = cbind(Z.const,p.date.lin)


The "sampling date" matrix giving the sampling dates associated with each survey 
must be constructed. Since the DATE covariate has been centered and standardized
in the "mallard.data" data set, the centering and standardizing values are needed
to obtained the raw sampling dates via back-transformation. The sampling design 
information provided by the authors is that the first sample occurred on date 1 and
the last sample occurred near date 90, with each sampling date a whole number. 
Therefore, the "sampling date" matrix is constructed using the following commands:  
  
    DATE.vec = c(DATE[,1],DATE[,2],DATE[,3])
  stds = sort(na.exclude(unique(DATE.vec)))
  diffs = numeric(length(stds)-1)
  for(i in 1:(length(stds)-1)){ diffs[i] = stds[i+1] - stds[i] }

  diffs.days = round(diffs/min(diffs))
  day.unique = numeric(length(diffs.days)+1)

  day.unique[1]=1
  for(i in 1:(length(diffs.days))){ day.unique[i+1] = day.unique[i] + diffs.days[i] }
  day.unique = c(day.unique,NA)
  DATE.mat = cbind(DATE.vec,seq(1:length(DATE.vec)))
  DATE.sort = DATE.mat[order(DATE.mat[,1]),]
  dup = duplicated(DATE.sort[,1])

  j=1
  for(i in 1:length(DATE.vec)){
     if(dup[i]==TRUE) DATE.sort[i,1]=DATE.sort[i-1,1]
     if(dup[i]!=TRUE) { 
         DATE.sort[i,1]= day.unique[j]
         j=j+1
       }
    }

  DATE.mat2 = DATE.sort[order(DATE.sort[,2]),]
  DATE.2 = matrix(DATE.mat2,nrow=239,ncol=3,byrow=FALSE)

DATE.2 is a matrix with the sampling date of every survey given for each site 
(1 is April 1 in this case).

Next, we need to construct a "date" matrix that provides the primary period number 
associated with each sample, under the alternative hypothesis of the population 
being open. For instance, given a length of primary period under the alternative 
hypothesis, if some site was not sampled during the first primary period, was 
sampled twice in the 2nd primary period and then once in the 3rd primary period, 
we would want the "date" record at that site to be [ 2  3  3 ].

The closure test depends on the length of the primary period chosen as the 
alternative to "the entire study duration". We will test the closure assumption 
here using the alternative duration of primary period length chosen as 30 days.

  DATE.3 = ceiling(DATE.2/30)


The function "nmix.mig" uses the date matrix, DATE.3, to construct the matrix of 
J.it values related to the number of samples obtained at each site during each
primary period, as well as the matrix of "Delta.it" values giving the length of
time (in primary periods) between sampling occasions at each site.

The optimization routine optim() will be used to find the MLEs of each model. The 
search method can be specified:
  
    Method = "BFGS"

More information related to optim() can be found by issuing the command

  help(optim)

The cut-off value K that will be used to approximate the infinite summations in 
both the original N-mixture model and the generalized model can also be specified:
  
    K.lim = 40

Note that while a cut-off value of 200 was used in the Dail and Madsen (2010) paper, 
using K=40 instead gives approximately the same answers and will require less 
computing time.

A list of the required inputs for "nmix.mig" is given below, along with a brief 
description:
  
  "vars" is a vector giving the parameter values. In general, this order is (lambda, 
                                                                             p, gamma, omega, over-dispersion), but "lambda" and/or "p" will be more than 
one value if more than an intercept is included in X and/or Z; gamma and omega 
are required only when the migration form requires them (i.e., when migration=
                                                           "constant"); and the over-dispersion parameter is only required when using the
negative binomial prior (i.e., prior="NB"). 
"n" is a matrix of the observed counts
"X" is a matrix of covariate values related to lambda       
"Z" is a matrix of covariate values related to p
"migration" is one of the following: "none" (for the original N-mixture model), 
"constant", "autoreg", "reshuf" for the different population dynamics models 
described by Dail and Madsen (2010)
"prior" is either "poisson" or "NB", and is the prior for initial abundance at 
each site
"Date" is a matrix of dates recording the primary period during which each 
observation was obtained 
"K" is the cutoff value for the infinite summation in obtaining the marginal 
likelihood.


The original N-mixture model with the negative binomial prior and no covariates
(ie, the null model) is fit with the following command:
  
    model.null = optim(c(0,0,0), nmix.mig, method=Method, hessian=TRUE, n=n.it,
                       X=X.const, Z=Z.const, migration="none", prior="NB", Date=DATE.3, K=K.lim)


The (0,0,0) in the above command gives the starting values in the search for the 
MLEs. The nmix.mig functions take a bit of computing time to find the MLEs, and it
is always a good idea to verify that the optim() procedure actually converged.

  model.null$conv

This will return "0" if the optim() is satisfied that it found the maximum of the 
likelihood (actually, optim() finds the minimum of the negative log likelihood).   

Furthermore, the stability of the MLEs can be checked by calculating the condition
number, which is the ratio of the largest eigenvalue of the hessian matrix to the 
smallest eigenvalue:
  
    ev.null = eigen(model.null$hessian)$values
  cn.null = max(ev.null)/min(ev.null)
  cn.null #the condition number

Values of the condition number near 0 or negative indicate a problem, possibly 
indicating fitting a model with too many parameters for the given data set.

The MLEs are obtained by the command:
  
    model.null$par

Because of a transformation that makes optim() run more smoothly, this gives the 
MLE for log(lambda), logit(p), and log(dispersion parameter) with the dispersion
parameter the overdispersion parameter of the negative binomial prior. The 
back-transformed MLEs for lambda and p can be obtained by the following commands:
  
    lambda.est = exp(model.null$par[1])
  p.est = expit(model.null$par[2])
  c(lambda.est, p.est)

The asymptotic standard errors of the parameter estimates are calculated with the 
Hessian matrix, and can be found with the following command:
  
    se = sqrt(diag(solve(model.null$hess)))
  se #the standard errors

As mentioned earlier, optim() is minimizing the negative log likelihood (nll). The 
minimum value of the nll can be retrieved by the command

  nll = model.null$val

With this, the AIC score for the null model can be calculated:
  
    aic = nll + 2*length(model.null$par)



In fitting the N-mixture model with covariates for lambda and p, the MLEs from the 
null model are used as initial values.

  model.closed = optim(c(model.null$par[1],0,0,model.null$par[2],0,model.null$par[3]),
                       nmix.mig, method=Method, hessian=TRUE, n=n.it, X=X.lam, Z=Z.p, migration="none",
                       prior="NB",Date=DATE.3,K=K.lim)

  model.closed$conv
  ev.closed = eigen(model.closed$hessian)$values
  cn.closed = max(ev.closed)/min(ev.closed)
  cn.closed #check condition number

The number of 0's inserted in the starting value vector is determined by the number
of covariates in the X and Z matrix for lambda and p (2 and 1, respectively, in 
this case)


Finally, the generalized model can be fit. We use the MLEs from the closed model as 
the starting point in the search for the generalized model MLEs: 

  model.open = optim(c(model.closed$par[1:5],-2,2,model.closed$par[6]), nmix.mig, 
method=Method, hessian=TRUE, n=n.it, X=X.lam, Z=Z.p, migration="constant", 
prior="NB", Date=DATE.3, K=K.lim)

  model.open$conv
  ev.open = eigen(model.open$hessian)$values
  cn.open = max(ev.open)/min(ev.open)
  cn.open #check condition number: this condition number is much larger than the others


The p-value from the closure test can then be obtained by using the following 
commands:

  t.stat= 2*(model.closed$val - model.open$val)
  obs.inf = -1*model.open$hess
  obs.inf.nn = obs.inf[1:5,1:5]
  obs.inf.np = obs.inf[1:5,6:7]
  obs.inf.pn = obs.inf[6:7,1:5]
  obs.inf.pp = obs.inf[6:7,6:7]
  I.tilda = obs.inf.pp - obs.inf.pn%*%solve(obs.inf.nn)%*%obs.inf.np
  prop = acos(I.tilda[1,2]/(sqrt(I.tilda[1,1]*I.tilda[2,2])))/(2*pi)
  prop.0 = 0.5 - prop
  prop.1 = 0.5
  prop.2 = prop
  p.value = prop.0*(0) + prop.1*(1-pchisq(t.stat,1)) + prop.2*(1-pchisq(t.stat,2))
  p.value


Even though this p-value is high enough (0.29) that there is no evidence against
the population closure hypothesis, we will get total abundance estimates (and 
asymptotic standard errors) for every primary period using both the open and the
closed models, as well as 95% confidence intervals for the population dynamics
parameters.

The first part is done using the "ests" function. The value for "T" that must be 
supplied indicates the duration of the study, in terms of the total number of 
primary periods. 

  ests.closed = ests(model.closed$par,model.closed$hess, migration="none", n=n.it,
X=X.lam, Z=Z.p, T=3, prior="NB")
  ests.closed

  ests.open = ests(model.open$par, model.open$hess, migration="constant", n=n.it, 
X=X.lam, Z=Z.p, T=3, prior="NB")
  ests.open

The first T (3 in this case) values are the estimated total abundance, and the
next T values are the asymptotic standard errors associated with each estimate, 
respectively.

Again, a back-transformation is required to retrieve the parameter estimates
from either the open or closed models with covariates. With covariates, each 
of the 239 sites has a different "lambda" value, and each (site, sampling 
occasion) combination can have a different detection probability (as in this
case). The mean of estimated lambda values and the mean of the estimated detection
probabilities are the next two values given as output in the ests function.

The last two values given as the output with the "ests" function are the 
back-transformed estimates of gamma and omega, the parameters controlling the 
population dynamics. 

  gamma = ests.open[9]
  gamma
  omega = ests.open[10]
  omega

Asymptotic standard errors of the dynamic parameter estimates (before 
back-transformation) can be obtained using the following commands:

  se = sqrt(diag(solve(model.open$hess)))
  se.gamma = se[6]
  se.omega = se[7]

Asymptotic 95% confidence intervals for gamma and omega can then be computed:

  gamma.ci = exp( c( model.open$par[6]- 1.96*se.gamma , model.open$par[6] +
1.96*se.gamma))
  gamma.ci
  omega.ci = expit( c( model.open$par[7] - 1.96*se.omega, model.open$par[7] +
1.96*se.omega))
  omega.ci


nmix.mig = function(vars,n,X,Z,migration="none",prior="poisson",Date=matrix(rep(1:length(n[1,]),length(n[,1])),nrow=length(n[,1]),ncol=length(n[1,]),byrow=TRUE), K=200) {
 # function returns the negative log likelihood of the generalized N-mixture model. 
 # David Dail, revised March 2, 2010
 #
 # migration is one these: "none", "constant", "autoreg", "trend1","reshuf".
 # migration="none" assigns r=0, w=1, so G~Pois(0) [ie, P(G=0)=1 ] and S~Binom(Nit-1,1) [ie, P(S=Nit-1) = 1 ];
 # migration="constant" has G~Pois(e^r0);
 # migration="autoreg" has G~Pois(e^r0*N.it-1).
 # migration="trend1" has r = (1-w)*lambda[1] ... needs no covariates for lambda
 # migration="reshuf" has r=lambda, w=0 ... needs no covariates for lambda
 #
 # X is the covariate matrix for log(lambda.i), including a column of 1's if an intercept is desired
 # Z is the covariate matrix for logit(p.it) [detection probability], which is a (R*T) X (no. of covariates) matrix.
 # The covariates for p.it need to be recorded across rows (sites): X11, X12, ..., X1T, X21, X22..., X2T, ..., XRT; 
 # need to include a column of 1's to include an intercept. 
 #
 # prior is either "poisson" or "NB".
 # 
 # w is survival probability, r is Poisson entering migration rate (assumed independent here)
 # Note that covariates for the migration parameters (r,w) are not accommodated in this program code.
 # 
 # vars is B0, B1,... Bk for each of: log(lambda.i), logit(p.it); then log(r0),logit(w) if migration != "none"; then log(dispersion).
 #
 # n is the matrix of number observed: sites are rows, and successive sampling occasions go across the columns,
 #  with NA (missing) entered at the end of row if less samples are obtained for that site (row) than the
 # maximum obtained for any site.
 #
 # Date is a matrix of the sampling dates (positive integers); it is best to have the first sampling occasion to be 1,
 # with NA listed (only when n.it is also missing) for missing data at the end of each row.
 #
 # K is the cut-off for the summations over N.it; the likelihood can be sensitive to this choice.
 # It is good practice to choose a value that is "high enough" so that increasing K does not change the resulting likelihood.
 # (Higher K's will take more computing power & time).

  n.it= n
  R = length(n[,1])
  nsite=R
#T.i is the number of sampling occasions per site
  T.i = numeric(R)
  for(i in 1:R){
    T.i[i] = sum( !is.na(n.it[i,]))
  }
 T=max(T.i)
#Delta.it is the time difference (no. of primary periods) between sampling occasions - could be 0
  Delta.it = matrix(0, nrow=R, ncol=(T-1))
  for(i in 1:R){
    for(t in 2:T){
      Delta.it[i,(t-1)] = Date[i,t] - Date[i,(t-1)]
    }
  }

  J.it= matrix(1,nrow=R,ncol=T)
  for(i in 1:R){
   for(t in 1:T){
     if(!is.finite(n.it[i,t])) J.it[i,t]=0
   }
  }

  T = max(T.i)
  X.lam.i=as.matrix(X)  # in null case this is a 1xR vector of ones
  X.theta.it=as.matrix(Z) # in null case this is a 1xRT vector of ones

  # obtain the parameter values for each site & sampling occasion

  lam.i = exp(X.lam.i%*%vars[1:length(X.lam.i[1,])])
  p.it.vec = expit(X.theta.it%*%vars[(1+length(X.lam.i[1,])):(length(X.lam.i[1,])+length(X.theta.it[1,]))])
  p.it = matrix(p.it.vec,nrow=R,ncol=T ,byrow=TRUE)
  
  if(prior!="poisson") {
    disp=exp(vars[length(vars)])
    disp2= max(disp, 0.0000001)
  }
  p.mat = matrix(p.it,nrow=R,ncol=T,byrow=TRUE)
  K.num = 0:K
  
  if(migration == "constant"){
    r0 = vars[1+(length(X.lam.i[1,])+length(X.theta.it[1,]))]
    w = expit(vars[2+(length(X.lam.i[1,])+length(X.theta.it[1,]))])
    r1 = 0
    r.it = exp(r0)*(0:K)^(r1)
  }
  if(migration == "autoreg") {
    r0 = vars[1+(length(X.lam.i[1,])+length(X.theta.it[1,]))]
    w = expit(vars[2+(length(X.lam.i[1,])+length(X.theta.it[1,]))])
    r1 = 1 # allows G.it ~ Pois(e^r0 * N.it-1)
    r.it = exp(r0)*(0:K)^(r1)
  }
  if(migration == "trend1") {
    r0 = 0
    w = expit(vars[1+(length(X.lam.i[1,])+length(X.theta.it[1,]))])
    r1 = 0 # allows G.it ~ Pois(e^r0 * N.it-1)
    r.it = (1-w)*lam.i[1] # requires no covariates for lambda
  }
  if(migration == "reshuf") {
    r0 = 0
    w = 0
    r1 = 0 # allows G.it ~ Pois(e^r0 * N.it-1)
    r.it = lam.i[1] # requires no covariates for lambda
  }

  if(migration == "none"){
    r0 = 0
    w = 1
    r.it = rep(0,(K+1))
  }

  # some preliminaries we'll need for each site, may as well build them now:
  
  Pois1 = rep(0:-K,(K+1)^2) + rep(0:K, each=(K+1)^2)
  r.it.rep = rep(rep(r.it,each=(K+1)),(K+1))
  Pois.mat = matrix(dpois(Pois1,r.it.rep),nrow=(K+1)^2,ncol=(K+1),byrow=TRUE)
  ####

  Bin1 = rep(0:(K),(K+1)^2)
  Bin2 = rep(rep(0:(K),each=(K+1)),(K+1))
  Bin.mat = matrix(dbinom(Bin1,Bin2,w),nrow=(K+1)^2, ncol=(K+1),byrow=TRUE)

  Pt.Agoal2 = matrix((Pois.mat*Bin.mat)%*%rep(1,(K+1)),nrow=(K+1),ncol=(K+1),byrow=FALSE)
  Pt.A2.array = array(0,dim=c(1+max(na.omit(Delta.it)),(K+1),(K+1)))
  Pt.A2 = diag((K+1))
  Pt.A2.array[1,,]=Pt.A2
  # using the Chapman-Kolmogorov equations:
  max.delta = max(max(na.omit(Delta.it)),max(na.omit(Date[,1]-1)))
  for(i in 2:(1+max.delta)){
    Pt.A2 = Pt.A2%*%Pt.Agoal2
    Pt.A2.array[i,,]=Pt.A2
  }

  probs.i = numeric(R)

  #Construct the likelihood, for each site (follows the supplemental appendix of Dail & Madsen.
  for (i in 1:R){

  #Make A
  
    T=T.i[i]
    p.t = p.it[i,]
    lam = lam.i[i]  
    delta.t = Delta.it[i,]
    n.t = n.it[i,]  
    J.t = J.it[i,]
    if(T 0){
      #make g_{1}(N_{iT})

      if(J.t[T]==0)break #this shouldn't happen...should only have observations in the data matrix!
      Bin.first = rep(c(n.t[(T-(J.t[T]-1)):T]), each=(K+1))
      Bin.index = rep(K.num,J.t[T])
      B=matrix(dbinom(Bin.first,Bin.index,p.t[T],log=TRUE),nrow=(K+1),ncol=(J.t[T]),byrow=FALSE)
      #B[,is.na(n[i,])]<-0
      g1.T = exp(B%*%rep(1,length(J.t[T])))
      g1.T.mat = matrix(rep(g1.T,(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE)
     # make g_{3}
     # if T=1 then what? skip this part, set A.goal = g1.T.
      if(T 1){
        if(delta.t[T-1] == 0) { 
          A.goal = matrix(g1.T,nrow=1,ncol=(K+1))
        }
        if(delta.t[T-1] != 0) {
          g3.T = Pt.A2.array[delta.t[T-1]+1,,]
          A.goal = matrix((g1.T.mat * g3.T)%*%rep(1,(K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
        }
      }
      if(T==1){
        A.goal = rep(1,(K+1))
      }

  #######
  # part B

      B.goal=A.goal

      hold.j = 1
      if(T 2) {
        for(reps in 1:(T-2)) {
          counter = T-1-reps
          Bin.first = rep(c(n.t[(T-hold.j-(J.t[counter+1]-1)):(T-hold.j)]), each=(K+1))
          Bin.index = rep(K.num,J.t[counter+1])
          B=matrix(dbinom(Bin.first,Bin.index,p.t[counter+1],log=FALSE),nrow=(K+1),ncol=(J.t[T]),byrow=FALSE)
          g1.T = apply(B,1,prod)
          g1.T.mat = matrix(rep(g1.T,(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE)
          hold.j = hold.j+1
      # make g_{3}
          if(delta.t[counter]==0){
            B.goal = B.goal*g1.T
          }      
          if(delta.t[counter]!=0){
            g3.T = Pt.A2.array[delta.t[counter+1]+1,,]
            #get g.star(N_{iT-1})

            pt.B3 = matrix(rep(B.goal,(K+1)),nrow=(K+1),ncol=(K+1),byrow=TRUE)
            B.goal = matrix((g1.T.mat * g3.T * pt.B3)%*%rep(1,(K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
          }
        }
      }


  ####
  # part C

      Bin.first = rep(c(n.t[(1:(T-hold.j))]), each=(K+1))
      Bin.index = rep(K.num,J.t[1])
      B=matrix(dbinom(Bin.first,Bin.index,p.t[1],log=FALSE),nrow=(K+1),ncol=(J.t[T]),byrow=FALSE)
      g1.T = apply(B,1,prod)

      # now we have g_{2}
      if(prior!="NB"){
        pt.C2.i = matrix(dpois((0:K),rep(lam.i[i],K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
      }
      if (prior=="NB"){
        pt.C2.i = matrix(dnbinom((0:K),mu=rep(lam.i[i],K+1),size=disp2),nrow=1,ncol=(K+1),byrow=TRUE)
      }

      pt.C2 = pt.C2.i/sum(pt.C2.i)
      pt.C3 = B.goal
      g1.T = as.vector(g1.T)
      pt.C2 = as.vector(pt.C2)
      pt.C3 = as.vector(pt.C3)
      # if the first sampling period for this site was the first overall period:
      if(Date[i,1] == 1){    
        probs.i[i] = sum((g1.T)*pt.C2*(pt.C3))
      } 

      # if the first sampling period for this site was NOT the first overall period:
      if(Date[i,1]!=1){
        g1.T.mat = matrix(rep(g1.T,(K+1)),nrow=K+1,ncol=K+1,byrow=TRUE)
        pt.B3 = matrix(rep(B.goal,(K+1)),nrow=(K+1),ncol=(K+1),byrow=TRUE)
        g4 = matrix((g1.T.mat* Pt.A2.array[Date[i,1],,]*pt.B3)%*%rep(1,(K+1)),nrow=1,ncol=(K+1),byrow=TRUE)
        g4 = as.vector(g4)
        probs.i[i] = sum(g4*pt.C2)
      }

    } # end of if(T 0)
    if((sum(J.t)==0)|(T==0)) probs.i[i]=1  

  } # end of site loop (i)

  q = -1*sum(log(probs.i))
  if(!is.finite(q)) return(1000000000)
  return(q)
}


ests = function(maxlikes,input.hess,migration="none", n, X, Z, T=length(n.it[1,]), prior="poisson"){
 #David Dail, revised March 2, 2010
 #
 #function to calculate the estimated values of N.t's (total abundances) for every sampling period
 #Will also provide the asymptotic SE's, using the Delta method
 #
 #maxlikes is the estimated MLE values of the parameters
 #input.hess is the hessian matrix of the likelihood evaluated at the MLEs (provided by optim())
 #migration="none" means r=0, w=1, closed popn assumption
 #migration="constant" means constant r (ie, log(r) = r0)
 #migration="autoreg" means r depends on N.it-1 (ie, log(r) = r0 + log(N.it-1)
 #if migration type is reshuffle or trend1, use migration="none" here.
 #
 #T is the number of primary sampling occasions
 #X and Z are the covariate matrices for lambda.i and p.it, respectively (constructed the same way
 # as in nmix.mig()
 #
 #prior="poisson" or "NB" (for negative binomial)

  X=as.matrix(X)
  Z=as.matrix(Z)
  betas = length(X[1,])
  phis = length(Z[1,])
  n.it=n
  R=length(n.it[,1])
  lam.i = exp(X%*%maxlikes[1:length(X[1,])])
  p.it.vec = expit(Z%*%maxlikes[(1+length(X[1,])):(length(X[1,])+length(Z[1,]))])
  p.it = matrix(p.it.vec,nrow=R,ncol=length(n.it[1,]),byrow=TRUE)
  ave.pit = mean(na.omit(p.it))

  if(prior=="NB") disp=exp(maxlikes[length(maxlikes)])
  
  if(migration=="none"){ #if there is not migration... a little simpler.
    inv.hess=abs(solve(input.hess)[1:betas,1:betas])
    lambda = lam.i
    par.ests = c(mean(lambda),ave.pit,0,1)
    N.hat = sum(lambda)
    partial.f.1 = N.hat
    partial.f = partial.f.1
    if(betas 1) {
      partial.f2.i = numeric(R)
      for (i in 1:R) {
        partial.f2.i[i] = X[i,2] * exp(X[i,]%*%maxlikes[1:betas])
      }
      partial.f2 = sum(partial.f2.i)
      partial.f = c(partial.f, partial.f2)
    }
    if(betas 2) {
      partial.f3.i = numeric(R)
      for (i in 1:R) {
        partial.f3.i[i] = X[i,3] * exp(X[i,]%*%maxlikes[1:betas])

      }
      partial.f3 = sum(partial.f3.i)
      partial.f = c(partial.f, partial.f3)
    }
    k.b=4
    while(k.b<=betas) {
      partial.f3.i = numeric(R)
      for (i in 1:R) {
        partial.f3.i[i] = X[i,k.b] * exp(X[i,]%*%maxlikes[1:betas])

      }
      partial.f3 = sum(partial.f3.i)
      partial.f = c(partial.f, partial.f3)
      k.b=k.b+1
    }
 


    var.Nhat = t(partial.f)%*%inv.hess%*%(partial.f)

    Nt = N.hat
    var.Nt = var.Nhat
    se.Nt = sqrt(abs(var.Nhat))
    N.ave = N.hat
    var.Nave = var.Nhat/T
    sd.Nave=sqrt(abs(var.Nave))

    trend=0
    var.trend = NA
    sd.trend=NA
    N.t = rep(Nt, T)
    seN.t = rep(se.Nt,T)
  }

  if(migration!="none") {  #if there is migration...
    Beta = maxlikes[1:betas]
    inv.hess=solve(input.hess)[c(1:betas,(betas+phis+1),(betas+phis+2)),c(1:betas,(betas+phis+1),(betas+phis+2))]
  
    w=expit(maxlikes[betas+phis+2])
    r=exp(maxlikes[betas+phis+1])
    lambda=lam.i
    par.ests = c(mean(lambda),ave.pit,r,w)
    ##hess and f.t need to be in this order: B0, B1, B2, ..., Bp, r, w

    N.t = numeric(T)
    seN.t = numeric(T)
    for(t in 1:T){   

      if( migration=="constant"){
        Nt.i = lambda*(w^(t-1)) + rep( r*(w^(t-1)-1)/(w-1), R)
        Nt = sum(Nt.i)
        St.i = lambda*((1-w^t)/(1-w)) - rep((r/(1-w))*(1-w^t)/(1-w),R) + rep(r*t/(1-w),R) 
        St = sum(St.i)
        N.ave = 1/t*St
        trend = rep(w,R) + rep(r,R)/lambda
        f.Nt.i = matrix(0,nrow=R,ncol=(betas+2))
        f.St.i = matrix(0,nrow=R,ncol=(betas+2))
        f.trend.i = matrix(0,nrow=R,ncol=(betas+2))

        for(i in 1:R){
          f.Nt.i[i,1]=lambda[i]*w^(t-1)
          if(betas 1) f.Nt.i[i,2] = X[i,2]* lambda[i]*w^(t-1)
          if(betas 2) f.Nt.i[i,3] = X[i,3]* lambda[i]*w^(t-1)
          k.b=4
          while(k.b<=betas) {
            f.Nt.i[i,k.b] = X[i,k.b]*lambda[i]*w^(t-1)
            k.b = k.b+1
          }
 
          f.Nt.i[i,betas+1] = r*(w^(t-1)-1)/(w-1)
          f.Nt.i[i,betas+2] = (lambda[i]*(w^(t-1)*(t-1))*(w-w^2)/w + 
                      r*w^(t-1)*(t-1)*(w-w^2)/(w*(w-1)) - r*(w^(t-1)-1)*(w-w^2)/(w-1)^2 )
      
          f.St.i[i,1]=lambda[i]*(1-w^t)/(1-w)
          if(betas 1) f.St.i[i,2] = X[i,2]*lambda[i]*(1-w^t)/(1-w)
          if(betas 2) f.St.i[i,3] = X[i,3]*lambda[i]*(1-w^t)/(1-w)
          k.b=4
          while(k.b<=betas) {
            f.St.i[i,k.b] = X[i,k.b]*lambda[i]*(1-w^t)/(1-w)
            k.b = k.b+1
          }
 
          f.St.i[i,betas+1] = -1*r*(1-w^t)/(1-w)^2 + r*t/(1-w)
          f.St.i[i,betas+2] = (-1*lambda[i]*w^t*t*(w-w^2)/(w*(1-w)) 
           - lambda[i]*(1-w^t)*(-w+w^2)/(1-w)^2 + lambda[i]*w^t*t*(w-w^2)/(w*(1-w)^2) 
           + 2*r*(1-w^t)*(-w+w^2)/(1-w)^3 - r*t*(-w+w^2)/(1-w)^2 )

          f.trend.i[i,1] = -r/(lambda[i]) 
          if(betas 1) f.trend.i[i,2] = -r*X[i,2]/(lambda[i])
          if(betas 2) f.trend.i[i,3] = -r*X[i,3]/(lambda[i])
          k.b=4
          while(k.b<=betas) {
            f.trend.i[i,k.b] = -r*X[i,k.b]/(lambda[i])
            k.b = k.b+1
          }
 

          f.trend.i[i,betas+1] = r/(lambda[i])
          f.trend.i[i,betas+2] = w-w^2
        } # end of site= 1,...,R loop
  
        f.Nt = colSums(f.Nt.i)
        f.St = colSums(f.St.i)
        #f.trend = colSums(f.trend.i)

        var.Nt = t(f.Nt) %*% inv.hess %*% (f.Nt)
        var.St = t(f.St) %*% inv.hess %*% (f.St)
        var.Nave = (1/t)^2*var.St
        se.Nt = sqrt(abs(var.Nt))
        se.Nave = sqrt(abs(var.Nave))
        trend.i = rep(w,R) + rep(r,R)/lambda
        var.trend.i = numeric(R)
      
        for(i in 1:R){
          var.trend.i[i] = t(f.trend.i[i,])%*%inv.hess%*%(f.trend.i[i,])
        }
        trend=sum(((trend.i-rep(1,R))^2)/var.trend.i)
      
        trend.i.se=sqrt(abs(var.trend.i))

      } #end of migration="constant"

    
      if( migration=="autoreg") {
        Nt.i = lambda*(w+r)^(t-1)
        Nt = sum(Nt.i)
        St.i = lambda*(((w+r)^t-1)/(w+r-1))
        St = sum(St.i)
        trend = w+r

        N.ave = 1/t*St
        f.Nt.i = matrix(0,nrow=R,ncol=(betas+2))
        f.St.i = matrix(0,nrow=R,ncol=(betas+2))
        for(i in 1:R){
          f.Nt.i[i,1]=lambda[i]*(w+r)^(t-1)
          if(betas 1) f.Nt.i[i,2] = X[i,2]* lambda[i]*(w+r)^(t-1)
          if(betas 2) f.Nt.i[i,3] = X[i,3]* lambda[i]*(w+r)^(t-1)
          k.b=4
          while(k.b<=betas) {
            f.Nt.i[i,k.b] = X[i,k.b]* lambda[i]*(w+r)^(t-1)
            k.b = k.b+1
          }
 
          f.Nt.i[i,betas+1] = lambda[i]*(w+r)^(t-1)*(t-1)*r/(w+r)
          f.Nt.i[i,betas+2] = lambda[i]*(w+r)^(t-1)*(t-1)*(w-w^2)/(w+r)

          f.St.i[i,1]=lambda[i]*(((w+r)^t )-1)/(w+r-1)
          if(betas 1) f.St.i[i,2] = X[i,2]*lambda[i]*(((w+r)^t )-1)/(w+r-1)
          if(betas 2) f.St.i[i,3] = X[i,3]*lambda[i]*(((w+r)^t )-1)/(w+r-1)
          k.b=4
          while(k.b<=betas) {
            f.St.i[i,k.b] = X[i,k.b]* lambda[i]*(((w+r)^t )-1)/(w+r-1)
            k.b = k.b+1
          }
 
          f.St.i[i,betas+1] =( lambda[i]*(w+r)^t*t*r/(w+r)/
                       ( w+r-1) - (lambda[i]*(w+r)^t - lambda[i])*r/(w+r-1)^2)

          f.St.i[i,betas+2] = (lambda[i]*(w+r)^t*t*(w-w^2)/(w+r)/(w+r-1) - 
            (lambda[i]*(w+r)^t - lambda[i])*(w-w^2)/ (w+r-1)^2 )
        } #end of site=1,...,R loop
        f.Nt = colSums(f.Nt.i)
        f.St = colSums(f.St.i)


        var.Nt = t(f.Nt) %*% inv.hess %*% (f.Nt)
        se.Nt = sqrt(abs(var.Nt))
        var.St = t(f.St) %*% inv.hess %*% (f.St)
        var.Nave = (1/t)^2*var.St
        se.Nave = sqrt(abs(var.Nave))
        f.trend = c(rep(0,betas), r, (w-w^2))
  
        var.trend=t(f.trend)%*%inv.hess%*%(f.trend)
        se.trend = sqrt(abs(var.trend))
      } #end of migration="autoreg"

      N.t[t] = Nt
      seN.t[t] = se.Nt
    } # end of loop through t= 1,...,T

  } # end of if migration!="none" 

  #could slightly modify to return N.ave, se.Nave, trend, se.trend also.

  ans.ret=c(N.t,seN.t,par.ests)
  return(ans.ret)
 
}

## Below is the Mallard data set that was originally published in the online Supplement
# of Kery, Royle, Schmid, Modeling avian abundance from replicated counts using 
# binomial mixture models, Ecological Applications, 2005.

"mallard.data" <-
structure(c(0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 
0, 1, 0, 1, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, NA, 4, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, NA, 
0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 2, 1, 0, 
0, 0, 1, 0, 0, 0, NA, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 
0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 
0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
0, 12, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
0, 0, 0, NA, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, NA, 0, 0, 0, 0, 1, 4, 0, 
0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, NA, 
0, 0, 1, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, NA, 
0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 1, 0, 0, 
0, 0, 0, NA, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, NA, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 1, 0, 0, 0, NA, 4, NA, NA, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, NA, 0, 0, 0, NA, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
1, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 1, 0, 
1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, NA, NA, 0, 0, NA, 0, 0, 
0, 3, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 3, 1, 0, 0, 0, 2, NA, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
0, 0, 0, 1, 0, 0, NA, 1, 0, 2, 0, NA, NA, NA, NA, 0, NA, NA, 
NA, 0, 1, 0, 0, 0, NA, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, -0.506, -0.934, 
-1.136, -0.819, 0.638, -1.329, -1.448, -0.321, -0.231, -1.097, 
-0.224, NA, 0.417, -1.354, -1.117, -0.278, -0.483, -1.07, -1.115, 
-0.862, -0.242, -0.462, -0.358, 0.417, -0.458, 1.866, -0.264, 
0.427, -0.42, 0.952, 0.904, -1.697, -1.067, -0.837, 0.568, 0.223, 
0.977, 0.952, 0.395, 0.338, -0.727, -0.819, 0.24, -0.139, -1.508, 
-0.101, -0.539, 0.897, -0.374, -0.308, -1.199, -0.984, -0.042, 
0.073, -1.097, 0.334, 0.662, -0.206, -1.568, 0.568, -1.097, -0.394, 
-1.212, 0.919, -0.019, -0.123, -0.735, -0.672, -0.476, 0.429, 
0.666, -0.417, -0.803, -0.172, -0.785, -0.792, -0.795, 0.647, 
-0.308, 0.298, 1.382, -0.394, -0.669, 0.476, -0.941, -0.601, 
0.735, 1.076, -0.212, -1.513, -0.958, 0.516, -0.692, -0.172, 
0.09, -0.731, 0.481, -1.428, -0.004, -0.443, -0.394, -0.366, 
-0.312, 0.152, 0.417, 0.558, -0.194, -0.431, -0.348, 0.984, -1.014, 
-1.148, -0.605, 0.88, -0.554, -0.355, -0.524, 0.355, -0.729, 
-1.393, -0.01, 0.006, -0.997, -1.108, -0.513, 1.973, 3.066, -0.413, 
-0.834, -1.377, -0.647, -0.984, -0.681, -1.148, -1.014, 0.104, 
1.168, -1.122, 0.557, -1.097, -0.652, 0.568, -0.469, -1.088, 
-0.086, NA, -0.735, -0.042, -1.249, -0.373, 1.013, -0.081, 0.342, 
-0.451, 1.225, -0.862, -1.196, -0.482, 0.806, 0.806, 0.735, -0.06, 
0.467, -1.485, 1.502, -0.264, 0.615, 3.066, -0.113, 2.119, 0.099, 
0.152, -0.35, 0.318, -0.51, 0.436, -0.086, -0.362, 0.502, 0.768, 
2.759, -0.451, 1.102, 0.066, 1.289, 0.028, 0.994, 1.734, 0.443, 
-0.139, 0.735, -0.448, 0.4, 0.659, 1.147, -0.194, 0.21, -0.745, 
-0.272, -0.344, -0.897, -0.819, 1.422, -0.643, 2.355, 1.379, 
-0.072, 1.193, 1.609, -0.291, 1.982, 0.323, -1.356, 0.107, -0.036, 
1.039, 0.493, 1.28, 0.201, 2.094, -0.203, 1.71, -0.506, 0.152, 
1.329, 0.251, -0.913, 0.092, 0.735, -1.228, 5.355, 0.066, 5.494, 
-0.378, -0.618, 0.264, -0.126, 3.713, 1.056, -0.506, -0.991, 
-1.339, -0.927, 0.88, -1.042, -1.562, -0.557, -0.231, -1.021, 
0.058, NA, 0.284, -1.014, -0.224, -0.182, -0.884, -1.258, -0.735, 
-1.115, -0.242, -0.615, -0.212, 0.682, -1.542, 1.883, 0.568, 
0.702, -0.264, 1.409, 0.913, -1.753, -1.543, -0.472, 0.152, -0.63, 
0.042, 1.123, -0.01, 0.586, -1.097, -0.738, 0.24, -0.139, -1.644, 
-0.228, -0.688, 0.491, -0.469, -0.308, -0.904, -0.589, -0.237, 
0.546, -1.216, 0.425, 0.735, 0.101, -1.52, 0.568, -1.051, -0.394, 
-1.274, 0.459, -0.362, 0.207, -0.228, -0.672, -0.61, 0.984, 0.666, 
-0.985, -0.903, -0.172, -0.889, -0.669, -1.014, 0.592, 0.305, 
0.298, 0.605, -0.212, -0.546, 0.476, -0.941, NA, 0.56, 0.863, 
0.334, NA, -1.143, 0.759, -0.615, 0.368, 0.214, -0.466, 0.199, 
-1.132, 0.256, -0.264, 0.061, -1.143, -0.047, 0.152, 0.735, 0.965, 
-0.293, -0.166, 0.035, 0.984, -1.055, -1.058, -0.264, 1.245, 
-0.025, -0.482, -0.472, 0.694, -0.729, -1.498, -0.01, 0.443, 
-0.505, -0.872, -0.666, 1.791, 3.066, -0.283, -0.655, -1.377, 
0.058, -0.539, -1.027, -1.148, -0.664, 0.104, 1.507, -0.334, 
0.314, -0.889, -0.652, 0.984, -0.182, -0.716, -0.026, NA, -0.482, 
-0.042, -1.417, -0.314, 1.145, -0.373, 0.215, -0.551, 1.225, 
-0.862, -0.935, -0.101, 0.806, 0.806, 0.851, -0.06, 0.073, -1.245, 
1.502, -0.264, 0.615, 3.267, -0.378, 1.609, 0.099, 0.298, -0.652, 
0.318, 0.086, 0.01, -0.264, -0.362, 0.385, 0.6, 1.379, -0.2, 
0.912, -0.019, 1.289, -0.902, 1.123, 1.234, 0.268, 0.152, 0.735, 
0.323, 1.268, 0.532, 0.507, -0.095, 0.385, -0.488, -0.855, -0.282, 
-1.208, -1.305, 1.497, -0.908, 2.355, 0.919, -0.072, 1.262, 1.344, 
0.025, 2.185, 0.323, -1.048, -0.117, 0.058, 1.039, 0.427, 0.908, 
0.201, 2.094, -0.203, 1.71, -0.882, 0.006, 1.161, 0.547, -0.577, 
0.092, -0.348, -0.615, 5.98, 0.323, 5.494, 0.615, -0.453, 0.264, 
0.152, 3.713, 2.262, -0.506, -1.162, -1.61, -1.197, 1.042, -0.899, 
-1.676, -0.636, -0.001, -0.832, -0.224, NA, 0.549, -1.159, -0.788, 
0.009, -0.819, -1.211, -0.482, -0.988, -0.399, -0.462, 0.735, 
0.639, -1.434, 0.923, -0.264, 0.482, -0.108, 1.237, 1.054, -1.585, 
-1.437, -0.639, 0.152, -0.203, NA, 1.18, 0.152, NA, -1.004, -0.981, 
0.682, 0.298, -1.474, -0.101, -0.737, 0.491, -0.66, -0.308, -1.157, 
-0.836, NA, 0.231, -1.037, 0.789, 0.516, -0.001, -1.615, 0.568, 
-1.189, -0.212, -1.491, NA, -0.819, 0.702, -0.228, -0.672, -0.52, 
0.984, NA, NA, -0.903, 0.152, -1.097, -0.546, -0.722, 0.207, 
-0.104, 0.006, 0.346, -0.394, -0.423, 0.638, -0.941, NA, 0.035, 
1.076, -0.03, NA, -0.264, 0.82, -0.615, -0.064, 0.214, -0.554, 
0.622, -1.28, 0.672, 0.271, -0.394, -1.273, NA, -0.2, -0.007, 
0.558, -0.144, -0.219, -0.098, 0.984, -1.014, -1.148, -0.491, 
1.245, 0.24, -0.608, -0.837, 1.033, -0.322, -1.217, 0.152, 0.37, 
-0.423, -0.99, -0.615, 1.609, NA, 0.108, -0.852, -1.424, 0.246, 
-0.539, -0.819, -1.148, -0.664, 0.582, 1.101, -0.161, 0.88, -0.847, 
-0.652, 2.649, -0.182, -0.902, 0.509, NA, NA, 0.152, -1.361, 
NA, 1.278, -0.442, 0.659, -0.551, NA, -1.115, -1.153, -0.101, 
0.806, 0.806, 1.201, 0.099, -0.084, -1.325, 1.076, -0.383, NA, 
2.965, -0.554, 1.609, 0.258, 0.103, -0.602, 0.735, -0.179, 0.436, 
-0.443, -0.362, 0.385, 1.105, 1.839, -0.099, NA, -0.448, 1.502, 
-0.282, 1.123, 1.234, 0.56, -0.431, 0.516, -0.448, 1.144, 0.279, 
1.289, 0.053, 0.735, -0.328, NA, -0.282, -0.936, -1.062, 1.557, 
NA, NA, NA, NA, 1.331, NA, NA, NA, 0.323, -1.528, 0.017, 0.105, 
1.039, NA, 0.747, 0.498, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, -1.761, -2.904, -1.69, 
-2.19, -1.833, -2.619, -2.69, -2.119, -2.047, -2.333, -1.69, 
-1.047, -1.833, -0.904, -2.547, -1.333, -1.833, -2.333, -2.404, 
-1.333, 0.024, -1.69, -1.976, -0.476, -2.404, -2.619, -2.761, 
-1.19, -0.619, -1.047, -2.19, -2.904, -1.761, -1.333, -2.119, 
-0.761, 1.096, 0.453, -2.119, 1.596, -2.19, -0.476, -1.833, -2.333, 
-1.833, 0.096, -1.833, -1.904, -1.833, -1.976, -2.119, -2.19, 
-0.547, -2.619, -1.69, -2.333, -2.404, -0.476, -2.69, -2.904, 
-2.404, -1.976, -2.404, 2.667, -2.261, -2.333, -2.404, -2.404, 
-2.547, -1.833, 1.596, 2.453, -1.904, -1.976, -2.047, -2.261, 
0.739, -1.833, -0.261, -2.19, -1.833, -2.904, -2.69, -2.404, 
-1.833, 0.453, -1.904, 0.524, -2.404, -2.904, -2.761, -2.333, 
-2.19, -2.404, -0.19, -2.261, -2.261, -1.904, -0.261, -0.976, 
-2.19, -2.404, 1.31, -0.761, -0.261, -0.404, -2.333, -2.547, 
-1.976, -0.833, -1.619, -1.19, -2.333, 1.239, -1.976, -2.333, 
-1.976, -2.333, -2.333, -2.333, -0.333, -2.333, -1.833, -2.333, 
-2.547, 2.167, 1.881, -2.547, -2.69, -1.976, -2.69, 1.096, -2.761, 
-2.761, -2.404, -1.976, -2.261, -2.547, 0.596, -2.69, -0.904, 
-1.119, -2.404, -2.333, -2.833, -2.333, 0.739, -1.047, -1.976, 
-0.761, 1.167, -0.19, -2.69, -1.904, 0.739, -2.547, -2.333, -1.833, 
-1.761, -1.761, -2.404, -2.404, -0.547, -2.761, 0.953, -2.333, 
-0.404, 0.096, -0.476, -0.547, -1.69, -2.333, -1.619, -2.119, 
-2.19, -2.261, -2.69, -1.619, -1.904, -2.904, 0.524, -2.333, 
1.596, -0.619, -1.047, -0.976, -1.833, -0.476, -2.904, -2.904, 
-1.833, -1.833, -0.047, -0.69, -1.619, -2.19, -2.904, -2.119, 
1.739, -1.833, -1.833, -2.19, -1.119, 0.667, 0.667, 0.667, 0.596, 
-1.976, 1.596, 1.667, 2.381, 0.667, -0.19, 0.453, 0.381, 0.096, 
1.453, -0.261, 0.81, 2.096, 0.453, 1.596, 0.167, 2.381, 2.31, 
2.096, 0.596, 1.81, 1.524, 2.381, 1.381, 2.239, 2.096, 1.167, 
0.667, 1.167, 0.667, 1.739, 1.81, 0.31, -1.047, -0.476, -0.69, 
0.167, 0.167, -1.19, -0.476, -0.547, -1.119, 0.453, 0.739, -0.333, 
1.096, -0.261, 0.31, -1.047, -0.976, -0.833, -0.261, 1.524, -0.261, 
-0.476, 0.596, -0.976, -0.547, -0.619, 0.453, 0.524, 0.453, -0.619, 
-0.976, -0.547, -0.261, -0.119, -0.19, 2.881, 1.381, -0.261, 
2.096, -1.261, 1.31, -0.761, -0.619, -0.833, 0.667, -1.047, -0.261, 
-0.119, -0.404, -0.619, -0.547, 0.81, -0.833, -1.19, -0.261, 
-1.333, 1.096, -0.761, -1.047, -1.261, -0.976, 0.096, 3.239, 
-0.619, -0.261, -0.976, -0.619, -1.19, -0.261, 2.667, 3.024, 
-0.976, -0.476, -0.547, -0.261, 1.739, -0.904, 1.596, 0.453, 
-0.761, -1.19, -1.261, -0.761, -0.833, 1.596, -0.833, 1.31, -0.976, 
-0.904, -1.404, -0.547, -0.547, -1.19, 1.167, -1.69, -0.833, 
-0.904, 1.096, -0.476, -1.19, -0.904, 3.239, 1.381, 0.381, 0.453, 
-0.833, -1.619, -0.547, 0.524, -1.119, -0.976, -0.261, 1.596, 
-0.19, -0.547, -0.547, -0.904, -0.976, -1.047, 2.667, -0.833, 
-0.261, -0.404, -0.547, 3.096, 3.024, -1.19, -0.547, 0.31, -0.833, 
1.596, -1.333, -1.047, -0.904, -0.619, -1.19, -0.547, 1.524, 
-0.19, -0.476, -0.119, -0.19, -0.619, -0.976, -0.476, 2.167, 
-0.476, -0.547, 0.096, 2.167, 1.31, -0.547, -0.619, 2.453, -0.261, 
-1.047, -0.476, -1.119, -1.119, -0.833, -0.833, 1.239, -1.761, 
2.381, -0.404, 1.596, 1.596, 0.596, -0.261, -0.619, 0.167, -0.261, 
-0.547, -0.619, -0.476, -1.69, -1.047, 0.096, -0.261, 2.167, 
-1.047, 3.024, 0.667, -0.261, 1.381, -0.619, 0.524, -1.047, -1.619, 
-0.404, -0.261, 1.667, 0.881, -1.047, -1.19, -0.547, 0.024, 2.31, 
0.096, -0.19, -0.904, -0.119, 1.667, 1.881, 1.596, 1.096, -0.904, 
2.596, 2.239, 2.881, 1.453, 0.596, 1.31, 0.524, 1.024, 2.453, 
0.881, 1.453, 3.239, 2.881, 2.167, 2.096, 3.453, 2.524, 3.096, 
1.453, 2.167, 1.739, 3.239, 1.81, 2.667, 2.667, 2.167, 1.596, 
3.167, 2.667, 2.81, 3.31, 1.381, 0.596, 1.453, 1.239, 1.381, 
1.381, 1.596, 1.453, 1.167, -0.261, 1.453, 2.596, 1.953, 1.239, 
1.31, 1.524, 0.524, 2.31, 0.667, 1.524, 3.239, 1.596, 1.239, 
1.667, 0.524, 0.524, 0.881, 2.596, 2.31, 2.096, 0.31, 1.453, 
-0.404, 1.31, 0.596, 0.596, NA, 2.167, 1.596, NA, 1.096, 1.953, 
0.596, 1.453, 0.596, 1.596, 0.596, 2.596, 2.167, 1.596, 0.524, 
1.524, NA, 0.596, -0.547, 0.524, 0.381, 1.667, 0.381, 0.096, 
0.453, -0.476, 1.667, NA, 0.167, 0.667, 1.024, -0.19, 0.596, 
0.667, NA, NA, 0.81, 1.31, 0.596, 1.524, 2.881, 0.167, 2.596, 
1.453, -0.119, -0.476, 0.31, 0.524, 0.453, NA, 0.596, 1.667, 
0.596, 1.096, 0.024, 1.381, 1.096, -0.261, 1.667, 0.381, 0.167, 
1.096, 1.739, 0.739, 0.524, -0.119, NA, 3.167, 1.596, 1.524, 
0.667, 0.524, 0.596, 1.667, 1.596, 0.096, 1.096, 1.953, 0.524, 
1.453, 0.596, 0.667, 0.667, 0.596, 2.739, 0.596, 0.667, 0.667, 
0.167, 3.524, NA, 0.881, -0.261, 0.81, 0.881, 2.096, 1.381, 0.524, 
1.096, 1.31, 0.524, 0.739, 1.667, 0.81, 0.524, 1.524, 2.024, 
0.881, 0.739, 0.667, NA, 1.096, 0.596, 1.31, 3.167, 2.167, 1.453, 
0.596, NA, 2.167, 0.667, 0.667, -0.261, -0.261, 1.596, 1.167, 
1.739, -0.119, 3.81, 0.596, NA, 3.381, 1.739, 0.596, 0.381, 1.31, 
0.596, 0.453, -0.261, 1.453, -0.619, -0.833, 1.881, 0.453, 2.596, 
-0.547, NA, 2.096, 0.453, 2.596, 0.096, 3.381, 0.096, -0.404, 
0.667, 1.953, 2.953, 1.739, -0.476, 0.81, 0.596, 1.739, NA, 1.596, 
0.667, 0.524, 1.31, NA, NA, NA, NA, 0.381, NA, NA, NA, 2.453, 
1.31, 1.453, 1.381, 2.024, NA, 1.596, 2.953, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
-1.173, -1.127, -0.198, -0.105, -1.034, -0.848, -0.91, -1.003, 
-0.058, -0.631, 0.066, -0.709, -1.111, -0.26, -0.647, 0.189, 
0.143, -1.235, -0.972, -1.219, 0.313, -1.003, 0.019, 0.484, -0.399, 
-0.863, -0.987, -0.554, 0.066, 0.004, -0.693, -0.925, -1.065, 
0.143, -0.383, -1.034, 1.01, 1.18, -0.136, 1.289, -0.6, 0.855, 
-0.089, -0.863, -0.91, 1.025, -0.383, -0.987, -1.173, 0.344, 
-0.445, 0.437, 1.041, 0.375, -0.894, -0.972, -0.662, 0.948, -0.77, 
-0.291, -1.018, -0.647, -1.25, 2.372, -0.817, -0.801, -0.693, 
-1.111, -1.266, -1.436, 1.18, 1.923, 0.22, -0.6, -1.096, -1.188, 
1.18, 0.607, 1.087, 0.174, 0.251, -0.461, -0.724, -0.848, -1.127, 
1.83, 0.731, 1.35, -0.6, -0.693, -1.065, -0.523, -1.127, -0.26, 
0.422, -0.507, 0.097, -0.941, 0.499, 0.375, -0.616, -1.096, 1.397, 
1.304, 0.886, 0.7, -0.832, -1.003, -0.801, 0.282, 0.174, -0.786, 
-1.08, 1.273, -0.585, -0.26, -0.77, -1.065, -1.204, -1.328, 0.344, 
-0.91, -1.157, -1.188, -1.034, 1.165, 0.592, -0.739, -1.142, 
-1.08, -1.065, 1.041, -1.065, -1.08, -1.142, -0.554, -0.786, 
-0.879, 0.809, -0.43, 0.112, -0.414, -0.987, -1.127, -1.096, 
-1.096, 1.087, 0.344, -1.096, -0.538, 1.196, 0.747, 0.05, -1.235, 
1.025, -0.352, -0.461, -1.08, -0.972, -0.972, -1.297, -0.972, 
0.22, -1.421, 1.289, -0.074, 1.134, 0.7, 0.561, 0.251, 0.004, 
-0.538, -0.987, 0.36, 0.22, -0.445, -0.925, -1.049, 0.066, -0.724, 
0.994, -1.018, 1.103, 0.731, 0.329, 0.329, 0.174, 0.979, -0.074, 
-1.08, -0.213, -1.188, 1.32, 1.087, -0.352, -1.127, -0.352, -0.167, 
1.242, 0.267, -0.987, -1.173, 0.375, 1.645, 1.242, 1.784, 1.35, 
0.004, 1.289, 1.614, 1.753, 1.567, 0.917, 1.01, 1.072, 1.025, 
1.923, 0.994, 1.041, 1.66, 1.66, 1.97, 1.815, 2.279, 2.434, 1.97, 
1.196, 1.505, 1.041, 1.35, 1.041, 1.35, 1.505, 1.35, 1.196, 1.196, 
1.66, 1.66, 1.815, 0.801, 0.115, -0.479, 0.315, -1.102, 0.741, 
0.115, -1.007, -0.913, 1.556, -1.626, 1.647, -0.399, 0.685, 0.801, 
0.741, -0.322, 0.801, 1.175, -0.245, -1.007, 0.503, -0.734, 0.381, 
-0.479, -1.301, -0.563, 0.252, 0.444, 0.115, 0.801, 0.185, 0.381, 
0.444, -0.563, -0.647, 0.252, 0.115, -1.102, -0.168, 0.857, -1.102, 
-1.406, -0.734, 1.944, 1.175, 0.626, -0.479, 0.741, 0.503, 1.175, 
0.626, -0.322, -1.007, -0.025, -1.514, -0.734, 0.503, 0.741, 
-0.563, 0.857, -0.095, 2.255, -3.336, 0.115, 0.252, 1.175, -0.245, 
0.965, -0.563, 0.115, -0.647, 0.566, 0.315, 1.86, 1.273, -0.734, 
0.252, 0.503, -0.734, -0.322, -0.095, 1.273, 0.315, -0.095, 2.063, 
0.045, -0.647, -0.095, 0.444, 0.857, -0.095, 1.51, 0.315, -0.168, 
-1.406, 0.801, 0.626, 0.444, -0.025, -1.514, -0.322, -0.399, 
0.566, 0.381, -0.479, 0.626, 0.381, -1.199, -1.199, 1.224, 0.965, 
1.556, -0.095, -1.406, -0.245, 0.444, -0.479, -0.479, 1.818, 
0.315, -0.734, 1.273, 1.416, 0.503, -0.095, -1.738, 1.07, 0.965, 
0.741, 0.801, 0.626, -0.563, 0.965, 0.045, 0.741, -0.479, 0.315, 
-1.102, 1.224, 0.566, -0.563, 0.741, -0.168, -0.025, 1.322, 1.175, 
-0.322, 0.185, 0.045, -0.399, 0.045, -0.245, 0.566, -0.913, 1.175, 
1.07, -0.245, -0.025, -0.025, 0.045, 0.381, -1.007, 1.371, -0.647, 
-0.025, -0.399, -1.857, -1.406, -0.734, 0.381, 0.685, 0.566, 
-1.199, -0.399, -0.647, -0.025, 0.115, 0.045, 0.185, -0.913, 
0.566, -0.245, -1.301, -0.647, -0.168, -0.322, -1.199, 0.045, 
0.685, -0.734, -1.301, -0.168, -0.245, -0.647, 0.626, 0.045, 
2.14, 0.381, -0.168, 1.465, 1.322, -0.822, 1.017, -0.647, -0.913, 
-0.822, -0.563, -0.399, -0.245, -0.479, 1.122, 1.902, 0.965, 
0.801, -0.245, 0.252, -0.168, 0.626, -1.102, -0.647, -0.479, 
2.217, -3.158, 0.185, 0.626, 0.185, -0.025, -1.199, 0.503, -4.406, 
-1.301, -4.945, -0.399, 0.252, 0.185, -0.563, -2.109, -1.857, 
-1.156, -0.501, -0.101, 0.008, -1.193, 0.917, -1.083, -0.792, 
0.553, 0.808, 1.79, -0.72, -0.647, 0.844, -0.32, -0.101, 1.135, 
0.372, -0.138, -0.974, 1.572, 0.081, -0.611, 0.699, -1.047, -1.156, 
-0.283, 0.772, 1.135, 0.299, 0.19, -1.083, -1.156, 0.626, 0.117, 
0.517, -0.72, 0.99, 1.717, -1.265, -1.083, 0.844, 2.299, -0.065, 
0.663, -0.647, -0.683, 0.917, -1.193, 0.626, 0.844, 0.553, -0.829, 
-0.21, -0.683, 1.754, 1.608, 0.59, -0.974, 0.372, -1.083, 1.063, 
-0.756, -1.265, -1.265, -0.174, -0.647, 1.099, -0.174, -1.265, 
-1.083, -1.265, 1.135, -0.974, 1.463, 1.208, 1.354, 2.19, 0.808, 
1.063, 0.335, 0.335, -0.029, -0.501, -0.429, -1.229, 0.299, -0.174, 
-0.32, -1.011, 1.899, 0.517, 0.153, 2.117, 0.481, -0.865, 1.39, 
-0.829, 0.262, -0.32, -0.792, -1.193, -0.72, 0.699, 1.426, 0.372, 
0.954, -0.065, 0.699, 1.245, 0.444, -1.083, 0.044, 0.917, 1.354, 
-0.21, -1.011, -1.011, -0.792, -1.193, -0.683, 0.99, -0.938, 
0.081, 1.317, 0.226, -0.574, -0.392, -0.902, -1.265, 1.281, 0.372, 
-0.611, -0.538, -1.156, -0.538, 0.044, -0.574, 0.008, 1.899, 
0.59, -0.501, -0.865, -0.21, 0.735, 0.699, -0.101, -0.101, -1.156, 
0.844, 0.044, 1.863, 1.463, -0.429, -1.265, 1.317, 1.245, -1.047, 
1.172, 1.172, 0.262, 1.354, 1.608, 0.444, -0.32, 0.772, -0.611, 
-0.065, -0.501, 0.153, 1.245, 0.262, 2.299, 0.626, 2.008, 0.553, 
0.699, -0.938, 0.59, 0.008, 0.808, -0.611, -0.829, -0.029, 1.208, 
1.026, 0.481, -0.065, 2.008, -0.065, 0.081, -0.974, 0.044, 1.135, 
0.772, -1.265, 0.226, 1.026, -1.265, 0.008, -0.756, -1.229, -0.283, 
-1.265, -1.047, -1.265, -1.229, 0.844, -0.938, -0.756, -1.193, 
-0.174, 2.117, 1.717, 1.681, 1.717, -1.265, 1.245, 0.663, -1.265, 
-1.265, -1.265, -1.265, -1.265, -1.265, -1.265, -1.265, -1.265, 
-1.229, -1.265, -1.265, -0.72, -1.265, -1.265, -1.193, -1.265, 
-1.265, -1.265, -1.265), .Dim = as.integer(c(239, 12)), .Dimnames = list(
    NULL, c("count1", "count2", "count3", "ivel1", "ivel2", "ivel3", 
    "date1", "date2", "date3", "elev", "length", "forest")))




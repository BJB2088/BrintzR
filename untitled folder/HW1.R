npts=65
japanese<-matrix(scan("http://www.lancaster.ac.uk/staff/diggle/pointpatternbook/datasets/japanese.txt"),npts,2,T)
plot(japanese)

x=seq(0,1,.2)


plots=data.frame(c1=rep(0,65),c2=rep(0,65))
for (k in 1:npts){
  for (i in 1:5){
    for (j in 1:5){
      for (l in 1:2){
      if(japanese[k,l] < x[2]) {plots[k,l] = 1
      } else if (japanese[k,l] < x[3]) {plots[k,l] = 2
      } else if (japanese[k,l] < x[4]) {plots[k,l] = 3
      } else if (japanese[k,l] < x[5]) {plots[k,l] = 4
      } else {plots[k,l] = 5}
      
}}}}
      

results=plots %.% group_by(c1,c2) %.% summarise(obs=n())
results=rbind(results,c(2,2,0))      
results$expected=rep(2.6,25)
results$chisqcontr=((results$obs-results$expected)^2)/results$expected
pchisq(sum(results$chisqcontr),24,lower.tail=FALSE)
# p-value is .8519242, no evidence to reject this is homogenous spatial poisson process


?ks.test
jap=as.data.frame(japanese)
colnames(jap)=c("x","y")
jap$numb=1:65
arb=sample(1:65,1)
plot(jap,type="n")
plot(jap$x,jap$y,type="n")
text(jap$x,jap$y)




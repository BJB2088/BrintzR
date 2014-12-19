

p=.2#p hidden proportion
theta_acc=.2
theta_hid=.6
yest=p*theta_hid + (1-p)*theta_acc
# read in data files
MCQ <- read.xport("MCQ.xpt")
MCQweights <- read.xport("MCQ_weight.xpt")
head(MCQ)
head(MCQweights)
#combine weights with outomes
MCQ$weight <- MCQweights$WTINT2YR[ match(MCQ$SEQN,MCQweights$SEQN)]

# take out important data
Canc = MCQ %.% filter(MCQ220 == 1 | MCQ220 == 2) %.% 
  filter(!is.na(weight)) %.%
  select(SEQN,MCQ220,weight)
table(Canc$MCQ220) #This matches the document! 


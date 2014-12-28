Ebola <- read.csv("Ebola.csv",header=TRUE)
unique(Ebola$Indicator)

# Which indiator should we look at? Probably number of confirmed Ebol cases
# in the last 21 days.Or last 7 days? Or could look at confirmed probable and 
# suspect...
# Start with confirmed in last 7 days

Ebola %.% filter(Indicator == "Number of confirmed Ebola cases in the last 7 days")

# Not a ton of data...Only have two dates for Guinea, Liberia, and Sierra Leone

Eb21=Ebola %.% filter(Indicator == "Number of confirmed Ebola cases in the last 21 days")

# That's better...I guess 
table(Eb21$Country) # Guinea, Liberia, and Sierra Leone each have 16 data points, others less
# Pulling out those countries 
Eb21 %.% filter(Country=="Guinea"|Country=="Liberia"|Country=="Sierra Leone")

# Unfortunately dates close enough that there should be overlap...
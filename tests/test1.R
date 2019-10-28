library(HCV)

data(sim1) # reads in d1, d2, d3, d4, d5 (i.e. Dataset_1.csv, Dataset_2.csv, ... etc)

# functions have now been substantially modified

# add (positive) jitter to break ties
# d1[d1>0] <- d1[d1>0] + rexp(sum(d1>0),1000)
# d2[d2>0] <- d2[d2>0] + rexp(sum(d2>0),1000)
#
# ml1 <- optimPL(d1)
# ml2 <- optimPL(d2)
#
# dList <- list(d1,d2)
# mlAll <- optimPLmulti(dList)
#
# p <- simulateProcess(param=exp(ml2$par),startcell=which(d2==max(d2)),nmax=sum(d2>0),nrow=nrow(d2),ncol=ncol(d2))
#
# btpar <- bootCI(exp(ml1$par),d1,5)

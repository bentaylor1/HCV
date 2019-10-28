library(HCV)

data(sim1) # reads in d1, d2, d3, d4, d5

# functions have now been substantially modified

# add (positive) jitter to break ties
# d1[d1>0] <- d1[d1>0] + rexp(sum(d1>0),1000)
# d2[d2>0] <- d2[d2>0] + rexp(sum(d2>0),1000)
# d3[d3>0] <- d3[d3>0] + rexp(sum(d3>0),1000)
# d4[d4>0] <- d4[d4>0] + rexp(sum(d4>0),1000)
# d5[d5>0] <- d5[d5>0] + rexp(sum(d5>0),1000)
#
# ml1 <- optimPL(d1)
# ml2 <- optimPL(d2)
# ml3 <- optimPL(d3)
# ml4 <- optimPL(d4)
# ml5 <- optimPL(d5)
#
# print(exp(ml1$par))
# print(exp(ml2$par))
# print(exp(ml3$par))
# print(exp(ml4$par))
# print(exp(ml5$par))

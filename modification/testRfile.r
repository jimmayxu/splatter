load(file="../data/Buettner.RData")
load(file="../data/Kolodziejczyk.RData")
load(file="../data/Pollen.RData")
load(file="../data/Usoskin.RData")


setwd("../../splatter_simulation/R")
sapply(list.files(), source)
library(checkmate)
library(MCMCpack)
library(scater)
library(BiocGenerics)

Test <- Kolodziejczyk
counts <- round(10^Test$in_X-1)
params <- splatEstimate(counts)

#getpara <- getParams(params,c("mean.rate","mean.shape","out.prob","out.facLoc","out.facScale"))

probs <- as.vector(rdirichlet(1,c(7,7,7)))
params <- setParams(params, group.prob = probs, de.prob = 0.2, dropout.type = "experiment")
params <- setParams(params, group.prob = probs, de.prob = 0.2, dropout.type = "none")

### Toy params ###
params <- newSplatParams()
probs <- as.vector(rdirichlet(1,c(7,7,7)))
params <- setParams(params, group.prob = probs, de.prob = 0.2, dropout.type = "experiment")


sim.groups <- splatSimulate(params, nGenes = 5000, method = "groups", batchCells = 200)
sim.groups <- splatSimulate(params, method = "groups")

sim.groups <- normalise(sim.groups)
plotPCA(sim.groups, colour_by = "Group")


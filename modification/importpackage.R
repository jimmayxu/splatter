setwd("../R")
sapply(list.files(), source)
setwd("../modification/")

setwd("../../SIMLR(John)/FZINB/")
source("packageimport.r")
setwd("../../splatter/modification/")
library(checkmate)
library(MCMCpack)
library(scater)
library(BiocGenerics)


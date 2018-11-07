setwd("../R")
sapply(list.files(), source)

setwd("../../FZINB_kernel/SIMLR/")
source("packageimport.r")

#dyn.load("../src/projsplx_R.so")
setwd("../../splatter/modification/")
library(checkmate)
library(MCMCpack)
library(scater)
library(BiocGenerics)


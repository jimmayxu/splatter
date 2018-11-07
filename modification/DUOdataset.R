
source("importpackage.R")
#summary <- readRDS("../../DuoClustering2018dataset/clustering_summary.rds")


#### True number of clusters obtained from supplementry figure 1 of Duoclustering
true_lab <- list(
  "Koh" = 9,
  "KohTCC" = 9,
  "Kumar" = 3,
  "KumarTCC" = 3,
  "SimKumar4easy" = 4,
  "SimKumar4hard" = 4,
  "Zhengmix8eq" = 8
)


Koh <-readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_Koh.rds")
KohTCC <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_KohTCC.rds")
Kumar <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_Kumar.rds")
KumarTCC <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_KumarTCC.rds")
SimKumar4easy <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_SimKumar4easy.rds")
SimKumar4hard <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_SimKumar4hard.rds")
Zhengmix8eq <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_Zhengmix8eq.rds")
SCE <- list(Koh,KohTCC,Kumar,KumarTCC,SimKumar4easy,SimKumar4hard,Zhengmix8eq)


X_counts <- assays(SimKumar4easy)$counts
sum(X_counts==0)/prod(dim(X_counts))
params <- splatEstimate(X_counts)
### Check for the library size
#hist(colSums(X_counts))
#sqrt(var(colSums(X_counts)))/mean(colSums(X_counts))


n_clusters<- 4

probs <- as.vector(rdirichlet(1,rep(7,n_clusters)))
params <- setParams(params, group.prob = probs, de.prob = 0.2, dropout.type = "experiment")
sim.groups <- splatSimulate(params, method = "groups", nGenes = 5000, batchCells = 200)

sim.groups <- normalise(sim.groups,return_log = FALSE)
assays(sim.groups)$logcounts <- log(assays(sim.groups)$normcounts+1,10)

plotPCA(sim.groups, colour_by = "Group")


sim.groupsKumarTCC <- sim.groups
sim.groupsZhengmix8eq <- sim.groups
sim.groupsSimKumar4easy <- sim.groups

SCE.data <- list("KumarTCC" = sim.groupsKumarTCC,
                 "Zhengmix8eq" = sim.groupsZhengmix8eq,
                 "SimKumar4easy" = sim.groupsSimKumar4easy)

sim.data <- list("KumarTCC" = assays(SCE.data$KumarTCC),
                 "Zhengmix8eq" = assays(SCE.data$Zhengmix8eq),
                 "SimKumar4easy" = assays(SCE.data$SimKumar4easy))
sim.lab <- list("KumarTCC" = as.numeric(factor(colData(SCE.data$KumarTCC)$Group)),
                 "Zhengmix8eq" = as.numeric(factor(colData(SCE.data$Zhengmix8eq)$Group)),
                 "SimKumar4easy" = as.numeric(factor(colData(SCE.data$SimKumar4easy)$Group)))


dropout.mid <- getParam(metadata(SCE.data$Zhengmix8eq)$Params,"dropout.mid")
dropout.shape <- getParam(metadata(SCE.data$Zhengmix8eq)$Params, "dropout.shape")

params <- metadata(SCE.data$SimKumar4easy)$Params
params <- setParams(params, update = list(dropout.mid = dropout.mid, dropout.shape = dropout.shape))
sim.groups <- splatSimulate(params, method = "groups", nGenes = 5000, batchCells = 200)
sim.groups <- normalise(sim.groups,return_log = FALSE)
assays(sim.groups)$logcounts <- log(assays(sim.groups)$normcounts+1,10)
SCE.data$SimKumar4easy <- sim.groups

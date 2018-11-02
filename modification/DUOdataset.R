
source("importpackage.R")
summary <- readRDS("../../DuoClustering2018dataset/clustering_summary.rds")
Koh <-readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_Koh.rds")
KohTCC <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_KohTCC.rds")
Kumar <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_Kumar.rds")
KumarTCC <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_KumarTCC.rds")
Kumar <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10/sce_filteredExpr10_SimKumar4easy.rds")
Kumar <- readRDS("../../DuoClustering2018dataset/sce_filteredExpr10")
SCE <- list()

params <- splatEstimate(assays(KumarTCC)$normcounts-1)

n_clusters<- 3

probs <- as.vector(rdirichlet(1,rep(7,n_clusters)))
params <- setParams(params, group.prob = probs, de.prob = 0.2, dropout.type = "experiment")
sim.groups <- splatSimulate(params, method = "groups")

sim.groups <- normalise(sim.groups,return_log = FALSE)
assays(sim.groups)$logcounts <- log(assays(sim.groups)$normcounts,10)

### Check for the library size
hist(colSums(counts(sim.groups)))
var(colSums(counts(sim.groups)))/mean(colSums(counts(sim.groups)))

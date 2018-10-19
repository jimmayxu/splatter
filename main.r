setwd("~/Documents/Git/Github/SIMLR(John)/FZINB/")
source("packageimport.r")
dyn.load("../src/projsplx_R.so")

load(file="../data/Buettner.RData")
load(file="../data/Kolodziejczyk.RData")
load(file="./data/Pollen.RData")
load(file="./data/Usoskin.RData")
#setwd("./R")
#sapply(list.files(), source)
#setwd("../")

library(splatter)
library(checkmate)
library(MCMCpack)
library(stringi)
#detach("package:splatter", unload=TRUE)

Test <- Buettner
counts <- round(10^Test$in_X-1)
params <- splatEstimate(counts)

getpara <- getParams(params,c("mean.rate","mean.shape","out.prob","out.facLoc","out.facScale"))


probs <- as.vector(rdirichlet(1,c(7,7,7)))
params <- setParams(params, group.prob = probs,de.prob = 0.2,dropout.present = TRUE)



sim <- splatSimulate(params, nGenes = 5000, method = "groups", batchCells = 200)


sim <- normalise(sim)
plotPCA(sim, colour_by = "Group")


##############
sim.groups <- splatSimulate(newSplatParams(), group.prob = c(0.2, 0.3, 0.5), method = "groups",
                            verbose = FALSE)
sim.groups <- normalise(sim.groups)
plotPCA(sim.groups, colour_by = "Group")

##############

X_counts <- counts(sim.groups)
X_counts <- X_counts[-which(rowSums(X_counts)==0),]
n_gene <- nrow(X_counts)
n_cell <- ncol(X_counts)
prior.zinb(X_counts)
X <- log(X_counts+1,10)
cores.ratio <- 6
maxK <- 51
cl <- start_cluster(cores.ratio)
FK <-FZINB.matrix(X_counts, cl = cl, calc.dists = FALSE, LogCOUNTS = FALSE)
stopCluster(cl)
cl <- start_cluster(2)
GK <- multiple.kernel(t(X), cl = cl , calc.dists = FALSE)[[maxK]]
stopCluster(cl)


true_lab <- as.numeric(sapply(sim$Group, function(x) stri_sub(x,-1,-1)))
result <- SIMLR(X = X, c = getParams(params, "nGroups")$nGroups, cores.ratio = cores.ratio)
compare(true_lab,result$y$cluster, method = "nmi")



nmi <- list()
n_clusters <- length(probs)
Kmeans_FZINB <- kmeans(FK, n_clusters, nstart = 30)$cluster
Kmeans_Gaussian <- kmeans(GK, n_clusters, nstart = 30)$cluster

nmi[[1]] <- compare(Kmeans_Gaussian,true_lab,method = "nmi")
nmi[[2]] <- compare(Kmeans_FZINB,true_lab,method = "nmi")

clustered <- sort(true_lab,index.return = TRUE)$ix

melted_FZINB <- melt(as.matrix(FK)[clustered,clustered])
melted_Gaussian <- melt(as.matrix(GK)[clustered,clustered])
pp<-0
for (melted in list(melted_Gaussian, melted_FZINB)){
  limits <- range(quantile(melted$value,c(0,0.99)))
  
   ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar", limits=limits)+
    labs(y ="",
         caption = paste("NMI ",signif(nmi[[1]],3) ) )+
    #theme_void()+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          aspect.ratio=1,
          legend.title=element_blank())+
    scale_y_continuous(trans = "reverse")
  
}






sim <- splatSimulate(params, nGenes = 5000)

source("../FZINB_kernel/R/functions.r")



Theta <- prior.zinb(X_counts)$Theta
mle <- MLE.zinb(X_counts, Theta = Theta)
MLE <- data.frame(mu=mle[,1],theta=mle[,2],pi=mle[,3])
xmax <- quantile(MLE$mu,.95)
ymax <- quantile(MLE$theta,.95)
nn <- length(intersect(which(MLE$mu < xmax),which(MLE$theta < ymax)))
size <- nrow(MLE)
ggplot(MLE,aes(x=mu,y=theta,colour = pi)) +
  geom_point() + 
  scale_colour_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar",limits=c(0,1)) +
  xlim(c(0,xmax)) +
  ylim(c(0,ymax)) +
  labs(title = "name",
       x = expression(hat(mu)),
       y = expression(hat(theta)),
       caption = paste(signif(nn/size*100,3),"% shown")) +
  theme(plot.title = element_text(hjust = 0.5))



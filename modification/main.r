source("SIMLRdataset.r")
#source("DUOdataset.R")
source("importpackage.R")
##############


plotPCA(sim.groups, colour_by = "Group")
X_counts <- counts(sim.groups)
if (any(which(rowSums(X_counts)==0))){
  X_counts <- X_counts[-which(rowSums(X_counts)==0),]
}
n_gene <- nrow(X_counts)
n_cell <- ncol(X_counts)
Theta0 <- prior.zinb(X_counts)$Theta
X <- log(X_counts+1,10)
cores.ratio <- 6

### MLE scattering

thet <- MLE.zinb(X_counts,Theta0,invisible =1)

MLE  <- data.frame(mu=thet[,1],theta=thet[,2],pi=thet[,3])
xmax <- quantile(MLE$mu,.95)
ymax <- quantile(MLE$theta,.95)
nn <- length(intersect(which(MLE$mu < xmax),which(MLE$theta < ymax)))
size <- nrow(MLE)
p1 <- ggplot(MLE,aes(x=mu,y=theta,colour = pi)) +
  geom_point() + 
  scale_colour_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar",limits=c(0,1)) +
  xlim(c(0,xmax)) +
  ylim(c(0,ymax)) +
  labs(title = "params estimated by Buettner",
       x = expression(hat(mu)),
       y = expression(hat(theta)),
       caption = paste(signif(nn/size*100,3),"% shown")) +
  theme(plot.title = element_text(hjust = 0.5))
plot(p1)

#plot mu against pi

xmax <- quantile(MLE$mu,.95)
ymax <- 1
nn <- length(which(MLE$mu < xmax))
size <- nrow(MLE)

p2 <- ggplot(MLE,aes(x=mu,y = pi)) +
  geom_point(shape=18, color="blue") + 
  xlim(c(0,xmax)) +
  ylim(c(0,ymax)) +
  labs(title = "name",
       x = expression(hat(mu)),
       y = expression(hat(pi)),
       caption = paste(signif(nn/size*100,3),"% shown")) +
  theme(plot.title = element_text(hjust = 0.5))
plot(p2)

### heatmap of similarity matrix

p <- list()
nmi <- list()

true_lab <- as.numeric(factor(colData(sim.groups)$Group))
result <- SIMLR(X = X, c = getParams(params, "nGroups")$nGroups, cores.ratio = cores.ratio)
nmi[[1]] <- compare(true_lab,result$y$cluster, method = "nmi")

setwd("../../SIMLR(John)/FZINB/")
cl <- start_cluster(cores.ratio)
FK <-FZINB.matrix(X_counts, cl = cl, calc.dists = FALSE, LogCOUNTS = FALSE)
stopCluster(cl)
setwd("../../splatter_simulation/modification/")

n_clusters <- length(probs)
Kmeans_FZINB <- kmeans(FK, n_clusters, nstart = 30)$cluster
nmi[[2]] <- compare(Kmeans_FZINB,true_lab,method = "nmi")

name <- list("SIMLR Gaussian kernel","FZINB kernel")
clustered <- sort(true_lab,index.return = TRUE)$ix

melted_FZINB <- melt(as.matrix(FK)[clustered,clustered])
melted_Gaussian <- melt(result$S[clustered,clustered])
pp<-1
for (melted in list(melted_Gaussian, melted_FZINB)){
  limits <- range(quantile(melted$value,c(0,0.99)))
  
   p[[pp]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar", limits=limits)+
    labs(y ="",
         title = name[[pp]],
         caption = paste("NMI ",signif(nmi[[pp]],3) ) )+
    #theme_void()+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          aspect.ratio=1,
          legend.title=element_blank())+
    scale_y_continuous(trans = "reverse")
   pp<-pp+1
  
}

p[[2]]




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



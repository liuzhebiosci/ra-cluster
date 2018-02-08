###############################################################################
# Readme
###############################################################################
# 1. install the relevant packages using the following codes
#    source("http://bioconductor.org/biocLite.R")
#    biocLite(c("ConsensusClusterPlus", "mclust", "fpc", "cluster"))
# 2. set the local directory (local.dir) and experiment file (exprs.f)
#    local.dir="/zheliu/UCB/workspace/hand-over code/gene-level/gse17755/"
#    exprs.f="/zheliu/UCB/workspace/hand-over code/gene-level/gse17755.exprs.txt"
# 
# A “gene level” directory containing the gene level analysis results is obtained after running the script. 
# 1. "clustering" folder includes the consensus clustering results usingusing hierarchical clustering, Partition Around Mediod algorithms and correlation coefficient, Euclidean distance. 
#    The interpretation of the results can be found from ConsensusClusterPlus package http://www.bioconductor.org/packages/2.12/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf. 
# 2. “silhouette plot” folder includes the silhouette plot for different algorithms and distance metrics from 2 to 20 clusters. 
#    The interpretation of the silhouette plot can be found from fpc package http://cran.r-project.org/web/packages/fpc/fpc.pdf.
# 3. “descriptive statistics” folder contains the descriptive statistics. 
#    The full description of these statistics can be found from fpc package http://cran.r-project.org/web/packages/fpc/fpc.pdf. 
# 4. “BIC plot” file contains the Bayesian Information Criteria analysis of the gene expression data. 
#    The interpretation of the results can be found from mclust package http://cran.r-project.org/web/packages/mclust/vignettes/mclust.pdf.


rm(list=ls())

# set the local directory
local.dir="/Users/zheliu/Desktop/R1/gene-level/"

# load the gene expression data
exprs.f="/Users/zheliu/Desktop/R1/gene-level/gse17755.exprs.txt"

setwd(local.dir)

dir.create(file.path(getwd(), "gse17755"), showWarnings = FALSE)
dir.create(file.path(getwd(), "gse17755/clustering"), showWarnings = FALSE)

res.f=paste(local.dir, "gse17755", sep="")

library(ConsensusClusterPlus)
library(mclust)
library(fpc)
library(cluster)

#use the genes with top 5000 MADs for clustering
d.raw=as.matrix(read.table(exprs.f, sep="\t"))

#normalize on gene level into zero mean and unit variance
d=t(apply(d.raw, 1, function(x) (x-mean(x))/sd(x)))
d=round(d,3)

#compute the distance metric
d_pearson = as.dist(1-cor(d,method="pearson"))
dt<-dist(t(d), method = "euclidean", p=2)

##############################################################
# The following code is used to perform consensus clustering 
# of the gene expression data
##############################################################

print("consensus clustering--running")

#perform hierachical clustering with correlation coefficient
title.hc.cor = paste(getwd(), "/gse17755/clustering/hc.cor", sep="")
res.hc.cor = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,title=title.hc.cor,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")
icl.hc.cor = calcICL(res.hc.cor,title=title.hc.cor,plot="pdf")

#perform hierachical clustering with euclidean distance
title.hc.euc = paste(getwd(), "/gse17755/clustering/hc.euclidean", sep="")
res.hc.euc = ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,title=title.hc.euc,clusterAlg="hc",distance="euclidean",seed=1262118388.71279,plot="pdf")
icl.hc.euc = calcICL(res.hc.euc,title=title.hc.euc,plot="pdf")

#perform pam with correlation coefficient
title.pam.cor = paste(getwd(), "/gse17755/clustering/pam.pearson", sep="")
res.pam.cor = ConsensusClusterPlus(d_pearson,maxK=20,reps=1000,pItem=0.8,pFeature=1,title=title.pam.cor,clusterAlg="pam",distance="pearson",seed=1262118388.71279,plot="pdf",verbose=F)
icl.pam.cor = calcICL(res.pam.cor,title=title.pam.cor,plot="pdf")

#perform pam with euclidean distance
title.pam.euc = paste(getwd(), "/gse17755/clustering/pam.euclidean", sep="")
res.pam.euc = ConsensusClusterPlus(dt,maxK=20,reps=1000,pItem=0.8,pFeature=1,title=title.pam.euc,clusterAlg="pam",distance="euclidean",seed=1262118388.71279,plot="pdf",verbose=F)
icl.pam.euc = calcICL(res.pam.euc,title=title.pam.euc,plot="pdf")

print("consensus clustering--done")

##############################################################
# The following code is used to perform silhouette plot of the
# different clustering results
##############################################################

setwd(res.f)
dir.create(file.path(getwd(), "/silhouette plot"), showWarnings = FALSE)

print("Silhouette plot--running")

#results using hierachical clustering with correlation coefficient
for(i in 2:20){
  pdf(file = paste(res.f, "/silhouette plot/hc.cor.Silhouette plot", i, ".pdf", sep=""))
  si <- silhouette(res.hc.cor[[i]][["consensusClass"]], dist=d_pearson)
  ssi <- summary(si)
  plot(si, main="Silhouette plot")
  dev.off()
}

#results using hierachical clustering with euclidean distance
for(i in 2:20){
  pdf(file = paste(res.f, "/silhouette plot/hc.eud.Silhouette plot", i, ".pdf", sep=""))
  si <- silhouette(res.hc.euc[[i]][["consensusClass"]], dist=dt)
  ssi <- summary(si)
  plot(si, main="Silhouette plot")
  dev.off()
}

#results using pam with correlation coefficient
for(i in 2:20){
  pdf(file = paste(res.f, "/silhouette plot/pam.cor.Silhouette plot", i, ".pdf", sep=""))
  si <- silhouette(res.pam.cor[[i]][["consensusClass"]], dist=d_pearson)
  ssi <- summary(si)
  plot(si, main="Silhouette plot")
  dev.off()
}

#results using pam with euclidean distance
for(i in 2:20){
  pdf(file = paste(res.f, "/silhouette plot/pam.eud.Silhouette plot", i, ".pdf", sep=""))
  si <- silhouette(res.pam.euc[[i]][["consensusClass"]], dist=dt)
  ssi <- summary(si)
  plot(si, main="Silhouette plot")
  dev.off()
}

print("Silhouette plot--done")

####################################################################
#plot the descriptive statistics for the consensus clustering results
####################################################################
setwd(res.f)
dir.create(file.path(getwd(), "descriptive statistics"), showWarnings = FALSE)

print("descriptive statistics--running")

tryCatch({
  hc.cor.list=list()
  hc.eud.list=list()
  pam.cor.list=list()
  pam.eud.list=list()
  
  hc.cor.res=cluster.stats(d_pearson, res.hc.cor[[2]]$consensusClass)
  stats.len=length(hc.cor.res)
  
  k=1
  
  for(j in 1:stats.len){
    if(length(hc.cor.res[[j]])==1){
      k=k+1
      
      hc.cor.list[[k]]=""
      hc.eud.list[[k]]=""
      pam.cor.list[[k]]=""
      pam.eud.list[[k]]=""
      
      names(hc.cor.list)[k]=names(hc.cor.res)[j]
      names(hc.eud.list)[k]=names(hc.cor.res)[j]
      names(pam.eud.list)[k]=names(hc.cor.res)[j]
      names(pam.cor.list)[k]=names(hc.cor.res)[j]
      
      for(i in 2:20){
        hc.cor.res=cluster.stats(d_pearson, res.hc.cor[[i]]$consensusClass)
        hc.eud.res=cluster.stats(dt, res.hc.euc[[i]]$consensusClass)
        pam.cor.res=cluster.stats(d_pearson, res.pam.cor[[i]]$consensusClass)
        pam.eud.res=cluster.stats(dt, res.pam.euc[[i]]$consensusClass)
        
        if(length(hc.cor.res[[j]])==1){
          
          #hc.cor consensus clustering
          hc.cor.list[[k]][[i]]=hc.cor.res[[j]]
          names(hc.cor.list[[k]])[[i]]=names(hc.cor.res)[j]
          
          #pam.eud consensus clustering
          pam.eud.list[[k]][[i]]=pam.eud.res[[j]]
          names(pam.eud.list[[k]])[[i]]=names(pam.eud.res)[j]
          
          #hc.eud consensus clustering
          hc.eud.list[[k]][[i]]=hc.eud.res[[j]]
          names(hc.eud.list[[k]])[[i]]=names(hc.eud.res)[j]
          
          #pam.cor consensus clustering
          pam.cor.list[[k]][[i]]=pam.cor.res[[j]]
          names(pam.cor.list[[k]])[[i]]=names(pam.cor.res)[j]
        }
      }
    }
  }
  
  for(j in 1:length(hc.cor.list)){
    if(length(hc.cor.list[[j]])>1){
      pdf(file=paste(res.f, "/descriptive statistics/consensus hc.cor clustering statistic.", names(hc.cor.list)[j],".pdf", sep=""))
      plot(1:20, hc.cor.list[[j]], xlab="number of clusters", ylab="value", main=names(hc.cor.list)[j])
      dev.off()
      
      pdf(file=paste(res.f, "/descriptive statistics/consensus pam.eud clustering statistic.", names(hc.cor.list)[j],".pdf", sep=""))
      plot(1:20, pam.eud.list[[j]], xlab="number of clusters", ylab="value", main=names(hc.cor.list)[j])
      dev.off()
      
      pdf(file=paste(res.f, "/descriptive statistics/consensus pam.cor clustering statistic.", names(hc.cor.list)[j],".pdf", sep=""))
      plot(1:20, pam.cor.list[[j]], xlab="number of clusters", ylab="value", main=names(hc.cor.list)[j])
      dev.off()
      
      pdf(file=paste(res.f, "/descriptive statistics/consensus hc.eud clustering statistic.", names(hc.cor.list)[j],".pdf", sep=""))
      plot(1:20, hc.eud.list[[j]], xlab="number of clusters", ylab="value", main=names(hc.cor.list)[j])
      dev.off()
    }
  }
}, error = function(err){
  print(names(hc.cor.list)[j])
  print("error here")
})

print("descriptive statistics--done")

#######################################
#model based clustering and compute BIC
######################################
setwd(res.f)
dir.create(file.path(getwd(), "BIC plot"), showWarnings = FALSE)

print("Bayesian Information Criterion plot--running")

d.clust<- list()
d.clust[[1]]<- mclustBIC(d, modelNames = "EII")
d.clust[[2]]<- mclustBIC(d, modelNames = "VII")
d.clust[[3]]<- mclustBIC(d, modelNames = "EEI")
d.clust[[4]]<- mclustBIC(d, modelNames = "VEI")
d.clust[[5]]<- mclustBIC(d, modelNames = "EVI")
d.clust[[6]]<- mclustBIC(d, modelNames = "VVI")
d.clust[[7]]<- mclustBIC(d, modelNames = "EEE")
d.clust[[8]]<- mclustBIC(d, modelNames = "EEV")
d.clust[[9]]<- mclustBIC(d, modelNames = "VEV")
d.clust[[10]]<- mclustBIC(d, modelNames = "VVV")

d.clust.plot<- do.call(cbind, d.clust)

model.name<- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV")

max_y <- max(na.omit(d.clust.plot))*0.9
min_y <- min(na.omit(d.clust.plot))*1.1

# plot the BIC plot
pdf(file = paste(res.f, "/BIC plot/mclust BIC plot.pdf", sep=""))

plot(seq(1,9), d.clust.plot[, 1], ylim=c(min_y,max_y), pch=0, xlab="Number of components", ylab="BIC", col=1)
lines(seq(1,9), d.clust.plot[, 1], lty=1)

points(seq(1,9), d.clust.plot[, 2], col=2, pch=1)
lines(seq(1,9), d.clust.plot[, 2], lty=1)

points(seq(1,9), d.clust.plot[, 3], col=3, pch=2)
lines(seq(1,9), d.clust.plot[, 3], lty=1)

points(seq(1,9), d.clust.plot[, 4], col=4, pch=3)
lines(seq(1,9), d.clust.plot[, 4], lty=1)

points(seq(1,9), d.clust.plot[, 5], col=5, pch=4)
lines(seq(1,9), d.clust.plot[, 5], lty=1)

points(seq(1,9), d.clust.plot[, 6], col=6, pch=5)
lines(seq(1,9), d.clust.plot[, 6], lty=1)

points(seq(1,9), d.clust.plot[, 7], col=7, pch=6)
lines(seq(1,9), d.clust.plot[, 7], lty=1)

points(seq(1,9), d.clust.plot[, 8], col=8, pch=7)
lines(seq(1,9), d.clust.plot[, 8], lty=1)

points(seq(1,9), d.clust.plot[, 9], col=9, pch=8)
lines(seq(1,9), d.clust.plot[, 9], lty=1)

points(seq(1,9), d.clust.plot[, 10], col=10, pch=9)
lines(seq(1,9), d.clust.plot[, 10], lty=1)

legend("bottomright", model.name,col=1:10,pch=0:9, ncol=2)
dev.off()

write.table(d.clust.plot, file="BIC plot/BIC.txt")

print("Bayesian Information Criterion plot--done")

#save the entire environment
save.image(file="gse17755.Rdata")

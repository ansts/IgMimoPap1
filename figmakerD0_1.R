
# Fig. 4 A 
require("lattice")
require(Biostrings)
require(plot3D)
pdf(file="genPWM.pdf", width=3, height=5)
par(mfrow=c(1,1))
generalPSSM=consensusMatrix(ODXsel[,1])
generalPSPM=apply(generalPSSM,2,function(i){i/224087})
generalPSSM=generalPSSM[names(aaPh7_ok),]
generalPSPM=generalPSPM[names(aaPh7_ok),]
generalaaLO=log2(rowSums(generalPSSM)/224087/7/aaPh7_ok)
generalPWM=log2(generalPSPM/matrix(rep(aaPh7_ok,7), nrow = 20))
generalPWM=generalPWM[order(generalaaLO),]
barplot(generalaaLO[order(generalaaLO)], horiz = TRUE, col=1,cex.names=0.5)
levelplot(t(generalPWM), cuts=100, col.regions=jet.col(1280),  xlab="Position", ylab="Amino Acid")
dev.off()

# Suppl. Fig. 2
require(gplots)
pdf(file="hmodx79030pwmx.pdf",width=14, height=14)
cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(odx79030pwmx,distfun = dist, hclustfun=hclwrd, scale="none",col=cpl1,  trace="none", cexRow = 0.01, cexCol =0.5)
dev.off()

# Tsne representations

require(Rtsne)
require(parallel)
zscore=read.csv(file="zsc.csv", header = TRUE)
rownames(zscore)=zscore$X
zscore=zscore[,-1]
ODXs790pzsc=t(apply(ODXs_790_30,1,function(l){
  p=strsplit(l[1], split="")
  unlist(zscore[unlist(p),])
}))
rownames(ODXs790pzsc)=ODXs_790_30[,1]

rndpep1u=unique(rndpep1)
rndp1zsc=t(sapply(rndpep1u,function(l){
  p=strsplit(l, split="")
  unlist(zscore[unlist(p),])
}))
rownames(rndp1zsc)=rndpep1u

mixpep=c(sample(ODXs_790_30[,1],50000),sample(rndpep1u, 50000))
mixpep=unique(mixpep)
mixpep=c(mixpep,pepnegrnd)
mixpep=unique(mixpep)
mixpeplab=1*(mixpep %in% ODXs_790_30[,1])+2*(mixpep %in% rndpep1u)+4*(mixpep %in% pepnegrnd)
mixpep=mixpep[mixpeplab %in% c(1,2,4)]
mixpeplab=1*(mixpep %in% ODXs_790_30[,1])+2*(mixpep %in% rndpep1u)+4*(mixpep %in% pepnegrnd)
mixpzsc=t(sapply(mixpep,function(l){
  p=strsplit(l, split="")
  unlist(zscore[unlist(p),])
}))

zsc=list(ODXs790pzsc,rndp1zsc, mixpzsc)

cl <- makeCluster(3)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterEvalQ(cl, {library(Rtsne)})
clusterExport(cl, ex)
clusterExport(cl, list("zsc"), envir=environment())
proct=proc.time()
tsne20zsc=parLapply(cl,zsc,function(Z){Rtsne(Z, initial_dims = 20, perplexity = 50, max_iter = 1500)}) 
stopCluster(cl)
names(tsnezsc)=c("ODXs790pzsc","rndp1zsc", "mixpzsc")
print(proc.time()-proct)

#
# Fig. 6A
#

pdf(file="tsneODXs790azsc.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$ODXs790pzsc$Y, pch=16, cex=0.1, col=rgb(0.5,0.5,0.5,0.2),xlim=c(-30,30), ylim=c(-30,30))
par(new=TRUE)
plot(tsnezsc$ODXs790pzsc$Y[fl790_pep5,], pch=16, col=reix790_5[fl790_pep5], cex=0.75, xlim=c(-30,30), ylim=c(-30,30))
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Selected Mimotopes with Five Highly Significant Clusters",side=3, cex=2, line=5)
par(opar)
dev.off()

#
# Fig. 6B
#

pdf(file="tsnerndp1zsc.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$rndp1zsc$Y, pch=16, cex=0.1, col=rgb(0.2,0.2,0.2,0.2),xlim=c(-30,30), ylim=c(-30,30), xlab="D1", ylab="D2")
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()

#
# Fig.7
#

pdf(file="tsneODXs790pepposzsc.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$ODXs790pzsc$Y, pch=16, cex=0.05, col=rgb(0.5,0.5,0.5,0.2),xlim=c(-30,30), ylim=c(-30,30))
par(new=TRUE)
plot(tsnezsc$ODXs790pzsc$Y[fl790_peppos,], pch=16, col=rgb(1,0,0,0.75), cex=0.75, xlim=c(-30,30), ylim=c(-30,30))
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Selected Mimotopes with the Optimized Library",side=3, cex=2, line=5)
par(opar)
dev.off()

#
# Fig. 8
#

pdf(file="tsnemixpazsc_RB.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$mixpzsc$Y, pch=16, col=rgb(0.1,0.1,0.1,0.3), cex=0.5, xlim=c(-40,40), ylim=c(-40,40), xlab="D1")
par(new=TRUE)
plot(tsnezsc$mixpzsc$Y[mixpeplab==1,], pch=16, col=rgb(1,0,0,0.3), cex=0.5, xlim=c(-40,40), ylim=c(-40,40), xlab="", ylab="",xaxt="n",yaxt="n")
par(new=TRUE)
plot(tsnezsc$mixpzsc$Y[mixpeplab==4,], pch=16, col=rgb(0,0,1,0.3), cex=0.5, xlim=c(-40,40), ylim=c(-40,40), xlab="", ylab="",xaxt="n",yaxt="n")
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mixture of Selected Mimotopes and Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()



#
# Suppl. Fig. 5
#

require(factoextra)
pdf(file="pcaODXs.pdf",width=4, height=4)
par(mai=c(1,1,1,1))
barplot(get_eig(pcaODxs)[,2], las=2, ylab="Percentage of Variance Explained",cex.axis=0.5,cex.lab=0.75)
par(new=TRUE)
plot(1:35,get_eig(pcaODxs)[,3], yaxt="n", pch=16, cex=0.3, ylim=c(0,100), ylab="", xlab="Number of Dimensions",cex.axis=0.5, cex.lab=0.75)
lines(1:35,get_eig(pcaODxs)[,3])
z=(1:10)*10
axis(4, at=z,labels=z, las=2,cex.axis=0.5)
par(new=FALSE)
dev.off()

#
# Suppl. Fig. 6
#

pdf(file="ClustersReprsnt.pdf",width=4, height=4)
opar=par(mai=c(1,1,0.5,0.5))
plot(mixtab[,1], mixtab[,2], col=mix2Nfls[,1]*1+1, pch=16, cex=0.5, cex.axis=0.5, cex.lab=0.75, xlab="Number of Mimotopes", ylab="Number of Random Peptides")  #xaxt="n", yaxt="n"
box(lwd=2)
#axis(1,lwd=2,lwd.ticks = 1,cex=0.25);axis(2,lwd=2, lwd.ticks = 1,cex=0.25);axis(3,lwd=2,lwd.ticks = 0, labels = FALSE);axis(4,lwd=2,lwd.ticks = 0, labels = FALSE)
par(opar)
dev.off()
pdf(file="Suppl_LOHIST.pdf", width=12, height=7.2)
#par(mfrow=c(2,2))

xr=range(odx79030pwmx)
xr[1]=floor(xr[1])
xr[2]=ceiling(xr[2])
brp=seq(from=xr[1], to=xr[2], length.out=60)
odx79030pwmx_hist=apply(odx79030pwmx,2,function(i){hist(i, breaks=brp, plot = FALSE)$density})
image(log(t(odx79030pwmx_hist[,order(diverg, decreasing = TRUE)])+0.1),col=gray(750:0/750,1),axes=FALSE, ylab="LO")
axis(1,at=(0:139)/139, labels = colnames(odx79030pwmx)[order(diverg, decreasing = TRUE)], las=2, cex.axis=0.5)
axis(2,at=(0:12)/12, labels = round(seq(xr[1],xr[2],length.out = 13),digits = 2), las=2, cex.axis=0.7)

dev.off()

pdf(file="3D_LO.pdf", width=8, height=8)
require(plot3D)
diverg=apply(odx79030pwmx,2,function(i){sum(abs(i))/790})
divergSD=apply(odx79030pwmx,2,sd)
divergMn=apply(odx79030pwmx,2,mean)
#scatter3D(divergMn,divergSD,diverg, phi=45, theta=45, cex=0, xlab="Mean",ylab="SD", zlab="SUM(|LO|)")
text3D(divergMn,divergSD,diverg, phi=5, theta=70, cex=0.75, xlab="Mean(LO)",ylab="SD(LO)", zlab="SUM(|LO|)", labels = names(diverg), colvar=divergMn, col=jet.col(32), clokey=TRUE, cex.text =0.2, clab="Mean(LO)", bty="u", col.panel=gray(0,.45))
dev.off()

require("lattice")
require(Biostrings)
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

pdf(file="hmodx79030pwmx.pdf",width=14, height=14)
cpl1=colorRampPalette(c("black","black","#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))(n=128)
heatmap.2(odx79030pwmx,distfun = dist, hclustfun=hclwrd, scale="none",col=cpl1,  trace="none", cexRow = 0.01, cexCol =0.5)
dev.off()


# Tsne presentations

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

proct=proc.time()

cl <- makeCluster(4)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterEvalQ(cl, {library(Rtsne)})
clusterExport(cl, ex)
clusterExport(cl, list("ODXs790pzsc"), envir=environment())
tsneODXs790zsc=parLapply(cl,c(50,60,70,80),function(n){Rtsne(ODXs790pzsc, initial_dims = 14, perplexity = n)}) #takes about 90 min at theta=0.5
stopCluster(cl)
print(proc.time()-proct)

fl790_peppos=ODXs_790_30[,1] %in% peppostive
fl790_pep5=ODXs_790_30[,1] %in% pep_top_5
fl790_pepother5=ODXs_790_30[,1] %in% pep_other_5
ixpep5=c(2,6,9,10,11)
ixpepother5=c(115,61,55,53, 258)
ixpep5=cbind(ixpep5,c(1,2,3,4,5))
ixpepother5=cbind(ixpepother5,c(1,2,3,4,5))
reix790_5=sapply(ODXs_790_30[,2], function(i){
  if ((i+1) %in% ixpep5[,1]) res=ixpep5[ixpep5[,1]==(i+1),2] else res=0
  return(res)
})
reix790_o5=sapply(ODXs_790_30[,2], function(i){
  if ((i+1) %in% ixpepother5[,1]) res=ixpepother5[ixpepother5[,1]==(i+1),2] else res=0
  return(res)
})
fl790_pepneg=ODXs_790_30[,1] %in% pepnegtve
fl790_pepneglo=ODXs_790_30[,1] %in% pepnegtlo
for (i in 1:4){
  plot(tsneODXs790zsc[[i]]$Y, pch=16, cex=0.05, col="gray",xlim=c(-25,25), ylim=c(-25,25), main=c(50,60,70,80)[i])
  par(new=TRUE)
  plot(tsneODXs790zsc[[i]]$Y[fl790_pep5,], pch=16, col=reix790_5[fl790_pep5], cex=0.5, xlim=c(-25,25), ylim=c(-25,25))
  plot(tsneODXs790zsc[[i]]$Y, pch=16, cex=0.05, col="gray",xlim=c(-25,25), ylim=c(-25,25))
  par(new=TRUE)
  plot(tsneODXs790zsc[[i]]$Y[fl790_pepneg,], pch=16, col=4, cex=0.5, xlim=c(-25,25), ylim=c(-25,25))
  par(new=TRUE)
  plot(tsneODXs790zsc[[i]]$Y[fl790_peppos,], pch=16, col=2, cex=0.5, xlim=c(-25,25), ylim=c(-25,25))
}
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
tsnezsc=parLapply(cl,zsc,function(Z){Rtsne(Z, initial_dims = 14, perplexity = 50, max_iter = 1500, theta=0.3)}) 
stopCluster(cl)
names(tsnezsc)=c("ODXs790pzsc","rndp1zsc", "mixpzsc")
print(proc.time()-proct)

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

tsnCLsmpl=c("1","62","104","138","407","447","481","550","706","778")
#fltsnClsmpl=sapply(kmtsneODX$cluster,function(c){if (c %in% tsnCLsmpl) rank(tsnCLsmpl)[tsnCLsmpl==c] else 0})
table(ODXs_790_30[fltsnClsmpl>0,2],fltsnClsmpl[fltsnClsmpl>0])

require(RFast)
bx=c(0,5:25,35)
n=nrow(ODXs790pzsc)

proct=proc.time()
cl <- makeCluster(7)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterEvalQ(cl, {library(matrixStats); library(Rfast)})
clusterExport(cl, ex)
clusterExport(cl, list("ODXs790pzsc","bx","n"), envir=environment())
dbin=parSapply(cl,1:n, function(i){     #nrow(ODXs790pzsc)
  di=dista(t(ODXs790pzsc[i,]),ODXs790pzsc)
  binCounts(di,bx=bx)/n 
})
stopCluster(cl)
print(proc.time()-proct)
dbinc=t(apply(dbin,2,function(c){sapply(1:22,function(r){log(sum(c[1:r]))})}))
dbins=rowSums(dbin)
cdbins=sapply(1:22, function(i){sum(dbins[1:i])})
line(log((5:11)),log(cdbins[1:7]))
adbinc=apply(dbinc,1,function(r){if (all(is.infinite(r))) 0 else {a=line(log(5:11),r[1:7]);a$coefficients[2]}})
hist(adbinc)

n=nrow(rndp1zsc)
proct=proc.time()
cl <- makeCluster(7)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterEvalQ(cl, {library(matrixStats); library(Rfast)})
clusterExport(cl, ex)
clusterExport(cl, list("rndp1zsc","bx","n"), envir=environment())
dbrin=parSapply(cl,1:n, function(i){     
  di=dista(t(rndp1zsc[i,]),rndp1zsc)
  binCounts(di,bx=bx)/n 
})
stopCluster(cl)
print(proc.time()-proct)
dbrinc=t(apply(dbrin,2,function(c){sapply(1:22,function(r){log(sum(c[1:r]))})}))
dbrins=rowSums(dbrin)
cdbrins=sapply(1:22, function(i){sum(dbrins[1:i])})
line(log(5:11),log(cdbrins[1:7]))
adbrinc=apply(dbrinc,1,function(r){if (all(is.infinite(r))) 0 else {a=line(log(5:11),r[1:7]);a$coefficients[2]}})
hist(adbinc)


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

pdf(file="tsnerndp1zsc.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$rndp1zsc$Y, pch=16, cex=0.1, col=rgb(0.2,0.2,0.2,0.2),xlim=c(-30,30), ylim=c(-30,30), xlab="D1", ylab="D2")
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()

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

cl <- makeCluster(7)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cl, ex)
clusterExport(cl, list("tsnezsc"), envir=environment())
kmtsneODXbtwn=parSapply(cl,seq(10,3000, by=10),function(i){
  k=kmeans(tsnezsc$mixpzsc$Y,i, iter.max = 100, nstart = 10)
  return(c(i,k$betweenss/k$totss))
})
stopCluster(cl)
kmtsneODXbtwn=t(kmtsneODXbtwn)

plot((kmtsneODXbtwn[50:300,]), xlab="N Clusters", ylab="Proportion of SS")
kmtsneODX=kmeans(tsnezsc$mixpzsc$Y,350, iter.max = 350, nstart = 50)


pdf(file="tsnemixp1000cl.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$mixpzsc$Y, pch=16, cex=0.3, col=kmtsneODX$cluster,xlim=c(-40,40), ylim=c(-40,40))
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mixture of Selected Mimotopes and Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()

mixrep=(table(kmtsneODX$cluster,mixpeplab)[,1]/table(kmtsneODX$cluster,mixpeplab)[,2])
mixR=sapply(seq_along(kmtsneODX$cluster), function(i){mixrep[kmtsneODX$cluster[i]]})
hist(mixrep, breaks = 50)
mixtab=(table(kmtsneODX$cluster,mixpeplab)[,1:2])

plot(mixtab[,1], mixtab[,2])
mix2N=mvnormalmixEM(mixtab, k=2)
mix2Nfl=apply(mix2N$posterior,2,">",0.5)
mix2Nfls=apply(mix2N$posterior,2,">",0.9)
plot(mixtab[,1], mixtab[,2], col=mix2Nfls[,1]*1+1)
plot(mixtab[mix2Nfl[,1],1], mixtab[mix2Nfl[,1],2])
abline(lm(mixtab[mix2Nfl[,1],2]~mixtab[mix2Nfl[,1],1]))
mixfs=sapply(seq_along(kmtsneODX$cluster), function(i){mix2Nfls[kmtsneODX$cluster[i],1]})
underrepTable=table(mixfs, mixpeplab)
undRclpep=data.frame(pep=mixpep[mixpeplab==2 & mixfs],cluster=as.integer(kmtsneODX$cluster[mixpeplab==2 & mixfs]),stringsAsFactors= FALSE)
undRclpepL=split(undRclpep$pep,as.factor(undRclpep$cluster))

pdf(file="ClustersReprsnt.pdf",width=4, height=4)
opar=par(mai=c(1,1,0.5,0.5))
plot(mixtab[,1], mixtab[,2], col=mix2Nfls[,1]*1+1, pch=16, cex=0.5, cex.axis=0.5, cex.lab=0.75, xlab="Number of Mimotopes", ylab="Number of Random Peptides")  #xaxt="n", yaxt="n"
box(lwd=2)
#axis(1,lwd=2,lwd.ticks = 1,cex=0.25);axis(2,lwd=2, lwd.ticks = 1,cex=0.25);axis(3,lwd=2,lwd.ticks = 0, labels = FALSE);axis(4,lwd=2,lwd.ticks = 0, labels = FALSE)
par(opar)
dev.off()

pdf(file="tsnemixp1000clre.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$mixpzsc$Y, pch=16, cex=0.3, col=(mixfs)*1+2*(mixpeplab==4)+1,xlim=c(-40,40), ylim=c(-40,40))
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mixture of Selected Mimotopes and Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()

plot(table(kmtsneODX$cluster,mixpeplab)[,1],table(kmtsneODX$cluster,mixpeplab)[,2])

pdf(file="tsnemixp1000clN.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$mixpzsc$Y, pch=16, cex=0, xlim=c(-40,40), ylim=c(-40,40))
text(tsnezsc$mixpzsc$Y, pch=16, cex=0.05,col=(mixfs)*1+2*(mixpeplab==4)+1, labels=kmtsneODX$cluster, xlim=c(-40,40), ylim=c(-40,40))
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mixture of Selected Mimotopes and Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()

adjop=c(239,88,222,115,111,158,67)
undrrpcl=c(162,59,242,212,228,86,201)
undrrp2cl=c(90,78,336,287,304,186,160)
adjundrcl=c(57,169,68,206,182,7,51)
adjundrcl=c(57,169,68,206,182,7,51)
adjopcl=c(239,88,222,115,111,158,67)
other2cl=c(3,130,131,135,197,270,286)
other3cl=c(58,163,147,107,75,40,45)

for (i in other3cl){
  fnm=paste("other3cl",i,".csv", sep="")
  write.csv(mixpep[kmtsneODX$cluster==i],file=fnm)
}

for (i in undrrp2cl){
  fnm=paste("undrrp2cl",i,".csv", sep="")
  write.csv(mixpep[kmtsneODX$cluster==i],file=fnm)
}


selcl=c(undrrpcl,adjundrcl,undrrp2cl,adjop,other2cl,othercl,other3cl)

pdf(file="tsnemixpselcl.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsnezsc$mixpzsc$Y, pch=16, col=rgb(0.9,0.9,0.9,0.5), cex=0.3, xlim=c(-40,40), ylim=c(-40,40))
text(kmtsneODX$centers[selcl,], pch=16, cex=2,col=((1:49)<22)*1+1, labels=(1:350)[selcl], xlim=c(-40,40), ylim=c(-40,40))
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mixture of Selected Mimotopes and Random Peptides",side=3, cex=2, line=5)
par(opar)
dev.off()

pdf(file="pcaODXs.pdf",width=4, height=4)
par(mai=c(1,1,1,1))
barplot(get_eig(pcaODxs)[,2], las=2, ylab="Percentage of Variance Explained",cex.axis=0.5,cex.lab=0.75)
par(new=TRUE)
plot(1:35,get_eig(pcaODxs)[,3], yaxt="n", pch=16, cex=0.3, ylim=c(0,100), ylab="", xlab="Number of Dimensions",cex.axis=0.5, cex.lab=0.75)
lines(1:35,get_eig(pcaODxs)[,3])
axis(4, at=z,labels=z, las=2,cex.axis=0.5)
par(new=FALSE)
dev.off()

require(bamboo)
data(bamboo.validation.astral30)
bamastrval7=apply(bamboo.validation.astral30,1,function(r){
  l=lapply(1:(nchar(r[2])-6), function(i){
      s=unique(unlist(strsplit(substr(r[3],i,i+6), spli="")))
      if (length(unique(s))==1) c(s,substr(r[2],i,i+6)) else NA 
    })
  l=l[!is.na(l)]
  })
bamastrval7=t(as.data.frame(unlist(bamastrval7, recursive = FALSE)))
rownames(bamastrval7)=NULL
bam7=rbind(bamastrval7[bamastrval7[,1]=="C",],
           bamastrval7[bamastrval7[,1]=="E",][sample(length(bamastrval7[bamastrval7[,1]=="E",1]),2000),],
           bamastrval7[bamastrval7[,1]=="H",][sample(length(bamastrval7[bamastrval7[,1]=="H",1]),2000),],
           bamastrval7[bamastrval7[,1]=="T",])
bam7=unique(bam7)
bam7pep=bam7[,2]


SSexmplzsc=t(sapply(bam7pep,function(l){
  p=strsplit(l, split="")
  unlist(zscore[unlist(p),])
}))
rownames(SSexmplzsc)=bam7pep

ODXSS=unique(rbind(ODXs790pzsc,SSexmplzsc))
mixSS=unique(rbind(mixpzsc,SSexmplzsc))

zscSS=list(ODXSS,mixSS)
cl <- makeCluster(2)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterEvalQ(cl, {library(Rtsne)})
clusterExport(cl, ex)
clusterExport(cl, list("zscSS"), envir=environment())
proct=proc.time()
tsneSSzsc=parLapply(cl,zscSS,function(Z){Rtsne(Z, initial_dims = 14, perplexity = 50, max_iter = 1500)}) 
stopCluster(cl)
names(tsneSSzsc)=c("ODX","mix")
print(proc.time()-proct)

sso=sapply(rownames(ODXSS), function(p){ 
      if (p %in% bam7pep) (bam7[bam7[,2]==p,1]=="H")*2+(bam7[bam7[,2]==p,1]=="E")*3+(bam7[bam7[,2]==p,1]=="C")*4+(bam7[bam7[,2]==p,1]=="T")*5 else 1
})
ssm=sapply(rownames(mixSS), function(p){
      if (p %in% bam7pep) (bam7[bam7[,2]==p,1]=="H")*2+(bam7[bam7[,2]==p,1]=="E")*3+(bam7[bam7[,2]==p,1]=="C")*4+(bam7[bam7[,2]==p,1]=="T")*5 else 1
})
sso=unlist(sso)
ssm=unlist(ssm)

mixsspep=rownames(mixSS)
mixpsslab=1*(mixsspep %in% ODXs_790_30[,1])+2*(mixsspep %in% rndpep1u)+4*(mixsspep %in% pepnegrnd)

pdf(file="tsneODXSSzsc.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsneSSzsc$ODX$Y, pch=16, cex=0.1+(sso %in% c(2,3,5))*0.4, col=rgb((sso==2)*1+(sso==3)*1+(sso %in% c(1,4))*0.8,(sso==3)*1+(sso %in% c(1,4))*0.8,(sso==5)*1+(sso %in% c(1,4))*0.8,0.67),xlim=c(-30,30), ylim=c(-30,30), xlab="D1", ylab="D2")
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mimotopes and Secondary structure Examples",side=3, cex=2, line=5)
par(opar)
dev.off()


pdf(file="tsnemixSSzsc.pdf",width=14, height=14)
opar=par(mai=c(2,2,1.5,0.5), ann=FALSE)
plot(tsneSSzsc$mix$Y, pch=16, cex=(ssm %in% c(2,3,5)*0.2)+0.2, col=rgb((ssm==2)*1+(ssm==3)*1+(ssm %in% c(1,4))*(0.4+(mixpsslab>1)*0.5),(ssm==3)*1+(ssm %in% c(1,4))*(0.4+(mixpsslab>1)*0.5),(ssm==5)*1+(ssm %in% c(1,4))*(0.4+(mixpsslab>1)*0.5),0.67),xlim=c(-40,40), ylim=c(-40,40), xlab="D1", ylab="D2")
mtext("D1",side=1, cex=2, line=5)
mtext("D2",side=2, cex=2, line=5)
mtext("Mimotopes and Secondary structure Examples",side=3, cex=2, line=5)
par(opar)
dev.off()

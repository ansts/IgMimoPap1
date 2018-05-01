# "Pipeline" for the c10 group of chips
# gpraddNames adds group names for the peptides in the Name column of .gpr
# filepr preprocesses .gpr files in lf (adjust path) using gprLNSVM to locally 
#   normalize them by SVR, the fitted background is recorded 
#   as B635 columns.
# Use pepStat package to get pSets (c10pre – before staining, c10post – after 
# 	staining, c10pSet – difference) and chemical propertes 
# 	normalized pnSet (c10pnSet). Flagged spots are removed in c10pSetn.
# 
# The rest of the steps are explained in-line. After aggregating the technical 
# replicates and filtering low quality data at 2 stages (flagged bad and those 
# with missing replicates or relative SEM > 30%). 
# The data matrix is pxn (p= N of peptides, n= N of patients).
# The between array normalization is done using Cyclic Loess (method "affy")
# from package limma. The number of iterations is determined by optimizing the 
# total variance (sum of the diagonal of the covariance matrix).
# 

pfx="c10"
ptm=proc.time()
require(pepStat)
require(limma)
require(Biobase)
require(matrixStats)
require(d3heatmap)
require(multcomp)
require(mixtools)

lf=list.files(path = "c10/")
lapply(lf, function(x){gpraddNames(x, p="c10")})
filepr("c10_class/", sh=F)
mappre="mapc10pre.csv"
mappost="mapc10post.csv"
dirToParse="proc_c10_class"
#cn=c("YPYDVPDYAG", "DYKDDDDKAS")
c10pre=makePeptideSet(path=dirToParse, mapping.file = mappre, use.flags = TRUE, log=TRUE, bgCorrect.method = "normexp")
c10post=makePeptideSet(path=dirToParse, mapping.file = mappost, use.flags = TRUE, log=TRUE, bgCorrect.method = "normexp")
c10pSet=c10post
exprs(c10pSet)=exprs(c10post)-exprs(c10pre)
#c10pSet=c10pSet[!is.na(c10pSet@featureRange@elementMetadata@listData$featureID),]
c10v0=c10pSet@assayData$exprs
c10cnm=paste(c10pSet@phenoData@data$ptid,c10pSet@phenoData@data$diag, sep  = "_")
colnames(c10v0)=c10cnm

narows=apply(c10v0, 1, function(x){all(!is.na(x))}) # identify peptides with "bad" spots in any chip
c10pSetn=c10pSet[narows,]
c10pepts0=c10pSetn@featureRange@elementMetadata@listData$peptide
c10gr=unique(data.frame(class=c10pSetn@featureRange@elementMetadata@listData$featureID,pep=c10pepts0, stringsAsFactors = FALSE))

c10dupp=c10gr[duplicated(c10gr$pep),2]     # identify sequences that participate in more than one class - altogether 34
c10duppp=c10gr[c10gr$pep %in% c10dupp,]    # and reassign them to one of the classes using the rules encoded in reassdup()
c10dupreas=reassdup(c10duppp)
c10gr=c10gr[!c10gr$pep %in% c10dupp,]
c10gr=rbind(c10gr,c10dupreas)
c10gr=unique(c10gr)

c10pnSet=normalizeArray(c10pSetn, centered = FALSE) 
c10v=c10pnSet@assayData$exprs
c10pepts=c10pnSet@featureRange@elementMetadata@listData$peptide
c10ncol=ncol(c10v)
colnames(c10v)=c10cnm

pp=lapply(c10pepts,function(x){x})
c10vndf=DataFrame(cbind(c10v, pp))
c10nms=c(c10cnm, "pp")
c10vndf=data.frame(matrix(unlist(c10vndf), nrow=nrow(c10vndf)),stringsAsFactors = FALSE)
colnames(c10vndf)=c10nms

for (i in 1:10) c10vndf[,i]=as.numeric(c10vndf[,i])             # calculate SEM for the duplicates of each peptide
c10va=aggregate(as.matrix(c10vndf[,1:10])~unlist(c10vndf[,11]),c10vndf,"mean")
c10vasd=aggregate(as.matrix(c10vndf[,1:10])~unlist(c10vndf[,11]),c10vndf,"sd")
c10peptsn=c10va[,1]
rownames(c10va)=c10peptsn
c10ns=aggregate(as.matrix(c10vndf[,1:10])~unlist(c10vndf[,11]),c10vndf,"length")
c10va=c10va[,-1]
c10vasd=c10vasd[,-1]
c10ns=c10ns[,-1]

c10vsdco=c10vasd # Correction of s for low n
for (i in 1:10){
  for (j in 1: nrow(c10vasd)){
    if (c10ns[j,i]==1) co=1
    if (c10ns[j,i]==2) co=sqrt(2/pi)
    if (c10ns[j,i]==3) co=sqrt(pi)/2
    if (c10ns[j,i]==4) co=2*sqrt(2/(3*pi))
    if (c10ns[j,i]==5) co=3/4*sqrt(pi/2)
    if (c10ns[j,i]==6) co=8/3*sqrt(2/(5*pi))
    if (c10ns[j,i]>6) co=1
    c10vsdco[j,i]=c10vasd[j,i]/co
  }
}

c10semco=c10vsdco/sqrt(c10ns) 
c10semcop=c10semco/c10va  # use relative standard error of the mean to...
flc10va=apply(c10semcop,1,function(x){all(!is.na(x))&(max(x)<0.3)}) # filter low  quality data
c10vaco=c10va[flc10va,]
c10pep=rownames(c10vaco)


c10vn0=as.matrix(c10vaco) # find optimal n iter cyclic loess
n=100
d=matrix(0,nrow=n, ncol=2)
d[1,1]=1
d[1,2]=sum(diag(cov(c10vn0)))
for (i in 2:n) {
  c10vn=normalizeCyclicLoess(c10vn0, iterations=1, method="affy")
  d[i,1]=mean((c10vn-c10vn0)^2)
  d[i,2]=sum(diag(cov(c10vn)))
  print(i)
  c10vn0=c10vn
}
d=log10(d)
plot(d)
lines(d)
nit=which(d[,2]==min(d[,2]))
print(nit)

c10vn=normalizeCyclicLoess(c10vaco, iterations=7, method="affy")

c10Mean=rowMeans(c10vn)
c10SD=rowSds(c10vn)
c10vagr=lapply(c10pep, function(x){c10gr[c10gr$pep==x,1]})
c10MnSD=data.frame(unlist(c10Mean), unlist(c10SD), unlist(c10vagr), stringsAsFactors = FALSE)
colnames(c10MnSD)=c("Mean","SD","Class")

mus=list(c(3.5, 0.3), c(4.2, 0.45), c(5, 0.8))
dumcol=mvnormalmixEM(c10MnSD[,1:2], k=3,mu=mus) # find k gaussians in Mean x SD
c10vnn=c10vn-dumcol$mu[[1]][1]                  # subtract the mean mean of the lowest group as a puative negative 
fl10post=dumcol$posterior[,1]<0.05             # select the pepts whose mean is definitely positive by this criterion
plotMDS(c10vnn[fl10post,],col=as.double(dgn), xlab="dim 1", ylab="dim 2")

c10Mean=rowMeans(c10vnn)
c10SD=rowSds(c10vnn)
c10vagr=lapply(c10pep, function(x){c10gr[c10gr$pep==x,1]})
c10MnSD=data.frame(unlist(c10Mean), unlist(c10SD), unlist(c10vagr), stringsAsFactors = FALSE)
colnames(c10MnSD)=c("Mean","SD","Class")
c10clpar=aggregate(c10MnSD[,c(1,2,5)],list(unlist(c10MnSD[,3])),"mean")
clnm=c10clpar$Group.1
clnot=clnm[1:8]

# 
# Normalized NND
#

n=dim(c10vnn)[2]
c10clnndist=lapply(clnm, function(c){
                            x=as.matrix(dist(scale(c10vnn[c10MnSD[,3]==c,], center = TRUE)))
                            sapply(1:nrow(x),function(l){min(x[l,-l])})
                          })
c10clN=sapply(clnot,function(i){nrow(c10vnn[c10MnSD[,3]==i,])}) # N of points (peptides) per library
c10vol=prod(2*colSds(c10vnn))*2*pi^(n/2)/(n*gamma(n/2)) # volume of the cloud of all the data as n-dimensional elipsoid with radii aproximated as the 2SD along each dimension
c10clspecdist0=(c10vol/c10clN)^(1/n)  # library specific mean distance between the points as the side of the hypercube  of the volume per point
c10clnndistNorm=lapply(1:8,function(i){c10clnndist[[i]]/c10clspecdist0[i]}) # the normalized NND per library
names(c10clnndistNorm)=clnot
c10clnndistNormn=sapply(c10clnndistNorm, median)
c10clnndistNorm=lapply(order(c10clnndistNormn, decreasing = TRUE),function(i){c10clnndistNorm[i]})
c10clnndistNorm=unlist(c10clnndistNorm, recursive = FALSE)
c10clnndistNorm_b=c10clnndistNorm[names(order(unlist(lapply(c10clnndistNorm,mean))))]

mc10clnndist=melt(c10clnndistNorm)
mc10clnndist$L1=factor(mc10clnndist$L1)
c10nnglm=glm(log(value)~L1, data=mc10clnndist)
summary(glht(c10nnglm, mcp(L1="Tukey")))

boxplot(c10clnndistNorm, notch=TRUE, ylab="NND", las=2, log="y", varwidth=TRUE)
boxplot(lapply(c10clnndistNorm,log), notch=TRUE, ylab="log (NND)", las=2,  varwidth=TRUE)

# Total correlation between the peptides in each library 

c10clKLD0=sapply(clnot, function(c){
  medKLD(c10vnn[c10MnSD[,3]==c,])
})
colnames(c10clKLD0)=clnot

# Correlations between patients as z score

require(psychometric)
require(reshape2)
require(multicomp)

c10clcorpat=lapply(clnm,function(c){X=cor(c10vnn[c10MnSD[,3]==c,]); X[lower.tri(X)]})
c10clcorpatz=c10clcorpat
c10clcorpatz=lapply(c10clcorpatz,r2z)
c10clcorpatzs=c10clcorpatz[1:8]
corpatzmn=lapply(c10clcorpatzs, mean)
c10clcorpatzso=sapply(order(unlist(corpatzmn), decreasing=TRUE),function(i){c10clcorpatzs[i]})
corptzDf=as.data.frame(c10clcorpatzs)
corptzm=melt(corptzDf)
corpatsLm=glm(value~variable,data=corptzm)
summary(glht(corpatsLm, mcp(variable="Tukey")))
                    


# Constructs fig 5 - criteria

pdf(file="fig5.pdf", width=8, height=10)
par(mfrow=c(2,2))
vw=TRUE
boxplot(scale(Mean[Class %in% clnot], center = TRUE, scale = FALSE)~reorder(Class[Class %in% clnot],Mean[Class %in% clnot], mean), data=c10MnSD, las=2, notch=TRUE, ylab=" Mean Log Intensity",  varwidth=vw )
boxplot(as.data.frame(c10clKLD0[,order(colMeans(c10clKLD0), decreasing = TRUE)]),notch=TRUE,las=2, ylab="Total Correlation/ KLD [bits]",  varwidth=vw)
boxplot(c10clcorpatzso,notch=TRUE, ylab="z Score", las=2,  varwidth=vw)
boxplot(lapply(c10clnndistNorm_b,log), notch=TRUE, ylab="log(NND)", las=2,  varwidth=vw)
dev.off()
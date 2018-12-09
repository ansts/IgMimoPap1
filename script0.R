# uniqueXX_XX.txt files are output of the scripts of Matocko et al. 
# applied to the FASTQ files from the sequencing of phages
# 3 experiments with suffix 1,2 or 3 - 1 and 2 are repetition of the main experiment
# 3 - the same with a preamplified library.
# The data for forward and reverse reads is pooled and then also the data for experiments 1 and 2 # to yield P12 - sequenes and copies of reads ordered in reverse by copy number, 
# data file ODX0.rda contains the same data renamed as ODX0

P1_QF=read.table("uniqueP1_QF.txt")
P1_QR=read.table("uniqueP1_QR.txt")
P1=rbind(P1_QR,P1_QF)
P2_QF=read.table("uniqueP2_QF.txt")
P2_QR=read.table("uniqueP2_QR.txt")
P2_QR=read.table("uniqueP2_QP.txt")
P2_QR=read.table("uniqueP2_QР.txt")
P2=rbind(P2_QR,P2_QF)
P3_QR=read.table("uniqueP3_QR.txt")
P3_QF=read.table("uniqueP3_QF.txt")
P3=rbind(P3_QR,P3_QF)

P1[,1]=as.character(P1[,1])
p1com=intersect(P1_QF[,1],P1_QR[,1])
P1o=P1[P1[,1] %in% p1com,]
P1o=aggregate(P1o[,2], by=list(P1o[,1]), FUN = sum)
P1=P1[-which(P1[,1] %in% P1o[,1]),]
colnames(P1o)=colnames(P1)
P1=rbind(P1,P1o)

P2[,1]=as.character(P2[,1])
p2com=intersect(P2_QF[,1],P2_QR[,1])
P2o=P2[P2[,1] %in% p2com,]
P2o=aggregate(P2o[,2], by=list(P2o[,1]), FUN = sum)
P2=P2[-which(P2[,1] %in% P2o[,1]),]
colnames(P2o)=colnames(P2)
P2=rbind(P2,P2o)

P12=rbind(P1,P2)
p12com=intersect(P1[,1],P2[,1])
P12o=P12[P12[,1] %in% p12com,]
P12o=aggregate(P12o[,2], by=list(P12o[,1]), FUN = sum)
P12=P12[!P12[,1] %in% P12o[,1],]
colnames(P12o)=colnames(P12)
P12=rbind(P12,P12o)
P12=P12[order(P12[,2], decreasing=TRUE),]
seqerr=aamut(P12[P12[,2]>99,1])
P12nm=P12[!P12[,1] %in% seqerr,]

# ODXnm - P12 with no "false" sequences arising from few reading errors of highly represented sequences  
# Taking only sequences with 3-10 copies

ODXsel=ODXnm[ODXnm[,2]>2 & ODXnm[,2]<11,]

# The peptides in ODXsel are subjected to GibbsCluster off-line 
# ODXs_790_30 -> the optimal clustering in 790 clusters, read from the GibbsCluster output files

ODXs_790_30=read.table("ODXs_790_30.core", stringsAsFactors = FALSE)[,-2]

# Constructing PWM for each cluster from the clusterred sequences

odx79030pwmx=cnstrPWMtab0("ODXs_790_30.core")
odx79030Lpwmx=apply(odx79030pwmx,1,function(x){list(matrix(x, nrow=7, byrow=TRUE))})


for (i in 1:790) {colnames(odx79030Lpwmx[[i]][[1]])=aa}

# Selection of representative sequences

ODXs79030_sc=sapply(1:790, function(i){y=strsplit(ODXs_790_30[ODXs_790_30[,2]==(i-1),1], split=""); median(sapply(y,function(x){profchk(odx79030Lpwmx[[i]][[1]],x,7)}))})
ODXs79030pwmsc=lapply(1:790,function(i){list(odx79030Lpwmx[[i]][[1]],ODXs79030_sc[i])})
rndpODXs79030pwmpv=findpvscore(ODXs79030pwmsc, rndpepchr)
rndpODXs79030pwmEv=sapply(1:790,function(i){1-pbinom(length(ODXs_790_30[ODXs_790_30[,2]==i,1])/2,224000,rndpODXs79030pwmpv[[i]])})
chosen79030=rndpODXs79030pwmEv<1e-8
ixchosen=which(chosen79030)
ixsuperhi=which(rndpODXs79030pwmEv==0)
ODXs79030pepch=lapply(ixchosen,function(i){ODXs79030pep[[i]]})
ODXs79030scall=sapply(1:790, function(i){y=strsplit(ODXs_790_30[ODXs_790_30[,2]==(i-1),1], split="");sapply(y,function(x){profchk(odx79030Lpwmx[[i]][[1]],x,7)})})
ODXs79030scchosen=ODXs79030scall[chosen79030]
names(ODXs79030scchosen)=ixchosen
pep79030_99=lapply(1:length(ixchosen), function(i){ODXs79030pepch[[i]][ODXs79030scchosen[[i]]>quantile(ODXs79030scchosen[[i]], probs=0.9873)]})
peppostive=unlist(pep79030_99)
pepnegtve=unlist(pep79030_01)


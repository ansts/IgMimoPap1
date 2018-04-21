#sketch of the algorithm with the exact steps taken to reach the peptide sets

#generate auxiliary array containing the aa letters in the predetermined frequency known for Ph.D-7
rndpepchr=rndpeppre(2300000)

ODxs790_30.out=readGibbs_out("ODXs_790_30.out") # 30 seeds - gave the highest KLD score
ODXs_790_30=read.table("ODXs_790_30.core", stringsAsFactors = FALSE)[,-2]
odx79030pwmx=cnstrPWMtab0("ODXs_790_30.core")
odx79030Lpwmx=apply(odx79030pwmx,1,function(x){list(matrix(x, nrow=7, byrow=TRUE))})
for (i in 1:790) {colnames(odx79030Lpwmx[[i]][[1]])=aa}

# median score of each cluster using its own pwm
ODXs79030_sc=sapply(1:790, function(i){y=strsplit(ODXs_790_30[ODXs_790_30[,2]==(i-1),1], split=""); median(sapply(y,function(x){profchk(odx79030Lpwmx[[i]][[1]],x,7)}))})
# binding the pwms and the median scores in the same container
ODXs79030pwmsc=lapply(1:790,function(i){list(odx79030Lpwmx[[i]][[1]],ODXs79030_sc[i])})
# the peptides by cluster
ODXs79030pep=lapply(1:790, function(i){ODXs_790_30[ODXs_790_30[,2]==(i-1),1]})
# one container for the four custering schemes achived with 1,15,20 and 30 seeds
ODXs=list(ODXs79030=ODXs_790_30,ODXs79001=ODXs_790,ODXs79015=ODXs_790_15,ODXs79020=ODXs_790_20)
# finding consensus sequences groupd together by at least 3 out of 4 clustering schemes - the concordant clustering by "30" and "1" is required since they have the best score
ODXs_x=clustensus(ODXs)
ODXsmx=lapply(ODXs_x,function(x){consensusMatrix(x)[2]})
ODXsxpwm=lapply(ODXs_x,function(x){cnstrPWM(x)})

# find the p of random occurrence of the profile score >= threshold
rndpODXs79030pwmpv=findpvscore(ODXs79030pwmsc, rndpepchr) # currently done with 3*224000
# calculate the p of finding as many peptides (half the cluster since this > of the criterion was the median) of profile score  >= threshold in the N of tested pepts
rndpODXs79030pwmEv=sapply(1:790,function(i){1-pbinom(length(ODXs_790_30[ODXs_790_30[,2]==i,1])/2,224000,rndpODXs79030pwmpv[[i]])}) 
# chose 790_30 clusters with binomial p of occurrence <0.01 - the idea of consensus clusters was abandoned due to lower number of selected 
# clusters than using the 790_30 alone
chosen79030=rndpODXs79030pwmEv<1e-4
ixchosen=which(chosen79030)
#... and those with p=0
ixsuperhi=which(rndpODXs79030pwmEv==0)
# the peptides of the chosen clusters
ODXs79030pepch=lapply(ixchosen,function(i){ODXs79030pep[[i]]})
# calculate the scores of all peptides according to the pwm of the cluster they belong to
ODXs79030scall=sapply(1:790, function(i){y=strsplit(ODXs_790_30[ODXs_790_30[,2]==(i-1),1], split="");sapply(y,function(x){profchk(odx79030Lpwmx[[i]][[1]],x,7)})})
# ... and those of the chosen
ODXs79030scchosen=ODXs79030scall[chosen79030]
names(ODXs79030scchosen)=ixchosen
# from the peptides of the chosen clusters, select those of the lowest scores...
p=1e-6
pep79030_01=lapply(1:length(ixchosen), function(i){ODXs79030pepch[[i]][ODXs79030scchosen[[i]]<quantile(ODXs79030scchosen[[i]], probs=p)]})
# ... and of the high scores
pep79030_99=lapply(1:length(ixchosen), function(i){ODXs79030pepch[[i]][ODXs79030scchosen[[i]]>quantile(ODXs79030scchosen[[i]], probs=(1-p))]})
peppostive=unlist(pep79030_99)
pepnegtve=unlist(pep79030_01)


# the highest scoring cluster and a smaller one
# rank the cluster Nos by the binomial p value so the first are the most significant
rank79030pvals=names(ODXs79030scchosen)[order(rndpODXs79030pwmEv[chosen79030])]
ODXs79030pep[as.double(rank79030pvals[1:6])]
lengths(ODXs79030pep[as.double(rank79030pvals[1:6])])
head(rank79030pvals)
# upper half of score distribution of the 5 best clusters
pep_top_5=lapply(rank79030pvals[1:5], function(x){i=as.double(x); ODXs79030pep[[i]][ODXs79030scall[[i]]>ODXs79030_sc[[i]]]}) 
pep_top_5=unlist(pep_top_5)


# pepnegtlo - sample of 590 from peplo, which are the lowest score peptides from:
# peplotrash, which are ODXs790trash with Corrected score<10 and Bg score<10
# ODXs790trash are those of ODXs790out with Corrected score< =Bg score or Corrected sc<12
# pepnegtve is pep79030_01 


top_5=as.double(rank79030pvals[1:5])
top_5_pred=smpep(ODXs79030pwmsc[top_5],rndpepchr)
top_5pred=unlist(sapply(1:5, function(i){sample(sapply(1:lengths(top_5_pred)[i], function(j){paste(top_5_pred[[i]][[j]],sep="",collapse="")}), 150)}))
top_5pred=c(top_5pred)



ODXs79030trash=as.data.frame(t(sapply(ODxs790_30.out,function(x){x[x$Corrected==min(x$Corrected),]})))
peplotrash=ODXs79030trash$Alignment[ODXs79030trash$Self<5]
pepnegtlo=unlist(peplotrash)


length(intersect(peplo,pepnegtve)) # very few?!?!

tu_out=readGibbs_out("tucalls.out")
tu_c=read.table("tucalls.core", stringsAsFactors = FALSE)[,-2]
pepX=BString(pep)
pepsaaX=AAStringSet(paste(peps, collapse=""))
aaTu=alphabetFrequency(pepsaaX)
aaTu=aaTu/nchar(paste(peps, collapse=""))
aaTu=unlist(aaTu[1:20])
tu_pwmx=cnstrPWMtab0("tucalls.core", aas=aaTu)
tu_Lpwmx=apply(tu_pwmx,1,function(x){list(matrix(x, nrow=7, byrow=TRUE))})
for (i in 1:11) {colnames(tu_Lpwmx[[i]][[1]])=aa}
tu_sc=sapply(1:11, function(i){y=strsplit(tu_c[tu_c[,2]==(i-1),1], split=""); median(sapply(y,function(x){profchk(tu_Lpwmx[[i]][[1]],x,7)}))})
tupwmsc=lapply(1:11,function(i){list(tu_Lpwmx[[i]][[1]],tu_sc[i])})

rndpepTuchr=rndpeppre(400000, aa=aaTu)
rndp_tu_pwmpv=findpvscore(tupwmsc, rndpepTuchr)
rnd_tu_pwmEv=sapply(1:11,function(i){1-pbinom(length(tu_c[tu_c[,2]==i,1])/2,4524*9,rndp_tu_pwmpv[[i]])})
chosen_tu=rnd_tu_pwmEv<0.05
tixchosen=which(chosen_tu)
tu_pep=lapply(1:11, function(i){tu_c[tu_c[,2]==(i-1),1]})
tu_pepch=lapply(tixchosen,function(i){tu_pep[[i]]})
tu_scall=sapply(1:11, function(i){y=strsplit(tu_c[tu_c[,2]==(i-1),1], split="");sapply(y,function(x){profchk(tu_Lpwmx[[i]][[1]],x,7)})})
tuscchosen=tu_scall[chosen_tu]
names(tuscchosen)=tixchosen
pep_tu_80=lapply(1:length(tixchosen), function(i){tu_pepch[[i]][tuscchosen[[i]]>quantile(tuscchosen[[i]], probs=0.80)]})
pep_tu_80=unlist(pep_tu_80)

turndprofsmpl=smpep(tupwmsc, rndpepTuchr)
pepturnd=unlist(sapply(1:10, function(i){sapply(1:lengths(turndprofsmpl)[i], function(j){paste(turndprofsmpl[[i]][[j]],sep="",collapse="")})}))
 
tm=proc.time()
bot_rnd=negsmpep(ODXs79030pwmsc,rndpepchr, n=2000)
print(proc.time()-tm)
pepnegrnd=unlist(lapply(bot_rnd, function(x){paste(x, collapse="")}))
pep2000worst=pepnegrnd
pepnegrnd=head(pepnegrnd, 753)

peptotest[,1]=sapply(peptotest[,1], function(x){paste(x,"GGGS", collapse = "", sep = "")})


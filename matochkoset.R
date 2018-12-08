# From:
# W. L. Matochko, R. Derda. "Error analysis of deep sequencing of phage libraries: peptides censored in 
# sequencing". Comput Math Methods Med, 2013, 491612. 2013.
#
# http://www.chem.ualberta.ca/∼derda/mathbiology/PhD7-Amp-30F.txt
#

matochkoF30=read.table(file = file.choose(), stringsAsFactors = FALSE)
match7mers=matochkoF30[,2]
alpbmatch=unique(unlist(strsplit(as.character(match7mers),split="")))
alpbmatch
require(stringi)
match7astrxi=which(!is.na(stri_locate_first_fixed(match7mers, "*")[,1]))
match7mers=match7mers[-match7astrxi]
length(match7mers)
write.table(match7mers, file="match7mers.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

require(Biostrings)
consmx0=consensusMatrix(as.character(match7mers))/length(match7mers)
consmxmim=consensusMatrix(mims)/length(mims)
PSSMim=log2(consmxmim/consmx0)
PSSMim=t(PSSMim)
write.table(PSSMim, file="PSSmim.txt", row.names = FALSE, quote = FALSE)

mtch7merbound=matochkoF30[matochkoF30[,3]>2&matochkoF30[,3]<11,2]
consmx0b=consensusMatrix(as.character(mtch7merbound))/length(mtch7merbound)
consmx0b=consmx0b[-1,]
PSSMimb=log2(consmxmim/consmx0b)
PSSMimb=t(PSSMimb)
write.table(PSSMimb, file="PSSmimb.txt", row.names = FALSE, quote = FALSE)


mtch7merbound3=matochkoF30[matochkoF30[,3]>6&matochkoF30[,3]<61,2]
mtch7merbound3=mtch7merbound3[-grep("\\*",mtch7merbound3)]
consmx0b3=consensusMatrix(as.character(mtch7merbound3))/length(mtch7merbound3)
PSSMimb3=log2(consmxmim/consmx0b3)
PSSMimb3=t(PSSMimb3)
write.table(PSSMimb3, file="PSSmimb3.txt", row.names = FALSE, quote = FALSE)
write.table(mtch7merbound3, file="matochk07bound3.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

mtch7mer10up=matochkoF30[matochkoF30[,3]>10,2]
mtch7mer10up=mtch7mer10up[-grep("\\*",mtch7mer10up)]
write.table(mtch7mer10up, file="matochko7mer10up.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

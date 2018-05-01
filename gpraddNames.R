# Adds group names for the peptides in the Name column of .gpr
#
gpraddNames=function(fin, p=NULL, Clf="Classes.csv", fout=NULL){
  require(stringr)
  f=paste(p,"/",fin,sep = "")
  cl= read.csv(Clf, stringsAsFactors = FALSE, header = FALSE)
  fcon=file(f)
  fl=readLines(con=fcon, n=2)
  nln=as.double(str_extract(fl[2], "[1-9]+"))
  fhead=readLines(con = fcon, n=nln+2)
  fl=read.delim(f, skip=nln+2, header=T, check.names = F, stringsAsFactors = FALSE)
  close(fcon)
 
  if (nrow(fl)==length(cl[,1])) fl$`Name`=cl[,1] else fl$`Name`=cl[cl[,1]!="cntr",1]
  #rpl=length(fl$`ID`[fl$`ID`==0])
  #if (rpl>0) {
  #  rpv=runif(rpl,min=min(fl$`F532 Median`),max=max(fl$`F532 Median`))
  #  fl$`F532 Median`[fl$`ID`==0]=rpv
  #  fl$`F532 Mean`[fl$`ID`==0]=rpv
  #}
  p=paste(p,"_class", sep="")
  if (!file.exists(p)) dir.create(p)
  fout=paste(p,"/",fin,sep = "")
  fconw=file(fout, 'w')
  writeLines(fhead,con=fconw)
  write.table(fl,file=fconw, row.names=F, sep = '\t')
  close(fconw)
}

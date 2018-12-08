aamut=function(p){ 
  # Mutate every position of a peptide to evry other AA one position at a time,
  # return a list of the muatants - this is a new function 15.09.17. The original is mutagen
  unlist(p)
  aa=c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V") 
  paa=unlist(strsplit(p, split=""))
  res=lapply(seq_along(paa),function(i){
      naa=aa[aa!=paa[i]]
      lapply(naa,function(x){
        a=paa
        a[i]=x
        return(a)})
  })
  res=unlist(res, recursive = FALSE)
  res=lapply(res,function(x){paste(x,sep="", collapse="")})
  return(unlist(res))
}
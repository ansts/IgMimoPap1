#
# Generates LO PWM from an alignment seq, aas are the background frequencies
#
LOPWMtab=function(seq, aas=aaPh7, aan=aa){
  require(Biostrings)
  mx=array()
  rnms=list()
  pssm=consensusMatrix(seq)
  l=nchar(seq[1])
  pssm0=matrix(0,nrow=20,ncol=l)
  rownames(pssm0)=aan
  for (rn in rownames(pssm)){
    pssm0[rn,]=pssm[rn,]
  }
  aas=aas[aa[order(rownames(aas))],2]
  pwm=LOpwm(pssm0, aas)
  
  return(pwm)
}
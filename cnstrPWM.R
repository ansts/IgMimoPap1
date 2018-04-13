
# Construct PSSM and LO PWM tables from sequnces aligned  
# n is the alignment, aas - the background frequencies and aan - the amino acid code
#
cnstrPWM=function(n, aas=aaPh7_ok, aan=aa){
  require(Biostrings)
      aln=n
      L=length(aln)
      l=nchar(aln[1])
      pssm0=matrix(0,nrow=20,ncol=l)
      rownames(pssm0)=aan
      pssm=consensusMatrix(aln)
      for (rn in rownames(pssm)){
          pssm0[rn,]=pssm[rn,]
        }
      pwm=LOpwm(pssm0, aas)

  return(list(pssm0/L,pwm))
}
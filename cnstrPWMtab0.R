# Constructs PWM tables from a .core file generated from gibbscluster 
#- one for each cluster defined by the second column
# Each PWM is vectorized and represented by a single row in the final matrix
#
cnstrPWMtab0=function(n, aas=aaPh7_ok, aan=aa){

  require(Biostrings)
  mx=array()
  rnms=list()
  #aan=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

      fnm=n
      tbl=read.table(fnm, stringsAsFactors = FALSE)[,-2]
      l=nchar(tbl[1,1])
      for (z in 0:max(tbl[,2])){
        pssm0=matrix(0,nrow=20,ncol=l)
        rownames(pssm0)=aan
        aln=tbl[tbl[,2]==z,1]
        pssm=consensusMatrix(aln)
        for (rn in rownames(pssm)){
          pssm0[rn,]=pssm[rn,]
        }
        pwm=LOpwm(pssm0, aas)
        #if (i==1 & j==2 & z==1) {pwmout=pwm}
        rnms=c(rnms,paste(c(n,"_",z), collapse=""))
        pwm=as.vector(pwm)
        mx=cbind(mx,pwm)
    }
  
  mx=as.matrix(mx[,-1])
  cnms=paste(aan,rep(1:l,each=20))
  mx=t(mx)
  rownames(mx)=rnms
  colnames(mx)=cnms
  return(mx)
}
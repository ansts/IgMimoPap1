# log odds pwm - p is a PSSM, aa is background freqs 
LOpwm=function(p, aa){ 
  N0=sum(p[,1])
  am=matrix(rep(aa,ncol(p)), nrow=length(aa))
  m=log2(pscnt(p, N0, aa)/am)
  return(m)
}
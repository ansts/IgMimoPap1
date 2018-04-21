findallscore=function(pwmx=NULL, pepl) {
  #screens a set of peptides from GibbsCluster .core file with a set of PWM-ces
  #
  require(parallel)
  ppl=strsplit(pepl[,1],split="")
  cln=pepl[,2]
  l=nchar(pepl[1,1])
  cl <- makeCluster(getOption("cl.cores", 6))
  clusterExport(cl,varlist=c("profchk", "pwmx", "ppl", "cln"), envir = environment())
  # clusterEvalQ(cl, library())
  x=1:nrow(pepl)
  res=parSapply(cl,x, function(i) profchk(pwmx[[cln[[i]]+1]],ppl[[i]],l))
  stopCluster(cl)
  return(res)
}
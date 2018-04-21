findpvscore=function(pwmxsc,rpep) {
  #screens a set of (random) peptide,sstrsplit to character arrays, with a set of PWM-ces with respective threshold scores attached 
  #and for each PWM returns the proportion of the rnd pep with score equal or higher then the score attached 
  require(parallel)
  l=length(rpep[[1]])
  cl <- makeCluster(getOption("cl.cores", 7))
  clusterExport(cl,varlist=c("ODXs_sc"), envir = environment())
  clusterExport(cl,"profchk")
  # clusterEvalQ(cl, library())
  res=parLapply(cl,pwmxsc, function(x){1-ecdf(sapply(rpep,function(p){profchk(x[[1]],p,l)}))(x[[2]])})
  stopCluster(cl)
  return(res)
}
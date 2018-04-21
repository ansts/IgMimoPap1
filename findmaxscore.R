findmaxscore=function(pwmx, rpep) {
  #screens a set of peptides (strsplit to character arrays) with a set of PWM-ces
  #and for each PWM returns the max score found
  require(parallel)
  l=nchar(rpep[1])
  cl <- makeCluster(getOption("cl.cores", 6))
  clusterExport(cl,"profchk")
  # clusterEvalQ(cl, library())
  res=parLapply(cl,pwmx, function(x){quantile(sapply(rpep,function(p){profchk(x,p,l)}), c(0.95,0.99,0.999))})
  stopCluster(cl)
  return(res)
}
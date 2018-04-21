findmxscore=function(pwmx, rpep) {
  #screens a set of peptides (strsplit to character arrays) with a set of PWM-ces
  #and for each PWM returns the max score found - ver. predict
  require(parallel)
  l=nchar(rpep[1])
  #cl <- makeCluster(getOption("cl.cores", 6))
  #clusterExport(cl,"profchk")
  res=lapply(pwmx, function(x){max(unlist(lapply(rpep,function(p){profchk(x,p,n=nrow(x))})))}) 
  #stopCluster(cl)
  return(res)
}
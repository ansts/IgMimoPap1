rndpeppre=function(n, l=7,aa=aaPh7_ok) {
  # generates n peptides of length l with predetermined frequency aa as chr arrays
  require(parallel)
  require(Biostrings)
  aas=AA_ALPHABET[1:20]
  cl <- makeCluster(getOption("cl.cores", 8))
  clusterExport(cl,varlist=c("n", "l", "aa", "aas"), envir = environment())
  res=parLapply(cl,1:n, function(i){sample(aas,7,replace = TRUE, prob = aa)})
  stopCluster(cl)
  return(res)
}
profchk=function(pwm, p, beg=1, n, aan=aa){
  #Check a peptide for relevance to a PWM profile
  #Yields a single value score equal to the sum of the LO for each residue
  require(parallel)
  if (nrow(pwm)!=n) return("Incompaible PWM!")
  if ((beg+n-1)>length(p)) return("Window incompatible with peptide!")
  if (n>length(p)) return("")
  return(sum(sapply(1:n,function(i){pwm[i,p[i+beg-1]]})))
}
pscnt=function(n, N0, bg){  # pseudocaounts based probability
 psc=(n+bg*sqrt(N0))/(N0+sqrt(N0)) 
 return(psc)
}
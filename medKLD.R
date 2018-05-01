#
# Calculate Total Correlation between the rows of mx 
# in the form of KLD with independent distribution 
# (random combination of the column values - maximal entropy) 
# as the model distribution (Q). By changing b to p
# the opposite KLD comparison is calculated (P and Q are swapped)
# The theoretical distribution is calculated by discretizing 
# the intensity values in only to categories  - low and high. 
# Thus, for 10 variables (patients) 1024 categories 
# are formed which is still more than the size 
# of some of the libraries. The threshold for the discretizing
# is md.
#
medKLD=function(mx,md,b="q"){
    require(combinat)
    require(optimx)
    k=ncol(mx)
    n=nrow(mx)
    Qi=sapply(1:k, function(c){combn(k,c)})
    
    Marg=apply(mx,2,function(co){
      x=length(co[co>md])/length(co)
      y=1-x
      c(x,y)
    })
    
    bns=matrix(rep(1,k*2^k), ncol=k)
    c=1
    for (i in seq_along(Qi)){
      if (is.null(dim(Qi[[i]]))) {
          bns[c,Qi[[i]]]=2
          c=c+1
      }
      else{
      for (j in 1:ncol(Qi[[i]])){
        bns[c,Qi[[i]][,j]]=2
        c=c+1
        }
      }
    }
    Q=apply(bns,1,function(r){X=Marg[1,]
                            X[r==2]=Marg[2,r==2]
                            prod(X)})
    th=(mx>md)*1+1
    P=apply(bns,1,function(r){X=apply(th,1,function(rt){all(rt==r)})
                            length(X[X])/n
    })
    if (b=="q") KLD=sum(P[P!=0]*log2(P[P!=0]/Q[P!=0])) else KLD=sum(Q[P!=0]*log2(Q[P!=0]/P[P!=0]))
    return(c(KLD))
}

cte=function(ct,mx,bn){
  require(entropy)
  n=nrow(mx)
  N=length(mx)
  Marg=apply(mx,2,function(co){
    x=length(co[co>ct])/length(co)
    y=-1-x
    c(x,y)
  })
  
}
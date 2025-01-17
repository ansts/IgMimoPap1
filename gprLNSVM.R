#
# Reads fin gpr file with p(ath) and writes  the fout with B(ackground) columns equal 
# to the corrected intensities of the empty chip and F columns equal to the corrected 
# intensities of the stained chip. The local normalization is based on SVR filtering
# of the background intensity determined on the bases of the lower 1/p0 quantile 
# followed by SVR filtering of the standard deviation.
# Flags: wr(ite), sh(o)w 
# and R(ed/Green)sw(ap) - to make the green channel 
# under red column names for pepStat. 
#
gprLNSVM<-function(fin, p=NULL, fout=NULL, wr=T, shw=F, Rsw=T) {
  require(parallel)
  require(reshape2)
  require(stringr)
  require(e1071)
  prt=proc.time()
  finp=paste(p,fin,sep = "")
  fcon=file(finp)
  f2l=readLines(con=fcon, n=2)
  nlns=as.double(str_extract(f2l[2], "[1-9]+"))
  fhead=readLines(con = fcon, n=nlns+2)
  f=read.delim(finp, skip=nlns+2, header=T, check.names = F, stringsAsFactors = FALSE)
  close(fcon)
  if (max(f$`Block`)>1) f=f[f$`Block`==3,]
  f=f[f$`ID`!="000000000000000",]
  dfim=data.frame(cbind(f$`Row`, f$`Column`, f$`F532 Median`))
  colnames(dfim)=c("R","C","V")
  N=nrow(dfim)
  nr=max(dfim[,1])
  nc=max(dfim[,2])
  if (shw==T){
    img0=acast(dfim, R~C, value.var = "V")
    image(img0, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))
    lr=readline()
  }

  r3=nr%/%10
  c3=nc%/%10
  btX=matrix(runif(nr*nc, min(dfim[,3]), max(dfim[,3])),nrow=nr)
  btX=melt(t(btX))
  x=btX[,1]
  btX[,1]=btX[,2]
  btX[,2]=x
  colnames(btX)=c("R","C","V")
  btX[1:nrow(dfim),]=dfim
  btD=dfim[dfim[,1] %in% (1:r3),]
  btD[,1]=-btD[,1]
  btU=dfim[dfim[,1] %in% ((nr-r3):nr),]
  btU[,1]=2*nr-btU[,1]+1

  btX=rbind(dfim,btD,btU)
  
  btL=btX[btX[,2] %in% (1:c3),]
  btL[,2]=-btL[,2]
  btR=btX[btX[,2] %in% ((nc-c3):nc),]
  btR[,2]=2*nc-btR[,2]+1
  
  btX=rbind(btX, btL, btR)
  nrX=max(btX[,1])
  ncX=max(btX[,2])
  NX=nrow(btX)
  
  co=FALSE 
  p0=3
  for (i in 1:NX){ 
      #print((i*100)%/%(N*2))
      rs=min(2,btX[i,1])
      cs=min(2,btX[i,2])
      p1 =min(p0,rs+cs)
      rb=max(btX[i,1]-rs,min(btX[,1]))
      rf=min(rb+2*rs+1, nr)
      cb=max(btX[i,2]-cs,min(btX[,2]))
      cf=min(cb+2*cs+1, nc)
      z0=btX[btX[,1] %in% rb:rf&btX[,2] %in% cb:cf,]
      lr=rank(z0[,3])
      z1=(z0[lr<length(lr)%/%p1,])                         # For each spot take the spots in a patch (rs.2+1)Rx(cs.2+1)C around it
      if (co==FALSE) {                                                            # and compare the spot value to the 100/p th percentile of the spots is the patch.
          btm=z1                                        # If less add this spot the btm set
          co=TRUE
      } 
      else {
      # if (dfim[i,3]<z1){                                                            
      #     btm=rbind(btm,dfim[i,])                 
      #}
          btm=rbind(btm,z1)
      }
  }
  print("bottom ready")
  frm=dfim[dfim[,1]==1|dfim[,1]==nr|dfim[,2]==1|dfim[,2]==nc,]
  btm=rbind(btm,frm)
  btm=unique(btm)

  
  lmd=svm(btm[,1:2], btm[,3], cost = 1000, gamma=3,epsilon = .001)     
  prd=predict(lmd,dfim[,1:2])
  lnew=dfim
  lnew[,3]=as.double(dfim[,3]-prd)
  print("bottom fit")
  co=FALSE 

  for (i in 1:N){ 
    #print(((i+N)*100)%/%(N*2))
    rs=2
    cs=2
    rb=max(lnew[i,1]-rs,min(lnew[,1]))
    rf=min(rb+2*rs+1, nr)
    cb=max(lnew[i,2]-cs,min(lnew[,2]))
    cf=min(cb+2*cs+1, nc)
    z0=lnew[lnew[,1] %in% rb:rf&lnew[,2] %in% cb:cf,]
    sd0=c(lnew[i,1],lnew[i,2],sd(z0[,3]))
                 
    if (co==FALSE) {
      sdmx=sd0
      co=TRUE
    } 
    else {
      sdmx=rbind(sdmx,sd0)
    }
  }
  print("sd ready")
  colnames(sdmx)=c("R","C","V")
  rownames(sdmx)=rownames(lnew)
  sdmx=as.data.frame(sdmx)
  lmsd=svm(sdmx[,1:2], sdmx[,3], cost = 1000, gamma=3,epsilon = .001)
  prsd=predict(lmsd,lnew[,1:2])
  mprsd=mean(sdmx[,3])
  lnew[,3]=lnew[,3]/prsd*mprsd

  lv=min(lnew[,3])-1                                               
#  lnew[,3]=lnew[,3]-lv
#  pred=dfim[,3]-lnew[,3]
  dfim[,3]=lnew[,3]+prd
  print(c("0 and less - ",length(dfim[dfim[,3]<=0,3])))
  pred=prd
  if (shw==T){
    imgX=acast(btm, R~C, value.var = "V")
    image(imgX, main=c(fin," X bottom"), zlim=c(0,65000), col=topo.colors(128))
    lr=readline()
  }
 if (shw==T){
    imgX=acast(sdmx, R~C, value.var = "V")
    image(imgX, main=c(fin," X SD"), zlim=c(0,6000), col=topo.colors(128))
    lr=readline()
  }
  if (shw==T){
#    for (i in 0:3) {
#      x=c((1+i*(N%/%4)),((1+i)*(N%/%4)))
#      plot(dfim[x[1]:x[2],3], cex=.3, ylim = c(0, 65000))
#      par(new=T)
#      plot(pred[x[1]:x[2]], cex=.1,col=2, ylim = c(0, 65000))
#      lines(pred[x[1]:x[2]], cex=.3,col=2, ylim = c(0, 65000))
#      par(new=F)
#      lr=readline()
#      plot(log2(lnew[x[1]:x[2],3]), cex=.4)
#      lr=readline()
#      plot(log2(dfim[x[1]:x[2],3]), cex=.4)
#      lr=readline()
#      plot(dfimnew[x[1]:x[2]], cex=.4)
#      lr=readline()
#      }
    hist(lnew[,3],breaks=300)
    lr=readline()
    
    bkg=matrix(pred, dim(img0), byrow = T)
    image(bkg, main=c(fin," Background"), zlim=c(0,65000), col=topo.colors(128))
    lr=readline()
 
    img1=matrix(lnew[,3]-lv, dim(img0), byrow = T)
    image(img1, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))
    
    #img1[img1<0]=0
    #image(img1, main=c(fin," Filtered"), zlim=c(0,65000), col=topo.colors(128))
    print(max(img1))
    lr=readline()
  }
  else print(fin)
  if (wr) {
    f$`B532`=pred
    f$`B532 Median`=pred
    f$`B532 Mean`=pred

    #f$`F532`=dfimnew
    f$`F532 Median`=dfim[,3]
    f$`F532 Mean`=dfim[,3]

    if (Rsw==T) {
      fcn=colnames(f)
      fcn=str_replace_all(fcn, "635","630")
      fcn=str_replace_all(fcn, "532", "635")
      colnames(f)=fcn
    }
    newp=paste("proc_",p,sep = "")                       # Write the new .grp files in 
    dir.create(newp)                                        # subfolder ...../.
    fout=paste(newp,fin,sep = "")
    fconw=file(fout, 'w')
    writeLines(fhead,con=fconw)
    write.table(f,file=fconw, row.names=F, sep = '\t')
    close(fconw)
    print(proc.time()-prt)
  }
}

TruncKin=function(NDL,SelCFF,Ofset=-1,mask=NULL) #Extract subset of data from NormDat.txt files
  # NDL: DataSet
  # SelCFF: selected CFF
  # Ofset: Number of timesteps to suppress at start of each OMMs (10 to suppress initial state) -1 fro no ofset
  #mask: Column of data to extract
  #return :  list of Tabf: filtered file and SIDX: index of start for each OMMM in the file
{
  SIDX=c(1)
  Tab=read.table(NDL)
  Tab=Tab[Tab[,"CFF_Id"]==SelCFF,]
  Idx=array(0,dim=c(dim(Tab)[1],1))
  Tab=cbind(Idx,Tab)
  a=Tab[1,2]
  cte=0
  n=0
  Inited=FALSE
  idx=1
  for(i in 1:dim(Tab)[1])
  {
    if(Tab[i,2]!=a)
    {
      n=0
      a=Tab[i,2]
      idx=idx+1
      SIDX=c(SIDX,i)
    } else n=n+1
    if(n<=Ofset) next
    Tab[i,1]=idx
    if(Inited)Tabf=rbind(Tabf,Tab[i,])
    else{
      Tabf=Tab[i,]
      Inited=TRUE
    }
    cte=cte+1
  }
  colnames(Tabf)=colnames(Tab)
  Tabf=Tabf[1:cte,]
  if(!is.null(mask)) Tabf=Tabf[,mask]
  return(list(Tabf=Tabf,SIDX=SIDX))
 
}

StatKin=function(NDL,Parm1,Parm2,SelCFF,Ofset=0,mask=NULL) #utilities functions for NormDataPlot
  #Parm1,Parm2: parameters to select
  # SelCFF: selected CFF
  # Ofset: Number of timesteps to suppress at start of each OMMs (10 to suppress initial state) -1 fro no ofset
  #mask: Column of data to extract
{
  TB=TruncKin(NDL,SelCFF,Ofset,mask)
  Tab1=TB$Tabf
  TLS=TB$SIDX
  x=Tab1[,Parm1]
  y=Tab1[,Parm2]
  mx=median(x)
  my=median(y)
  Mt=cor(Tab1,method="pearson")
  test=cor.mtest(Tab1, method="pearson",conf.level = 0.95)$p
  testL=test[test==0]=1e-50
  dd=dim(test)[1]
  re=-array(sapply(test,log10),dim=c(dd,dd))
  cp=corrplot(Mt, p.mat = test, sig.level = 0.01,type = 'lower',number.cex=0.5,insig = "blank",col = COL2('PiYG')) # p-value 0.01 for 95% confidence
  text(cp$corrPos$x, cp$corrPos$y, round(cp$corrPos$corr, 2),cex=0.5,font=2,col="black")
  return (list(medParm1=mx,medParm2=my,Cortast=test,Parm1=Tab1[,Parm1],Parm2=Tab1[,Parm2],TLS=TLS))
}

valf= function(RC, D2, D4, A2)# geometric function
{
  A2=A2*pi/180
  fu=sqrt(RC*(RC-2*D4*cos(A2))+D4^2+D2^2)
  return(sqrt(RC*(RC-2*D4*cos(A2))+D4^2+D2^2))
} 

NormDataPlot=function(NDL,Parm1,Parm2,Ofset,SelCFF,mask, line=TRUE)
  #evaluate direct plots and correlations between mask items
{
  RC=3.233081 # CFF cen to C2 distance
  #NDL=DataSet #dataset
  rr=StatKin(NDL,Parm1,Parm2,SelCFF,Ofset,mask) #correlation plot
  cat("Clic on graph to continue")
  locator(1)
  Border=rr$TLS
  col=c()
  ncol=rainbow(length(Border))
  par(mfrow=c(1,1),mar = c(5, 5, 1, 1))
  for(i in 1:(length(Border)-1))
  {
    for(j in Border[i]:Border[i+1]) col[j]=ncol[i]
  }
  plot(rr$Parm1,rr$Parm2, col=col, xlab=paste(Parm1,"( )"),ylab=paste(Parm2,"( )"),cex.axis=2,cex.lab=2) # bipartite parm plot
  if (line) lines(lowess(rr$Parm1,rr$Parm2,f=0.1), col="black",lwd=2)
}

PlotCloud=function()
{
  mask=c("A0","A1","A3","A2","D2","D4",'IDENA','IDENB','dmin','imin','imin1')
  mask=c("A0","A1","A3","A2","D2","D4",'dmin','imin','imin1')
  Ofset=25#Ofset=25  Skipped time steps for ech OMM (-1: none)
  SelCFF=1 #SelCFF=2  selected CFF (1:effector,2: producer)
  NDL='MixLN/NormData.txt' #'MixLL/NormData.txt'
  NormDataPlot(NDL,"D4","imin",Ofset,SelCFF,mask, FALSE)
}

t_col= function(color, percent = 50)
{
  rgb.val=col2rgb(color)
  acol=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,
           alpha = (100 - percent) * 255 / 100)
  return(acol)
}

MDSProj=function(NDL,Parm1,Parm2,Ofset,SelCFF,mask) #Projection MDS reduced coordinate
{
  library(FactoMineR)
  library(ggplot2)
  library (dplyr)
  library(plotrix)
  mask=c("A0","A1","A3","A2","D2","D4",'IDENA','IDENB','dmin','imin','imin1')
  Ofset=25 #Ofset=25  Skipped time steps for ech OMM (-1: none)
  SelCFF=2 #SelCFF=2  selected CFF (1:effector,2: producer)
  NDL='MixLL/NormData.txt' #'MixLL/NormData.txt' #'MixLN/NormData.txt'
  
  mask1=c("A0","A1","A2", "A3","D2","D4")
  mask1=c("A1","A3","D2")
  mask1=c("A0","D2","A3")
  
  TB=TruncKin(NDL,SelCFF,Ofset,mask)
  #TB1=TruncKin('MixLN/NormData.txt',1,Ofset,mask)
  Tab1=TB$Tabf[,mask1]
  TLS=TB$SIDX
  f=scale(Tab1)
  g=dist(f)
  Cmd=cmdscale(g, k = 2)
  write.table(Cmd, file="temp.txt") 
  
  Cmd=read.table("temp.txt")
  #rb=rainbow(4) 
  rb=c("blue","green","red","black","purple","cyan","magenta","orange","brown","grey","darkgreen","darkblue","darkred","darkcyan","darkmagenta","darkorange","darkgrey","darkviolet","darkturquoise")
  rb=sapply(rb,t_col,0)
  crit=TB$Tabf$imin
  dmin=TB$Tabf$dmin
  col=rb[crit]
  Frame=array(TRUE,dim=dim(Cmd)[1])
  sel=c(1:9) #to set
  ThDmin=2.5 #to set
  ThDmax=10 #to set
  
  titre="MDS(A0, A1, A2, A3, D2, D4) dmin= 1 - 15 Å range"
  sel1=crit%in%sel & dmin<=ThDmax & dmin>ThDmin & Frame
  plot(Cmd[sel1,1],Cmd[sel1,2],col=col[sel1],xlim=c(-5,5),ylim=c(-3,4),xlab="AX1",ylab="AX2",main=titre,cex.lab=1.5,cex.axis=1.5)
  legend(4.3,4,legend=as.character(c(1:9)),col=rb[1:9],pch=20)
  
  range=11164:dim(Cmd)[1]
  Frame=array(FALSE,dim=dim(Cmd)[1])
  Frame[range]=TRUE
  col[1:dim(Cmd)[1]]="green"
  col[range]=rgb(1,0,0,0.10)#"red" #
  sel1=crit%in%sel & dmin<=ThDmax & dmin>ThDmin & Frame
  plot(Cmd[sel1,1],Cmd[sel1,2],col=col[sel1],xlim=c(-5,5),ylim=c(-4,3),xlab="AX1",ylab="AX2",main="MDS(A0, A1, A3, D2, D4)",cex.lab=2,cex.axis=2)
  
  write.table(Cmd, file="Hybrid_LN(1-11164)LL(11164-end).txt") 
  #Cmd=read.table("temp.txt")    to reload
}
  


# additional bloc of codes for reuse (do not run)

ToCircleClusters=function()
{
  NClust=3
  for(i in sel )
  {
    sel=i
    sel1=crit%in%sel & dmin<=ThDmax & dmin>ThDmin
    coord=as.data.frame(Cmd[sel1,])
    a=kmeans(coord,centers=NClust)
    Ncount=table(a$cluster)
    SdClust=array(0,dim=NClust)
    for(j in 1:NClust)
    {
      diff1=coord[a$cluster==j,1]-a$centers[j,1]
      diff2=coord[a$cluster==j,2]-a$centers[j,2]
      diff=cbind(diff1,diff2)
      SdClust[j]=median((apply(diff,1,norm_vec)))
      symbols(a$centers[j,1],a$centers[j,2],circles=SdClust[j], inches=F, add=T)
    }
  }
}

BuildBloc1=function() # clusterisation
{
library(FactoMineR)
library(ggplot2)
library (dplyr)
coord=as.data.frame(Cmd[sel1,])
a=kmeans(coord,centers=1)
table(a$cluster)
names (coord)= paste ("Axe", 1:2, sep = "")
coord= mutate (coord, classe = as.factor (a$cluster))
ggplot (coord, aes (x = Axe1, y = Axe2, color = classe)) + geom_point()
}

BuildBloc2=function()
{
  
  CorrAng=c(254.58,0,143.06,93.12) #angle to C2 of C1,c2,C3 C8 imin1-4
  CorrDist=c(3.23,3.23,3.23,2.65) #Dist cen to Cx
  val=expression(sqrt(RC*(RC-2*D4*cos(A2))+D4^2+D2^2))
  TT=TruncKin(NDL,SelCFF,Ofset,mask)
  Tab=TT$Tabf
  SIDX=TT$SIDX
  RR=array(0,dim=length(Tab$A2))
  AA=array(0,dim=length(Tab$A2))
  for(i in 1:length(Tab$A2))
  {
    AA[i]=Tab$A2[i]+CorrAng[Tab$imin1[i]]
    RR[i]= CorrDist[Tab$imin1[i]]
  }
  
  q=valf(RR,Tab$D2,Tab$D4,AA)
  q1=valf(RC,Tab$D2,Tab$D4,Tab$A2)
  m=max(Tab$A2[Tab$imin1==2])
  sel=Tab$imin1==1
  sel=TRUE
  qf=q1[sel]
  vf=Tab$A2[sel]
  plot(qf,vf)
  plot(Tab$Idx[SIDX],Tab$A2[SIDX])
  x=Tab$Idx[SIDX]
  y=Tab$A2[SIDX]
  y[9:13]=360-y[9:13]
  plot(30*(x-1),y-y[1])
  
  plot(q,Tab1$dmin)
  plot(Tab1$imin1,Tab1$imin)
  dD2=D(val,'D2')
  dD4=D(val,'D4')
  dA2=D(val,'A2')
  var_dmin=var(Tab1$dmin)
  var_d4=var(Tab1$D4)
  var_d2=var(Tab1$D2)
  var_a2=var(Tab1$A2)
  m_d2=mean(Tab1$D2)
  m_d4=mean(Tab1$D4)
  m_a2=mean(Tab1$A2)
  val_dD2_vec <- eval(dD2, data.frame(RC = RC, D4 = Tab1$D4, D2 = Tab1$D2,A2 =Tab1$A2))
  val_dD4_vec <- eval(dD4, data.frame(RC = RC, D4 = Tab1$D4, D2 = Tab1$D2,A2 =Tab1$A2))
  val_dA2_vec <- eval(dA2, data.frame(RC = RC, D4 = Tab1$D4, D2 = Tab1$D2,A2 =Tab1$A2))
  dif_D2 <- Tab1$D2-m_d2 *as.vector(val_dD2_vec)
  dif_D4 <- Tab1$D4-m_d4* as.vector(val_dD2_vec)
  dif_A2 <- Tab1$A2-m_a2 *as.vector(val_dA2_vec)
  V_d2=sum(dif_D2*dif_D2)/length(dif_D2)
  V_d4=sum(dif_D4*dif_D4)/length(dif_D4)
  V_a2=sum(dif_A2*dif_A2)/length(dif_A2)
  s=V_d2+V_d4+V_a2
  V_d2=V_d2/s
  V_d4=V_d4/s
  V_a2=V_a2/s
  plot(q-Tab1$dmin)
}
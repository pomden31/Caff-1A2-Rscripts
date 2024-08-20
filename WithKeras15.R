library(keras)
library(tensorflow)
#use_condaenv("r-reticulate")

#"MixLL/" c(1:12,62,63)
#"MixLN/" c(15,17:20,22,25:26,30:32,49:55,57:60,70,75)
#"MixRN/" c(16,24,33:41,43:48,61,71:74)
#"MixRL/" c(13,14)
#"MixRR/" c(64:69)
#"Serie669/" c(100:123)
#"Serie671/" c(76:99) 
#"Serie669_1CFF/" c(136:147)
#"Serie671_1CFF/" c(124:135)


Test=function()
{
#DataSet="MixLN/NormData.txt"
#DataSet=c("Serie671_1CFF")
#DataSet=c("MixRN_OXY")
#DataSet=c("MixLN_OXY","MixRN_OXY")
#DataSet=c("Serie671_1CFF","Serie671")
DataSet=c("Serie671")
#DataSet=c("Serie669_1CFF","Serie669")
#DataSet=c("Serie669_1CFF")
#DataSet=c("Serie669")
#DataSet=c("Serie669_1CFF","Serie669","Serie671_1CFF","Serie671","MixLL/","MixLN/")
#DataSet=c("Serie669_1CFF","Serie669","Serie669","Serie671_1CFF","Serie671","Serie671","MixLL/","MixLN/")
#DataSet=c("Serie669_1CFF","Serie671_1CFF","MixLN")
#DataSet=c("MixRN")
#DataSet=c("MixLN","Serie669_1CFF","Serie671_1CFF")
#DataSet=c("Serie669_1CFF","Serie669","Serie669","Serie671_1CFF","Serie671","Serie671","MixLL/","MixLL/","MixLN/")
#c(1,1,2,2,1,2,1,2,1) 
#DataSet=c("Serie669","Serie669","Serie671","Serie671","MixLL","MixLL","MixLN")
#DataSet=c("Serie669","Serie671")
#DataSet=c("MixRL","MixRL","MixRN","MixRR","MixRR","MixLL","MixLL","MixLN")
#DataSet=c("MixRN","MixLN")
#Mask=c(1,2,3,5,6) Mask=c(4,5,6) #Obs=4
# Case settings #SelCFF=c(1,2,2,2,2,1) c(1,2,1,2,2,1,2,1) c(1,2,1) SelCFF=c(1,2,1,2,1,2,1) 

#Global settings 
DataSet=c("MixRN","MixLN")
SetNames=c("RN000","LN000")
SelCFF=c(1,1) #SelCFF=c(1,1) # selected CFF (1:effector, 2: producer) 
Ofset=25 # Skipped frames at the start of each  OMM (-1: none) to allow relaxation
Type=c(2,2,2,2,1,1,0,0) # type of column 2 for angle 1 for distance 0 for codes (not to be normalized
Mask=c(1,2,4,5,6) #parameter used for learning and prediction (X)
Obs=3 #parameter used for learning and prediction (Y
id=c('A0','A1','A2','A3','D2','D4','dmin','imin')
Sampling=c(0.7,0.3) # learning and test ratios 
FilterOnIronState=TRUE
FilterTarget="OXY" 

RandomizeOrder=TRUE  #FALSE# TRUE: randomize order of data
RandomizeEnv=FALSE # TRUE: randomize environment parameters
RandomizeResp=FALSE  # TRUE: randomize response=f(environment)
PonderatePools=TRUE  # TRUE to ponderate pools weight by inverse of size during learning
ColorOnBatch=FALSE # TRUE to color learning and test set on batch

Batch_Id=c(2) # c(2,3)#NULL #c(1) # ID of batches to include in pool 1
SelectPoo0=1#1 # 0 : All batch valid,  1: as defined by Batch_Id, 2: all but defined by Batch_Id
LearnTestFilter=1 # 0: no Selection, 1 Select learning set, 2: select Test set ,3 Select both 

SortTSon=NULL #A3"#"dmin"#"A0" #sorting order f: see id above for possible values or time series NULL=no sorting
action="Learn"# Learn Save Load keras model

#-----------------------------------------------------------Load dataset
Tabx=GetNewTab(DataSet,TRUE,Ofset,SelCFF)
#-----------------------------------------------------------reload dataset (when already loaded)
Tab0=Tabx[[1]]
Tax=Tabx[[2]]
limits=Tabx[[3]]

#plot(Tab0$A3,Tab0$D4,col=rgb(0,0,0,0.3))
#plot(Tab0$A1,Tab0$A2,col=rgb(0,0,0,0.3))
#plot(Tab0$A3,Tab0$A2,col=rgb(0,0,0,0.3))
#plot(Tab0$D2,Tab0$A2,col=rgb(0,0,0,0.3))
#plot(Tab0$D4,Tab0$A2,col=rgb(0,0,0,0.3))
#Rand = sample(nrow(Tab0))
#Tab0[,"A2"]=Tab0[Rand,"A2"]
#Tab0[,"D4"]=Tab0[,"D4"]+rnorm(nrow(Tab0),0,1)
#Tab0[,"A2"]=60

sdn=c('E','P')
nname=paste(DataSet,sdn[SelCFF],sep="_")
ofl=0
u=1
Shift=c()
limits=c(limits,1)
ulimit=c(1)
for(i in 2:length(limits)) 
{
  if(limits[i]==1)
  {
    ofl=ofl+limits[i-1]-1
    Shift=c(Shift,ofl)
    next
  }
  ulimit=c(ulimit,limits[i]+ofl)
}
limits=ulimit
Shift=c(1,Shift)
valShift=rep(0,dim(Tab0)[1])
for(i in 1:(length(Shift)-1)) valShift[Shift[i]:(Shift[i+1]-1)]=i
valShift[dim(Tab0)[1]]=valShift[dim(Tab0)[1]-1]


sh=!is.null(Batch_Id)
MaskL=rep(TRUE,dim(Tab0)[1])

if(SelectPool!=0 && sh) 
{
  if(SelectPool==1)
  {
    MaskL[1:length(MaskL)]=FALSE
    for(i in 1:(length(Batch_Id))) MaskL[Shift[Batch_Id[i]]:(Shift[Batch_Id[i]+1]-1)]=TRUE
  }
  if(SelectPool==2)
  {
    MaskL[1:length(MaskL)]=TRUE
    for(i in 1:length(Batch_Id)) MaskL[Shift[Batch_Id[i]]:(Shift[Batch_Id[i]+1]-1)]=FALSE
  }
}




F1=which(Tax[,3]=="OXY") # allow separate ordering for OXY et FE OMMs
F2=which(Tax[,3]=="FE") # allow separate ordering for OXY et FE OMMs
TB1=Tab0[F1,]
TB2=Tab0[F2,]

#optional sorting for time series plot
if(!is.null(SortTSon))
{
  Ord=order(Tab0[ ,SortTSon])
  Tab0=Tab0[Ord,]
  Tax=Tax[Ord,]
  MaskL=MaskL[Ord]
  valShift=valShift[Ord]
}


#end sorting

#--------------------------------------------------

if(RandomizeOrder)
{
  Rand = sample(nrow(Tab0))
  Tab0=Tab0[Rand,]
  Tax=Tax[Rand,]
  MaskL=MaskL[Rand]
  valShift=valShift[Rand]
}
FiltFE=Tax[,3]==FilterTarget
if(FilterOnIronState) MaskL=MaskL & Tax[,3]==FilterTarget

if(FALSE)#optional plots
{
A0=Tab0[MaskL,"A0"]
A1=Tab0[MaskL,"A1"]
A2=Tab0[MaskL,"A2"]
A3=Tab0[MaskL,"A3"]
D2=Tab0[MaskL,"D2"]
D4=Tab0[MaskL,"D4"]
DM=Tab0[MaskL,"dmin"]
col=rep("black",length(x1))
col[DM<Thd1]=rgb(1,0,0,0.5)
col[DM<Thd2 & x7>=Thd1]=rgb(0,1,0,0.5)
col[DM>=Thd2 & x7<Thd3]=rgb(0,0,1,0.5)
col[DM>=Thd3 ]=rgb(0.5,0.5,0.5,0.5)
plot(A2,A0,col=col,xlim=c(-180,180),ylim=c(-100,100))
plot(A2,A1,col=col,xlim=c(-180,180),ylim=c(50,180))
plot(A2,A3,col=col,xlim=c(-180,180),ylim=c(-100,100))
plot(A2,D2,col=col,xlim=c(-180,180),ylim=c(-10,10))
plot(A2,D4,col=col,xlim=c(-180,180),ylim=c(0,15))
plot(D2,D4,col=col,xlim=c(-10,10),ylim=c(0,15))
plot(D2,D4)
}



res=GetData(Tab0, Type,Mask,Obs) 
res1=res
n=dim(res[[1]])[1]
Env=as.matrix(res[[1]])
Resp=as.matrix(res[[2]])
FullDat=as.matrix(res[[3]])

#f=scale(Env)
#g=dist(f)
#Cmd=cmdscale(g, k = 1)

if(RandomizeResp)
{
  Rand = sample(nrow(Resp))
  Resp=Resp[Rand,]
}
if(RandomizeEnv)
{
  Rand = sample(nrow(Env))
  Env=Env[Rand,]
  for(i in 1:ncol(Env))
  {
    Rand = sample(nrow(Env))
    #Env[,i]=Env[Rand,i]
  }
}

PoolSz=rep(0,length(Shift)-1)
for(i in 1:length(Shift)-1)
{
  PoolSz[i]=sum(valShift==i)
}
Szm=max(PoolSz)
PoolWeight=Szm/PoolSz
GSampleWeight=rep(1,length(valShift))
if(PonderatePools)
{
GSampleWeight=PoolWeight[valShift]
GSampleWeight[length(valShift)]=1
}

ss=sample(x=2, size=dim(Env)[1],replace=TRUE, prob=Sampling)
ss1=(ss==1) #no filter
ss2=(ss==2)
if(LearnTestFilter==1)
{
  ss1=(ss==1 & MaskL)  # train 
  ss2=(ss==2)  # test
}
if(LearnTestFilter==2)
{
  ss1=(ss==1 )  # train 
  ss2=(ss==2 & MaskL) # select for filtered test set
}
if(LearnTestFilter==3) 
{
  ss2=(ss==2  & MaskL)
  ss1=(ss==1  & MaskL) 
}



Nlearn=sum(ss1)
Ntest=sum(ss2) 
X_train=Env[ss1,]
X_test=Env[ss2,]
Y_train=Resp[ss1,]
SampleWeight=GSampleWeight[ss1]
TestWeight=GSampleWeight[ss2]
Y_test=Resp[ss2,]
X_full=Env
Y_full=Resp


custom_loss=function(y_true, y_pred, weights=1)
{
  mse_loss1 =((y_true[,1] - y_pred[,1])^2+(y_true[,2] - y_pred[,2])^2)
  s1=sqrt(y_true[,1]^2+y_true[,2]^2)+0.01
  s2=sqrt(y_pred[,1]^2+y_pred[,2]^2)+0.01
  mse_loss2=0
  mse_loss2 =(acos(y_pred[,1]/s2)-acos(y_true[,1]/s1))^2+(asin(y_pred[,2]/s2)-asin(y_true[,2]/s1))^2
  mse_loss=mean(weights*( mse_loss1 + 0*mse_loss2))
  return (mse_loss) 
}

action='Learn2'
if(action=='Learn1')
{
  # Refference Modèle Keras
  model <- keras_model_sequential() %>%
    layer_dense(units = 288, activation = "relu",input_shape = ncol(X_train))%>% 
    layer_dense(units = 72, activation = "relu")%>%
    layer_dropout(rate = 0.1) %>% 
    layer_dense(units = 24, activation = "relu") %>% 
    layer_dropout(rate = 0.1) %>% 
    layer_dense(units = 2)
  model %>% compile(optimizer ='Adam' , loss ='mse',weighted_metrics='mse') #mse
  
  history <- model %>% fit(X_train, Y_train, epochs = 75,batch_size=32,sample_weight=SampleWeight,validation_data = list(X_test, Y_test,TestWeight))
  model %>% evaluate(X_test, Y_test)
}

if(action=='Learn3')
{
  # Modèle Keras
  model <- keras_model_sequential() %>%
    layer_dense(units = 9, activation = "relu",input_shape = ncol(X_train)) %>% #288 #128
    layer_dense(units = 128, activation = "relu")%>%
    layer_dropout(rate = 0.15) %>% 
    layer_dense(units = 64, activation = "relu")%>%   #72 #64
    layer_dense(units = 12, activation = "relu") %>% #18 #16
    layer_dense(units = 2)
  
  model %>% compile(optimizer ='Adam' , loss = custom_loss,weighted_metrics=custom_loss) #"mse" "adam"
  history <- model %>% fit(X_train, Y_train, epochs = 75,batch_size=32,sample_weight=SampleWeight, validation_data = list(X_test, Y_test,TestWeight))
  model %>% evaluate(X_test, Y_test)
  Y_pred <- model %>% predict(X_test)
  Y_predL <- model %>% predict(X_train)
  Y_pred_full <- model %>% predict(X_full)
  ypV=ModNorm(Y_pred[,2],Y_pred[,1])
  ypLV=ModNorm(Y_predL[,2],Y_predL[,1])
  ypFV=ModNorm(Y_pred_full[,2],Y_pred_full[,1])
  SampleWeight=SampleWeight*ypLV 
  TestWeight=TestWeight*ypV
}

if(action=='Learn2')
{
  # Modèle Keras
  model <- keras_model_sequential() %>%
    layer_dense(units = 9, activation = "relu",input_shape = ncol(X_train)) %>% #288 #128
    layer_dense(units = 128, activation = "relu")%>%
    layer_dropout(rate = 0.15) %>% 
    layer_dense(units = 64, activation = "relu")%>%   #72 #64
    layer_dense(units = 12, activation = "relu") %>% #18 #16
    layer_dense(units = 2)
  
  model %>% compile(optimizer ='Adam' , loss = custom_loss,weighted_metrics=custom_loss) #"mse" "adam"
  history <- model %>% fit(X_train, Y_train, epochs = 75,batch_size=32,sample_weight=SampleWeight, validation_data = list(X_test, Y_test,TestWeight))
  model %>% evaluate(X_test, Y_test)
  
}


#Result display

SelSet=1 #Show 1: train 2:test 3:train+test 4:all o:raw data
SelDist=2 # Mode of distance calculation 1: distance to projection 2: distance to closest CFF atom
Dist_filter=c(0,7) # distance filter for match display
ColorTrh=c(3.5,5,7,Dist_filter[2]) #distance threshold 1,2,3,mask
FitTol=15 # Tolerance on fit %
UseCSPond=0.5 # -1 no ponderation, 0: proportional square  , >0 et <1 use threshold (0.75)
BuildSelfCtr=FALSE # control for self -match

# Prédictions
Y_pred <- model %>% predict(X_test)
Y_predL <- model %>% predict(X_train)
Y_pred_full <- model %>% predict(X_full)

#plot(Y_pred_full,Y_full)
#sum(abs(Y_pred_full-Y_full))/length(Y_full)

#angle back conversion
yp=atan2(Y_pred[,2],Y_pred[,1])*180/pi #pred test
ypV=ModNorm(Y_pred[,2],Y_pred[,1])
yd=atan2(Y_test[,2],Y_test[,1])*180/pi #data test
ypL=atan2(Y_predL[,2],Y_predL[,1])*180/pi #pred train
ypLV=ModNorm(Y_predL[,2],Y_predL[,1])
ydL=atan2(Y_train[,2],Y_train[,1])*180/pi #data train
ydF=atan2(Y_full[,2],Y_full[,1])*180/pi #data full
ypF=atan2(Y_pred_full[,2],Y_pred_full[,1])*180/pi #pred full
ypFV=ModNorm(Y_pred_full[,2],Y_pred_full[,1])
#plot(ypF,ydF)




if(BuildSelfCtr)
{
rr=sample(length(ypF))
ypF=ydF[rr]
rr=sample(length(yp))
yp=yd[rr]
rr=sample(length(ypL))
ypL=ydL[rr]
}

Thd1=ColorTrh[1]#<=RED
Thd2=ColorTrh[2] # <green
Thd3=ColorTrh[3] # <BLUE
Thd4=ColorTrh[4] # <mask
if(SelSet==2) tb=FullDat[ss2,] else tb=FullDat[ss1,]
dset=FullDat
if(SelDist==1) dy=sqrt(dset[,'D2']^2+dset[,'D4']^2) #define color based on distannce ofactive Oxy to projection on CFF plane
if(SelDist==2) dy=dset[,"dmin"] #alternate option: define color based OF CLOSEST cff ATOM TO of active Oxy t
dyD=dy[ss2]
dyP=dy[ss1]

col=rep(rgb(0,0,0),length(ydF))
if(!ColorOnBatch)
{
  Trps=0.5

col[dy<Thd1]=rgb(1,0,0,0.5)
col[dy<Thd2 & dy>=Thd1]=rgb(0,1,0,0.5)
col[dy>=Thd2 & dy<Thd3]=rgb(0,0,1, 0.5)
col[dy>=Thd3]=rgb(1,0.65,0.5,0.5)
col[dy>=ColorTrh[4]]=rgb(1,1,1,0)
}else{
  col=rep(rgb(0,0,1),length(ydF))
  col[MaskL]=rgb(1,0,0)
}

if(BuildSelfCtr==FALSE)
{
  for(i in 1:length(yp))
  {
    sl=c(yp[i],yp[i]+360,yp[i]-360)
    sw=which(abs(sl-yd[i])==min(abs(sl-yd[i])))
    yp[i]=sl[sw]
  }
  for(i in 1:length(ypL))
  {
    sl=c(ypL[i],ypL[i]+360,ypL[i]-360)
    sw=which(abs(sl-ydL[i])==min(abs(sl-ydL[i])))
    ypL[i]=sl[sw]
    
  }
  for(i in 1:length(ypF))
  {
    sl=c(ypF[i],ypF[i]+360,ypF[i]-360)
    sw=which(abs(sl-ydF[i])==min(abs(sl-ydF[i])))
    ypF[i]=sl[sw]
  }
}

sdp=sqrt(sum((yp-yd)^2)/length(yp))
sdl=sqrt(sum((ypL-ydL)^2)/length(ypL))
cat("Standard deviation for test and learn data fits",sdp,sdl,"\n")

pp=c(yp,ypL)
dd=c(yd,ydL)
yyc=c(ypV,ypLV)
dyc=c(dyD,dyP)
dycf=(dyc<=Dist_filter[2] & dyc>Dist_filter[1])
pps0=cbind(pp,pp+360,pp-360)
for(i in 1:length(pp))
{
  y=which(abs(pps0[i,]-dd[i])==min(abs(pps0[i,]-dd[i])))
  pp[i]=pps0[i,y]
}
ppf=pp[dycf] #learn & train DL
ddf=dd[dycf] #learn & train MD

mss=array(0,dim=3)
mss[1]=sum(abs(pp-dd)<FitTol)/length(pp) #match not filtered for distance
mss[2]=sum(abs(ppf-ddf)<FitTol)/length(pp) #match filtered for distance

ppr=pp[sample(length(pp))]
mss[3]=sum(abs(ppr-dd)<FitTol)/length(pp) # self coincidence

mssb=array(0,dim=c(length(Shift)-1,4))
mcco=mssb
sz=length(Shift)-1
dyf=(dy<=Dist_filter[2] & dy>Dist_filter[1])
hww=array(0,dim=c(sz,10))
for(i in 1:sz)
{
  mm=rep(FALSE,length(ypF))
  mm[valShift==i]=TRUE
  pp1=ypF[FiltFE & mm] #filt FE & batch
  dd1=ydF[FiltFE & mm ]
  dyf1=dyf[FiltFE & mm]
  ww=ypFV[FiltFE & mm]
  hh=hist(ww,plot=FALSE,breaks=10)$count
  wws=sort(ww)
  clz=seq(0,1,by=0.1)
  for(j in 2:11)
  {
    hww[i,j-1]=sum(wws>clz[j-1] & wws<=clz[j])
  }
  pp1n=pp1
  dd1n=dd1
  pp1=pp1[dyf1] #filt distance
  dd1=dd1[dyf1]
  ww1=ww[dyf1]
  if(UseCSPond<0)
  {
    mssb[i,1]=sum(abs(pp1n-dd1n)<FitTol)/length(pp1n)
    if(length(pp1)>0) mssb[i,2]=sum(abs(pp1-dd1)<FitTol)/length(pp1)
  }
  if(UseCSPond==0) 
  {
    rzn=ww[abs(pp1n-dd1n)<FitTol]
    rz=ww1[abs(pp1-dd1)<FitTol]
    mssb[i,1]=sqrt(sum(rzn*rzn)/sum(ww*ww)) #not distance filtered
    mssb[i,1]=sqrt(sum(rzn*rzn)/length(ww)) #not distance filtered
    mssb[i,2]=sum(rz*rz)/sum(ww1*ww1) # distance filtered
    mssb[i,2]=sqrt(sum(rz*rz)/length(ww1)) # distance filtered
  }
  if(UseCSPond>0)
  {
    rzn=rep(0,length(pp1n))
    rz=rep(0,length(pp1))
    rzn[ww>UseCSPond]=1
    rz[ww1>UseCSPond]=1
    mssb[i,1]=sum(rzn[abs(pp1n-dd1n)<FitTol])/sum(rzn)
    mssb[i,1]=sum(rzn[abs(pp1n-dd1n)<FitTol])/length(rzn)
    mssb[i,2]=sum(rz[abs(pp1-dd1)<FitTol])/sum(rz)
    mssb[i,2]=sum(rz[abs(pp1-dd1)<FitTol])/length(rz)
  }
}
#barplot(t(hww),beside=TRUE, col=rainbow(10))
plot(y=hww[1,],x=seq(0.05,1,by=0.1),type="b",col=1,ylim=c(0,max(hww)),pch=16,xlab="Radius",ylab="count")
for(i in 2:sz)
{
  lines(y=hww[i,],x=seq(0.05,1,by=0.1),type="b",col=i,pch=16)
}
legend("topleft", legend=c(paste("Batch ",1:sz)), col=1:sz, lty=1:1, cex=1)
       
SumRes=cbind(mss[1],mss[2],mss[3],mssb[,1],mssb[,2])
colnames(SumRes)=c("GMatch","FMatch","CMatch","PMatch1","PFMatch2")
rownames(SumRes)=nname
rownames(SumRes)=SetNames
SumRes[SumRes<0]=0
print(SumRes)
sd=sort(SumRes[,4],decreasing=TRUE,index.return=TRUE)
colr=rainbow(sz)


par(mar=c(8,3,2,0))

barplot(SumRes[,4], col=colr,ylim=c(0,1),las=2,cex.name=1.3) #not distance filtered  
legend(x="topright",legend = rownames(SumRes),cex=1,fill=colr)  

barplot(SumRes[,c(5)],beside = TRUE, col=colr,ylim=c(0,1),las=2,cex.name=1.3) #distance filtered 
legend(x="topright",legend = rownames(SumRes),cex=1,fill=colr) 

# self coincidence from randomization
#barplot(SumRes[,6],beside = TRUE, col=colr,ylim=c(0,1))
#legend(x="topright",legend = rownames(SumRes),cex=1,fill=colr)  


#--------------------------------------------end plot phase I
LimQ=c(0,Thd1,Thd2,Thd3,15)
tp0=hist(dd,breaks=17,plot=FALSE)
ycc0=rep(1,length(dd))
if (UseCSPond==0) ycc0=yyc*yyc
if (UseCSPond>0) ycc0[yyc<UseCSPond]=0 

for(i in 1:4)
{
  filt=(dyc<=LimQ[i+1] & dyc>LimQ[i])
  dd1=dd[filt]
  pp1=pp[filt]
  yyc1=ycc0[filt]
  tt=dd1[abs(pp1-dd1)<FitTol]
  tty=yyc1[abs(pp1-dd1)<FitTol]
  cut1=cut(tt,breaks=tp0$breaks,ordered_result = TRUE,labels=FALSE)
  cut2=cut(dd1,breaks=tp0$breaks,ordered_result = TRUE,labels=FALSE)
  tp1=rep(0,length(tp0$mids))
  tp2=rep(0,length(tp0$mids))
  for(j in 1:length(tp0$mids))
  {
    #if(is.na(sum(cut1==j))) next
    tp1[j]=sum(tty[cut1==j])
    #if(is.na(sum(cut2==j))) next
    tp2[j]=sum(yyc1[cut2==j])
  }
  if(i==1)
  {
    success=tp1
    trials=tp2
  }else{
    success=cbind(success,tp1)
    trials=cbind(trials,tp2)
  }
}
res=success/trials
res[is.na(res)]=0
ssm=apply(t(trials)[1:4,],2,sum)
tt=t(success)[1:4,]
for(i in 1:dim(tt)[2]) tt[,i]=tt[,i]/ssm[i]*100
tt[is.na(tt)]=0


#Success rate plot by dmin (color threshold)
barplot(tt[1,],names=tp0$mids,ylim=c(0,100),las=2,col=c("red"))
barplot(tt[2,],names=tp0$mids,ylim=c(0,100),las=2,col=c("green"))
barplot(tt[3,],names=tp0$mids,ylim=c(0,100),las=2,col=c("blue"))
barplot(tt[4,],names=tp0$mids,ylim=c(0,100),las=2,col=c("grey"))
#success rate single plot
barplot(tt,names=tp0$mids,las=2,col=c("red","green","blue","grey"), ylim=c(0,100))
#population single plot
barplot(t(trials),names=tp0$mids,las=2,col=c("red","green","blue","grey"))
#---------------------------------------------------------------------------End plot phase II
dyp1=(dyP<=Dist_filter[2] & dyP>Dist_filter[1])
dyd1=(dyD<=Dist_filter[2] & dyD>Dist_filter[1])
if(SelSet==2 ) plot(yp,yd,col=col[ss2],xlim=c(-180,180),ylim=c(-180,180))
if(SelSet==1 ) plot(ypL,ydL,col=col[ss1],xlim=c(-180,180),ylim=c(-180,180)) 

if(SelSet==4 )
{
  dyf=(dy<=Dist_filter[2] & dy>Dist_filter[1])
  pp=ypF[dyf] #all predicted 
  dd=ydF[dyf] #all data
  ccol=col[dyf]
  vv=ypFV[dyf]
}
if(SelSet==3 )
{
  pp=c(yp[dyd1],ypL[dyp1]) #all predicted 
  dd=c(yd[dyd1],ydL[dyp1]) #all data
  ccol=c(col[ss2][dyd1],col[ss1][dyp1])
  vv=c(ypV[dyd1],ypLV[dyd1])
}
if(SelSet==2 )
{
  pp=yp[dyd1] 
  dd=yd[dyp1] 
  ccol=col[ss2][dyd1]
  vv=ypV[dyd1]
}
if(SelSet==1 )
{
  pp=ypL[dyp1]
  dd=ydL[dyp1]
  ccol=col[ss1][dyp1]
  
}
  rr=sample(length(pp))
  xx=pp[rr]
  yy=dd[rr]
  cc=ccol[rr]
  vx=1:length(xx)
  if(UseCSPond>0) vx=vv[rr]>UseCSPond
  plot(xx,yy,col=cc,xlim=c(-180,180),ylim=c(-180,180),pch=1,cex=0.4)
  plot(xx[vx],yy[vx],col=cc[vx],xlim=c(-180,180),ylim=c(-180,180),pch=1,cex=0.4)
  
  mean(xx)
  sd(xx[abs(xx)<50])
  plot(xx[col=="#FF0000"],yy[col=="#FF0000"],col="red",xlim=c(-180,180),ylim=c(-180,180),pch=1,cex=0.4)
  plot(xx[col=="#0000FF"],yy[col=="#0000FF"],col="blue",xlim=c(-180,180),ylim=c(-180,180),pch=1,cex=0.4)

#------------------------------------------------- End plot phase III

ccl=c("red","green","blue","orange","black")
plot(c(),xlim=c(-180,180),ylim=c(-180,180))
for(i in 1:(length(Shift)-1))
{
  mm=rep(FALSE,length(ypF))
  mm[Shift[i]:Shift[i+1]]=TRUE
  points(ypF[FiltFE & mm],ydF[FiltFE &mm],col=ccl[i])
 # points(ypF[FiltFE &mm],ydF[FiltFE &mm],col=col[FiltFE &mm])
}

#------------------------------------------------
# mapping by overlay
xs0=1:dim(Env)[1]
xs=xs0[ss2]
xt=xs0[ss1]
xf=1:length(ydF)
plot(x=NULL,xlim=c(1,dim(Env)[1]),ylim=c(-180,180))
if(SelSet==2 || SelSet==3)
{
  for(i in 1:max(length(yd),length(ypL)))
  {
    #if(i<=length(yd)) points(xs[i],yd[i],col=col[ss2][i])
    if(i<=length(yd)) points(xs[i],yd[i],col="black")
    if(i<=length(ypL)) points(xt[i],ypL[i],col=col[ss1][i])
  }
}

if(SelSet==1 || SelSet==3)
{
  for(i in 1:max(length(yd),length(ypL)))
  {
    if(i<=length(yp)) points(xs[i],yp[i],col=col[ss2][i])
    if(i<=length(ydL)) points(xt[i],ydL[i],col=col[ss1][i])
  }
}


if(SelSet==0)
{
  plot(x=NULL,xlim=c(1,dim(Env)[1]),ylim=c(-180,180))
  points(xf,ydF,col=col)
}
for (i in 1:length(limits)) lines(c(limits[i],limits[i]),c(-180,180), col="blue",lty='dashed')

 
#-------------------------------------------------------------------------

  
nclass=13
ccutD=cut(yd,breaks=nclass,ordered_result = TRUE,labels=FALSE)
ccutP=cut(yp,breaks=nclass,ordered_result = TRUE,labels=FALSE)
acutD=array(0,dim=nclass)
acutP=array(0,dim=nclass)
for(i in 1:length(ccutD))
{
  acutD[ccutD[i]]=acutD[ccutD[i]]+1
  if(ccutD[i]==ccutP[i]) acutP[ccutD[i]]=acutP[ccutD[i]]+1
}
r=acutP/acutD*100
lv=seq(-180,180,length=nclass+1)
lb=array(0,dim=nclass)
for(i in 1:nclass) lb[i]=paste(trunc(lv[i]),trunc(lv[i+1]),sep=":")
barplot(r,names=lb,las=2)
locator(1)
barplot(acutD,names=lb,las=2)
am=mean(r[acutD>0.1*max(acutD)])
plot(acutD,r)

col1=rep(rgb(0,1,0),length(yp))
col1[dy<Thd]=rgb(1,0,0)
col1[dy>10]=rgb(0,0,1)
plot(lowess(xd,yd,f=0.02,iter=2),ty='p')
points(lowess(xd,yp,f=0.02,iter=2),col=col1,ty='p')

plot(X,Ym[,1], ylim=c(min(Ym[,1]),max(Ym[,1])),ty='p')
points(X,Y_pred[,1], ylim=c(min(Y_pred[,1]),max(Y_pred[,1])), col="red",ty='p')

plot(lowess(X,Y,f=0.01,iter=2), ylim=c(min(Y),max(Y)),ty='l')
points(lowess(X,Y_pred,f=0.01,iter=2), ylim=c(min(Y_pred),max(Y_pred)), col="red",ty='l')

}



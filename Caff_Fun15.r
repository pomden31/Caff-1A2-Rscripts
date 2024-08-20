library(bio3d)
library(pracma)
library(stringr)
library(plyr)
library(testit)
library(R2HTML)
library(fitdistrplus)
library(corrplot)
library(rmarkdown)
library (mixtools)
library(vioplot)
library(FactoMineR)
library(ggplot2)
library (dplyr)
library(plotrix)

source("Caff_UtilA15.R")
source("Caff_UtilB15.R")


SetMainParm=function(Prefix,Range,UseCountMode=FALSE)#Description File location
{
  WorkSetFile<<-"Files/WorkSet.txt"
  IndexFile<<-"Files/Index.txt"
  VarListFile<<-"Files/VarList.txt"
  SubvarListFile<<-"Files/SubVarList.txt"
  IndexFile<<-"Files/Index.txt"
  TimeSeriesFile<<-"Files/TimeSeries.txt"
  
  
  SetPrefix<<-Prefix
  SumResLoc<<-paste(SetPrefix,"sumres.txt",sep="")
  SumFramLoc<<-paste(SetPrefix,"sumfram.txt",sep="")
  SumSumloc<<-paste(SetPrefix,"sumsum.txt",sep="")
  SumLoc<<-paste(SetPrefix,"sum.txt",sep="") 
  NormDataLoc<<-paste(SetPrefix,"NormData.txt",sep="") 
  AngRec<<- array(0,dim=c(100,3))
  #end
  
  #running parameters for phase I
  setrange<<-Range
  FullMode<<-TRUE # When true override the frame windows defined in "worset.txt and MaskOverlay setting
  MaskOverlay<<-"50:800" # (Negative values are interpreted as relative to last frame). Use NULL for setting using "Workset.txt" file values
  #AnalyseCFFInteract<<-FALSE # When TRUE analyse interaction between two CFF when present
  DoPrint<<-TRUE #Turn on detailed statistic printing
  PrintMode1<<-1 # NULL: No detailed Print, Normal print for any numeric value 
  SetReadRefPDB<<-TRUE  # FALSE en standard. When TRUE the 10 first frame are initial condition before running MD/ mask limits  are not taken in account
  if(SetReadRefPDB) MaskType<<-0 else  MaskType<<-NULL # must be 0 when ReadRefPDB=TRUE or NULL for standard mode 
  AltMode<<-UseCountMode # When TRUE use the "alternative "count mode" for CFF analysis of orientations
  AltProcSet<<-c(22,23,24) #IDENA,IDENB,IDENC #Index of experimental parameters reported by values (1-9) and not by statistics
  #end
  return
}


#start of phase I
RunPhaseI=function(setrange,DoWrite=TRUE)
{
NewBatch<<-TRUE
SumResDat=data.frame()
SumResTest=data.frame()
SumResFrame=data.frame()
Sumsum=data.frame()
NorData=data.frame()
NormData=data.frame()
ws=read.table(WorkSetFile,header=TRUE,sep='\t')
idx=read.table(IndexFile,header=TRUE,sep='\t')

for (u in setrange)
{
  #skip invalid frames
  if(as.numeric(ws$CFF[u])==0) next
  if(idx$Frames[as.numeric(ws$Id[u])]==0) next
  
  #setting PDB read when used
  tt=file.exists(paste(ws$OMM[u],"RefPDB.pdb",sep="/"))
  if(SetReadRefPDB && tt) ReadRefPDB<<-TRUE else ReadRefPDB<<-FALSE
  if(ReadRefPDB) MaskType<<-0 else MaskType<<-NULL
  
  # detecting 2 CFF bound
  code=idx$Struct[ws$Id[u]]
  if(length(grep("c",code))==1 && length(grep("C",code))==1 && AnalyseCFFInteract) twoCFF=TRUE else twoCFF=FALSE
  
  #setting name
  FN=paste(u,ws$OMM[u],"_Id_",ws$Id[u],"CFF_",ws$CFF[u],".pdf",sep="")
  FN1=paste(u,ws$OMM[u],"_Id_",ws$Id[u],"CFF_",ws$CFF[u],".txt",sep="")
  pdf(file=FN, width=7,height=5)
  plot(1:100,1:100,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  text(1,90,paste("Sample Id=",ws$Id[u]),adj = c(0,0))
  text(10,80,paste("Sample name:",FN),adj = c(0,0))
  if(length(grep("c",code))==1 && length(grep("C",code))==1)
  {
    text(10,70,"CFF at active and effector sites",adj = c(0,0))
    NbCFF=2
  }
  else{
    text(10,70,"Single CFF",adj = c(0,0))
    NbCFF=1
  } 
  if(ws$CFF[u]==1) text(10,60,"Effector CFF selected",adj = c(0,0))
  else text(10,60,"Productor CFF selected",adj = c(0,0))
  if(ws$Target[u]=="FE") text(10,50,"Heme reference atom iron",adj = c(0,0))
  else text(10,50,"Heme reference active oxygen",adj = c(0,0))
  text(10,40,paste ("Flavin redox state :",ws$FMN[u]),adj = c(0,0))
  text(10,30,paste ("OMMId :",ws$OMM[u]),adj = c(0,0))
  
  
  #Setting mask
  Mmask=as.numeric(unlist(str_split(MaskOverlay,":")))
  if(sum(Mmask)<0 && FullMode==FALSE)
  {
    TFR=idx[idx$ID==ws$Id[u],"Frames"]
    Mmask[1]=max(Mmask[1],1)
    Mmask[2]=max(Mmask[2],2)
    Mmask=as.character(Mmask)
    MaskOverlay=paste(Mmask[1],":",Mmask[2], sep="")
  }
  if(!is.null(MaskOverlay)) mmode=MaskOverlay
  if(FullMode) mmode="1:9999" else mmode=ws$mode[u]
  # end setting mask
  
  #res=AnalyseCFF(ws$Id[u],ws$Target[u],ws$CFF[u],mmode, FN,twoCFF,MaskType,PrintMode1, ReadRefPDB,DoPrint, FN1)
  if((SinglePDBRef&& u==setrange[1]) || !SinglePDBRef)
  {
    table=read.table(IndexFile, sep='\t')
    colnames(table)=table[1,]
    table=table[2:dim(table)[1],]
    aSet=c(paste("Set",ws$Id[u], sep='_'),table[ws$Id[u],2])
    if(!ReadRefPDB) PDBRef=read.pdb(paste(aSet[2],"/data.pdb", sep="")) else{
      PDBRef=read.pdb(paste(aSet[2],"/RefPDB.pdb", sep=""))}
  }
  res=AnalyseCFF(ws$Id[u],ws$Target[u],ws$CFF[u],mmode, FN,twoCFF,VarListFile,SubvarListFile,IndexFile,TimeSeriesFile,
                 PDBRef,MaskType,PrintMode1, ReadRefPDB,DoPrint, FN1, AltMode)
  Tp=res[[8]]
  Id=array(c(ws$Id[u],ws$Target[u],ws$CFF[u],NbCFF),dim=c(4,dim(Tp)[1]))
  rownames(Id)=c("Id","Fe/O","CFF_Id","CFF_Nb")
  TP=cbind(t(Id),Tp)
  NormData=rbind(NormData,TP)
  ss=PrintRes(res, paste(ws$OMM[u],"_CFF",ws$CFF[u],".txt", sep=""))
  if(AltMode) AltProc=AltProcSet else AltProc=NULL 
  ssm=SumSS(ss,paste(ws$OMM[u],"_caff",ws$CFF[u], sep=""), AltProc)
  SumResFrame=rbind(SumResFrame,data.frame(id=u,OMM=ws$Id[u],ws$mode[u],CFFTyp=ws$CFF[u]))
  dev.off()
  SumResDat=rbind(SumResDat,ss) 
  if(u==setrange[1]) Sumsum=ssm else Sumsum=cbind(Sumsum,ssm)
  cat(u)
}
if(DoWrite)
{
write.table(SumResDat,SumResLoc)
write.table(SumResFrame,SumFramLoc) 
write.table(Sumsum,SumSumloc) 
write.table(NormData,NormDataLoc)
}

#Fuze previous table into a summary table ResFile

Sumsum=read.table(paste(SetPrefix,"sumsum.txt",sep=""))
SumResFrame=read.table(paste(SetPrefix,"sumfram.txt",sep=""))
expname=paste(ws$OMM[setrange],"CFF_",ws$CFF[setrange],sep="")
expname=expname[as.numeric(ws$CFF[setrange])!=0]
if(AltMode)
{
  AltProc=AltProcSet 
  select=c(3,4,5) # 3: Gauss fit multiple main extended
}else{
  AltProc=NULL
  select=c(3,4,5) # 3: Gauss fit multiple main
} 
rr=FuzeResultSet2(Sumsum,SumResFrame, expname,SumLoc,select, AltProc) #sum.txt
}
#end of phase I
#sink(NULL)

SetPhaseIIParm=function(NumbCFF,SelectCFF)
{
  #start of phase II (independent of phase I when result files of phase I are available)
  #Setting parameters for phase II
  NbCFF<<-NumbCFF #Number of CFF ( 1 or 2)
  SelCFF<<-SelectCFF #selected CFF (1: effector, 2 producer)
  WorkSetFile<<-"Files/WorkSet.txt" 
  IndexFile<<-"Files/Index.txt"     
  PdfSet<<-paste(SetPrefix,"sum.txt" ,sep="") #pdf file name
  SortSelection<<-FALSE # activate data selection sorting on criteria
  Bar<<-TRUE   # activate use of bar plot
  UseRainbowPlot<<-TRUE # select TRUE for rainbow and violon plots
  AltMode<<-TRUE # When TRUE use the alternative mode for CFF analysis of orientations
  AltProcSet<<-c(22,23,24) #IDENA,IDENB,IDENC #Index of experimental parameters reported by values (1-9) and not by statistics
  ParmFile<<-"Files/PlotParms.txt"
  #end of phase II
  return
}

RunPhaseIIA=function()
{
InitPos<<-NULL
InitOxDist<<-NULL
WorkSet<<-WorkSetFile
if(UseRainbowPlot)
{
  MeanType=c(1,2,3) # Select record of sum txt used as mean (or else)
  ErrorType=c(4,5,6)  # Select record of sum txt used as Sd (or else)
  CountType=c(7,8,9)  # Select record of sum txt used as Sd (or else)
  OrientType=c(1,2,3,4,5,6,7,8,9) # Select record of counts instead of stats (NULL for basic mode)
  if(AltMode)
  {
    InitType=11
    OrientType=c(1,2,3,4,5,6,7,8,9) # Select record of counts instead of stats (NULL for basic mode)
  }else{
    InitType=11  # Select record of sum txt used as Sd (or else)
    OrientType=NULL
  }
}else {
  MeanType=1 # Select record of sum txt used as mean (or else)
  ErrorType=4  # Select record of sum txt used as Sd (or else)
  CountType=7  # Select record of sum txt used as Sd (or else)
  InitType=11
  OrientType=NULL
}

#extract data subset from tab
tab=as.matrix(read.table(paste(SetPrefix,"sum.txt",sep=""), header = TRUE))
asel=selectRec(PdfSet, WorkSet, NbCFF, SelCFF, InitPos, InitOxDist, MeanType, ErrorType, CountType, InitType,OrientType) # sel OMM on set of criteria

if(SortSelection) #optional sort order of selecton on the given criteria
{
  ss=sort(tab[17,asel$sel], index.return=TRUE,decreasing=TRUE) # optional sort selected criteria
  Sort=ss$ix
  tab=tab[,Sort]
}
#optional selection modes : selectRec1=function(PdfSet, WorkSet, OMMId, SelCFF=1) # select for specific OMM
#                         SelectOMMSet(): select for specific OMM group (manual) #SelectOMMSet()

ws=read.table(WorkSet, sep='\t', header = TRUE)
tab1=tab[,asel[,1]]
tab2=tab1[1:4,tab1[1,]==1]
Stab1=tab[,asel[,2]]
Ntab1=tab[,asel[,3]]
Itab1=tab[,asel[,4]]
ds=dim(tab1)[1]-1
mat=tab1[5:ds,] #mat and matSD are the extracted data subset
matSD=Stab1[5:ds,]#mat and matSD are the extracted data subset
matNb=Ntab1[5:ds,]#mat and matNb are the extracted data subset
tabs=tab[c(26,27,28),asel[,5]]
colnames(tabs)=rep(c(1,2,3,4,5,6,7,8,9),dim(tabs)[2]/9)
if(!is.null(dim(Itab1))) matInit=Itab1[5:ds,] #raw matrix of initial values of parameters
else matInit=Itab1[5:ds]
if(!is.null(dim(matInit))) MatInits=matInit[c(22,23,24),] # matrix of initial values of parameters
else MatInits=matInit[c(22,23,24)] 

# set and printlist of used OMMs
re=c()
for (i in 1:length(tab2[3,]))
{
  t=which(ws[,2]==tab2[3,i])
  if(as.numeric(ws$CFF[min(t)])==0) next #skip OMM not considered in scan
  re=c(re,ws[min(t),3])
}
print(re) # list of used OMMs
return(list(mat,matSD,matNb,matInit,MatInits,tabs,re))
}

SetPhaseIIBParm=function(MyItem,SetPrefix,SelectCFF)
{
Item=MyItem # select parameter to compare between OMM runs
CffId=SelectCFF #1:effector, 2:producer USE RESTRICTED TO PRINT PARM SETTING (not md selection)
MapSize<<-500
Dset=c("Serie671/","Serie669/", "Serie671_1CFF/","Serie669_1CFF/","MixLL/","MixLN/","MixRL/","MixRN/","MixRR/")
Dataset=which(Dset==SetPrefix)
MyParms=GetPlotParms(ParmFile,Dataset,CffId,Item)
if(!is.null(MyParms) && dim(MyParms)[1]>0)
{
  AsStaking=MyParms$AsStaking
  Limits=MyParms$Limits
  ScaleY=MyParms$ScaleY
  AsAngle=MyParms$AsAngle
  Ofset=MyParms$Ofset
  ylab=MyParms$ylabs
  AsCOC=MyParms$AsCOC
  DOfset=eval(parse(text = MyParms$DOfset)) #ofset on specific angles
  ROfset=eval(parse(text = MyParms$Rofset)) #ofset on specific angles
}else{
  YlabV=c("Distance  ( Å )","Angle  ","Staking extend (%)","Orientation code","")
  Limits="c(-500,500)" 
  ScaleY="c(0,360)"
  Ofset=0
  AsAngle=FALSE
  AsCOC=FALSE
  AsStaking=FALSE
  if(AsAngle)
    {
    ylab=YlabV[2]
    ScaleY="c(0,360)"
    }else{
    ylab=YlabV[1]
    ScaleY="c(0,15)"
    Limits="c(0,20)" 
    } 
  if(AsCOC)
    {
    ylab=YlabV[4]
    ScaleY="c(0,18)"
    }
  if(AsStaking)
    {
    ylab=YlabV[3]
    ScaleY="c(0,100)"
    Limits="c(0,100)" 
    }
  DOfset=c(0,0,0,0,0,0,0,0,0,0,0,0)
  ROfset=c(0,0,0,0,0,0,0,0,0,0,0,0)
}
return(list(Limits=Limits,ScaleY=ScaleY,Ofset=Ofset,AsAngle=AsAngle,AsCOC=AsCOC,AsStaking=AsStaking,ylab=ylab,DOfset=DOfset,ROfset=ROfset,Index=MyParms$Index ))
}


RunPhaseIIB=function(Datamat,res1, Item,Sync=FALSE)
{
  AsCOC=res1[[5]]
  AsAngle=res1[[4]]
  ScaleY=res1[[2]]
  AsCOC=res1[[5]]
  Limits=res1[[1]]
  DOfset=res1[[8]]
  mat=Datamat[[1]]
  matSD=Datamat[[2]]
  matNb=Datamat[[3]]
  matInit=Datamat[[4]]
  MatInits=Datamat[[5]]
  tabs=Datamat[[6]]
  if(Sync)
  {
    ival=Datamat[[4]]
    of=c(0,360,-360)
    for (i in 1:dim(mat)[1])
    {
      for (j in 1:dim(ival)[2])
      {
        a=mat[i,1+3*(j-1)]
        c=ival[i,j]
        if(c<0)
        {
          c=c+360
          Datamat[[4]][i,j]=c
        }
        if(c>360)
        {
          c=c-360
          Datamat[[4]][i,j]=c
        }
        b=c(a,a+360,a-360)
        d=which(abs(b-c)==min(abs(b-c)))
        mat[i,1+3*(j-1)]=mat[i,1+3*(j-1)]+of[d]
        mat[i,2+3*(j-1)]=mat[i,2+3*(j-1)]+of[d]
        mat[i,3+3*(j-1)]=mat[i,3+3*(j-1)]+of[d]
      }
    }
  }
  if(!Item%in%c(22,23,24) || !AltMode)
  {
    val=PrepStatplot(mat,matSD,matNb, Item,MapSize,AsAngle,AsCOC,DOfset)
  }else{
    nr=dim(tabs)[2]/9
    val1=array(NA,dim=c(MapSize,3*nr))
    for(i in 1:nr)
    {
      val=array(NA,dim=c(MapSize,3))
      tt=tabs[Item-21,(1+9*(i-1)):(9*i)]
      nt=sum(tt)
      tmat=c()
      for (j in 1:9)
      {
        tmat=c(tmat,rep(j,trunc(tt[j]*MapSize/nt)))
      }
      val1[1:length(tmat),i]=tmat
    }
    val=val1
  } 
  if(AsCOC)
  {
    lab=c("C2(C12N3)","C2N4","N4(N9)","C8(C8)","C3(C14N7)","C3O2","O2(O13)","C1(C10N1)","O1(O11)")
  }else lab=NULL
  sc=eval(parse(text = ScaleY)) 
  li=eval(parse(text = Limits))
  return(list(val,sc,li,lab,Datamat))
}

RunPhaseIIBPlot=function(Mode,Item,res,res1,res2)
{
  val=res2[[1]]
  li=res2[[3]]
  sc=res2[[2]]
  ylab=res1[[7]]
  mat=res[[1]]
  matInit=res[[4]]
  lab=res2[[4]]
  Ofset=res1[[3]]
  Rofset=res1[[9]]
  DOfset=res1[[8]]
  par(mfrow=c(1,1))
  if(Mode==1) RainbowPlot(val,li,sc,"Rotation angle  ",ylab,paste("Parm",rownames(mat)[Item], sep=" "), Item,matInit,lab, Ofset,Rofset )
  if(Mode==2)ViolonPlot(val,li,sc,"Rotation angle  ",ylab,paste("Parm",rownames(mat)[Item], sep=" "),Item,matInit,lab,Ofset,Rofset)
  if(Mode==3)
  {
    par(mfrow=c(1,1))
    RainbowPlotSP(val,li,sc,"Rotation angle  ",ylab,paste("Parm",rownames(mat)[Item], sep=" "), Item,matInit,lab,Ofset,DOfset,0)
  }
  par(mfrow=c(1,1))
  dev.flush()
}

SavePhaseIIBPlotParm=function(Item,cffId,res,res1,res2)
{
  Rofset=res1[[9]]
  DOfset=res1[[8]]
  Limits=res1[[1]]
  ScaleY=res1[[2]]
  AsCOC=res1[[5]]
  AsAngle=res1[[4]]
  Ofset=res1[[3]]
  ylab=res1[[7]]
  
  MatInits=res[[5]]
  tabs=res[[6]]
  Dset=c("Serie671/","Serie669/", "Serie671_1CFF/","Serie669_1CFF/","MixLL/","MixLN/","MixRL/","MixRN/","MixRR/")
  Dataset=which(Dset==SetPrefix)
  par(mfrow=c(1,2))
  RO= paste("c(",paste(Rofset,collapse=","),")",sep="")
  DO= paste("c(",paste(DOfset,collapse=","),")",sep="")
  if(is.null(res1$"MyParms$Index")) res1$"MyParms$Index"=0
  nn=data.frame(Index=res1$"MyParms$Index",DataSet=Dataset,Parm=Item,CFF=SelCFF,
                  Limits=Limits,ScaleY=ScaleY,
                  AsAngle=AsAngle,AsCOC=AsCOC,
                  Ofset=Ofset,ylabs=ylab,DOfset=DO,Rofset=RO)
  SetPlotParms(ParmFile,Dataset,cffId,Item,nn)
  return
}

RunPhaseIII=function(res,mode,Multi,ShowCorrelation, xlabel,SortParm=NULL, LabScale=1)
{
  mat=res[[1]]
  matSD=res[[2]]
  re=res[[7]]
  
  if(LabScale==1 || LabScale==3) par(mfrow=c(2,2))
  if(LabScale==2) par(mfrow=c(3,3))
  Sets=DefOrientSet(mat, 5,180,Multi)
  CorData=array(dim=c(21,dim(mat)[2]/Multi))
  Cols=c(8,5,11,7,6,9,12,10,1,13,2,14,3,15,4,16,17,18,19,20,21) 
  colnames(res[[1]])=paste(xlabel,colnames(res[[1]]),sep="-")
  rownames(CorData)=rownames(mat)[Cols]
  Multi=3 # nombre de colonnes "means"
  DS=SortParm # sort Parm NULL: no sort 1:individual sort array: sort order
  if(!is.null(SortParm) && SortParm>0)
  {
   
    t=SplitPrintStat(mat,matSD,xlabel,SortParm, Sets,0, DS, FALSE, "Sorting run","","ANG",Multi,DS,0)
    DS=t$SortIdx
  }
  if(mode[1])
  {
    R1=SplitPrintStat(mat,matSD,xlabel,8, Sets,0, DS, FALSE, "Angle (A0) of LI and heme N2->N4 axis","Angle (A0)","none",Multi,DS,LabScale)
    CorData[1,]=R1[[8]]
  }
  if(mode[2])
  {
    Ofset=0
    R2=SplitPrintStat(mat,matSD,xlabel,5, Sets,0, DS, FALSE, "Angle (A1) of CFF and heme planes","Angle A1  ","none",Multi,DS,LabScale)
    CorData[2,]=R2[[8]]
  }
  if(mode[3])
  {
    Ofset=0
    R3=SplitPrintStat(mat,matSD,xlabel,11, Sets,0, DS, FALSE, "Angle CFF_C2-M3 with M2-M3 (A2) ","Angle A2  ","none",Multi,DS,LabScale)
    CorData[3,]=R3[[8]]
  }
  if(mode[4])
  {
    Ofset=0
    R4=SplitPrintStat(mat,matSD,xlabel,7, Sets,0, DS, FALSE, "Angle M2-MP1 with M3-M2 (A3) ","Angle A3  ","none",Multi,DS,LabScale)
    CorData[4,]=R4[[8]]
  }
  
  if(mode[5])
  {
    R5=SplitPrintStat(mat,matSD,xlabel,6, Sets,0, DS, FALSE, "Distance (D2) of FE/O to CFF plane","Distance D2  ( Å )","none",Multi,DS,LabScale)
    CorData[5,]=R5[[8]]
  }

  if(mode[6])
  {
    R6=SplitPrintStat(mat,matSD,xlabel,9, Sets,0, DS, FALSE, "Distance (D4) of CFFcen to FE/O proj.","Distance D4 ( Å )","none",Multi,DS,LabScale)
    CorData[6,]=R6[[8]]
  }
  if(mode[7])
  {
    R7=SplitPrintStat(mat,matSD,xlabel,12, Sets,0, DS, FALSE, "Smaller dist. (dmin) of CFF atom to FE/O.","Distance dmin  ( Å )","none",Multi,DS,LabScale)
    CorData[7,]=R7[[8]] 
  }
  if(mode[8])
  {
    R8=SplitPrintStat(mat,matSD,xlabel,10, Sets,0, DS, FALSE, "CFF orientation code (imin)","","ORI",Multi,DS,LabScale)
    CorData[8,]=R8[[8]]
  }
  
  if(mode[9])
  {
    R9=SplitPrintStat(mat,matSD,xlabel,1, Sets,0, DS, FALSE, "Stacking coverage Phe202(226)","% overlap","none",Multi,DS,LabScale)
    CorData[9,]=R9[[8]]
  }
  if(mode[10])
  {
    R10=SplitPrintStat(mat,matSD,xlabel,13, Sets,0, DS, FALSE, "Stacking coverage Phe236(260)","% overlap","none",Multi,DS,LabScale)
    CorData[10,]=R10[[8]]
  }
  if(mode[11])
  {
    R11=SplitPrintStat(mat,matSD,xlabel,2, Sets,4, DS, FALSE, "Angle between CFF and Phe202(226) planes","Angle  ","ANG",Multi,DS,LabScale)  
    CorData[11,]=R11[[8]]
  }
  if(mode[12])
  {
    R12=SplitPrintStat(mat,matSD,xlabel,14, Sets,4, DS, FALSE, "Angle between CFF and Phe236(260) planes","Angle  ","ANG",Multi,DS,LabScale)  
    CorData[12,]=R12[[8]]
  }
  if(mode[13])
  {
    R13=SplitPrintStat(mat,matSD,xlabel,3, Sets,0, DS, FALSE, "CFFcen to Phe202(226) plane distance", "Distance  ( Å )","none",Multi,DS,LabScale) 
    CorData[13,]=R13[[8]]
  }
  if(mode[14])
  {
    R14=SplitPrintStat(mat,matSD,xlabel,15, Sets,0, DS, FALSE, "CFFcen to Phe226(260) plane distance", "Distance  ( Å )","none",Multi,DS,LabScale) 
    CorData[14,]=R14[[8]]
  }
  if(mode[15])
  {
    R15=SplitPrintStat(mat,matSD,xlabel,4, Sets,0, DS, FALSE, "Distance between CFF and Phe202(226)\n centromer projections on Phe202 plane", "Distance  ( Å )","none",Multi,DS,LabScale)
    CorData[15,]=R15[[8]]
  }
  if(mode[16])
  {
    R16=SplitPrintStat(mat,matSD,xlabel,16, Sets,0, DS, FALSE, "Distance between CFF and Phe236(260)\n centromer projections on Phe236 plane", "Distance  ( Å )","none",Multi,DS,LabScale)
    CorData[16,]=R16[[8]]
  }
  if(mode[17])
  {
    R17=SplitPrintStat(mat,matSD,xlabel,17, Sets,0, DS, FALSE, "Dist. cenA-cenB","Distance  ( Å )","none",Multi,DS,LabScale)
    CorData[17,]=R17[[8]]
  }
  if(mode[18])
  {
    R18=SplitPrintStat(mat,matSD,xlabel,18, Sets,0, DS, FALSE, "Angle of heme and CFF1 planes","Angle  ","ANG",Multi,DS,LabScale)
    CorData[18,]=R18[[8]]
  }
  if(mode[19])
  {
    R19=SplitPrintStat(mat,matSD,xlabel,19, Sets,0, DS, FALSE, "Angle of heme and CFF2 planes","Angle  ","ANG",Multi,DS,LabScale)
    CorData[19,]=R19[[8]]
  }
  if(mode[20])
  {
    R20=SplitPrintStat(mat,matSD,xlabel,20, Sets,0, DS, FALSE, "Angle C2>CFF1 cen/cen>Proj. Fe/O","Angle  ","none",Multi,DS,LabScale)
    CorData[20,]=R20[[8]]
  }
  if(mode[21])
  {
    R21=SplitPrintStat(mat,matSD,xlabel,21, Sets,0, DS, FALSE, "Angle C2>CFF2 cen/cen>Proj. Fe/O","Angle  ","none",Multi,DS,LabScale)
    CorData[21,]=R21[[8]]
  }

   if(ShowCorrelation)
   {
    par(mfrow=c(1,1))
    rownames(CorData)=c("A0","A1","A2","A3","D2","D4","dmin","imin","SCP1","SCP2","ACP1","ACP2","DCP1","DCP2","DCP1P","DCP2P","DC12","AHC1","AHC2","ORC1O","ORC2O")
    vx=apply(CorData, 1,sd)
    CorData=CorData[vx!=0,]
    TCorData=t(CorData)
    Mt=cor(TCorData,method="kendall")
    par(mfrow=c(1,1))
    test=cor.mtest(TCorData, method="kendall",conf.level = 0.95)
    corrplot(Mt, p.mat = test$p, sig.level = 0.01,addCoef.col ='white',type = 'lower',insig = "blank",number.cex = 0.5) # p-value 0.01 for 95% confidence
    corrplot.mixed(Mt, p.mat = test$p, sig.level = 0.01,bg="grey",insig = "blank",number.cex = 0.5,tl.cex=0.5) # p-value 0.01 for 95% confidence
   }
    return(CorData)
}

  RunPhaseIIIB=function(res,mode,Multi,ShowCorrelation,xlabel,SortParm=NULL)
  {
    mat=res[[1]]
    matSD=res[[2]]
    re=res[[7]]
    
    par(mfrow=c(2,2))
    Sets=DefOrientSet(mat, 5,180,Multi)
    CorData=array(dim=c(24,dim(mat)[2]/Multi))
    Cols=c(5,6,8,9,12,10)
    colnames(res[[1]])=paste(xlabel,colnames(res[[1]]),sep="-")
    #rownames(res)=rownames(mat)[sels]
    rownames(CorData)=rownames(mat)
    Multi=3 # nombre de colonnes "means"
    DS=SortParm # sort Parm NULL: no sort 1:individual sort array: sort order
    if(!is.null(SortParm) && SortParm>0)
    {
      t=SplitPrintStat(mat,matSD,xlabel,SortParm, Sets,0, DS, FALSE, "Sorting run","","ANG",Multi,DS,LabScale)
      DS=t$SortIdx
    }
    if(mode[1])
    {
      R1=SplitPrintStat(mat,matSD,xlabel,1, Sets,0, DS, FALSE, "Stacking coverage Phe202(226)","% overlap","none",Multi,DS,LabScale)
      CorData[1,]=R1[[8]]
    }
    if(mode[2])
    {
      R13=SplitPrintStat(mat,matSD,xlabel,13, Sets,0, DS, FALSE, "Stacking coverage Phe236(260)","% overlap","none",Multi,DS,LabScale)
      CorData[13,]=R13[[8]]
    }
    if(mode[3])
    {
      R2=SplitPrintStat(mat,matSD,xlabel,2, Sets,4, DS, FALSE, "Angle between CFF and Phe202(226) planes","Angle  ","ANG",Multi,DS,LabScale)  
      CorData[2,]=R2[[8]]
    }
    if(mode[4])
    {
      R14=SplitPrintStat(mat,matSD,xlabel,14, Sets,4, DS, FALSE, "Angle between CFF and Phe236(260) planes","Angle  ","ANG",Multi,DS,LabScale)  
      CorData[14,]=R14[[8]]
    }
    if(mode[5])
    {
      R3=SplitPrintStat(mat,matSD,xlabel,3, Sets,0, DS, FALSE, "CFFcen to Phe202(226) plane distance", "Distance  ( Å )","none",Multi,DS,LabScale) 
      CorData[3,]=R3[[8]]
    }
    if(mode[6])
    {
      R15=SplitPrintStat(mat,matSD,xlabel,15, Sets,0, DS, FALSE, "CFFcen to Phe226(260) plane distance", "Distance  ( Å )","none",Multi,DS,LabScale) 
      CorData[15,]=R15[[8]]
    }
    if(mode[7])
    {
      R4=SplitPrintStat(mat,matSD,xlabel,4, Sets,0, DS, FALSE, "Distance between CFF and Phe202(226)\n centromer projections on Phe202 plane", "Distance  ( Å )","none",Multi,DS,LabScale)
      CorData[4,]=R4[[8]]
    }
    if(mode[8])
    {
      R16=SplitPrintStat(mat,matSD,xlabel,16, Sets,0, DS, FALSE, "Distance between CFF and Phe236(260)\n centromer projections on Phe236 plane", "Distance  ( Å )","none",Multi,DS,LabScale)
      CorData[16,]=R16[[8]]
    }
    if(mode[9])
    {
      R5=SplitPrintStat(mat,matSD,xlabel,5, Sets,2, DS, FALSE, "Angle of CFF and heme planes (A1)", "Angle  ","ANG",Multi,DS,LabScale) 
      CorData[5,]=R5[[8]]
    }
    if(mode[10])
    {
      R7=SplitPrintStat(mat,matSD,xlabel,7, Sets,2, DS, FALSE, "Angle M2-MP1 with M3-M2 (A3)", "Angle  ","ANG",Multi,DS,LabScale) 
      CorData[7,]=R7[[8]]
    }
    if(mode[11])
    {
      R8=SplitPrintStat(mat,matSD,xlabel,8, Sets,2, DS, FALSE, "Angle CFF/heme plane ints.\n with [heme N2-N4] (A0)", "Angle  ","ANG",Multi,DS,LabScale) 
      CorData[8,]=R8[[8]]
    }
    if(mode[12])
    {
      R6=SplitPrintStat(mat,matSD,xlabel,6, Sets,0, DS, FALSE, "Active O/Fe to CFF plane Dist. (D2)", "Distance ( Å )","none",Multi,DS,LabScale) 
      CorData[6,]=R6[[8]]
    }
    if(mode[13])
    {
      R9=SplitPrintStat(mat,matSD,xlabel,9, Sets,0, DS, FALSE, "Dist. of active FE/O  proj. \n and CFFcen (D4) ", "Distance  ( Å )","none",Multi,DS,LabScale) 
      CorData[9,]=R9[[8]]
    }
    if(mode[14])
    {
      R10=SplitPrintStat(mat,matSD,xlabel,10, Sets,0, DS, FALSE, "CFF orient. code (imin)", "","ORI",Multi,DS,LabScale)  #verifier
      CorData[10,]=R10[[8]]
    }
    if(mode[15])
    {
      R11=SplitPrintStat(mat,matSD,xlabel,11, Sets,2, DS, FALSE, "CFF orientation Angle", "Angle ","ANG",Multi,DS,LabScale) #verifier
      CorData[11,]=R11[[8]]
    }
    if(mode[16])
    {
      R12=SplitPrintStat(mat,matSD,xlabel,12, Sets,0, DS, FALSE, "Shorter dist. from CFF oxid. atom to O/Fe","Distance  ( Å )","none",Multi,DS)
      CorData[12,]=R12[[8]]
    }
    if(mode[17])
    {
      R17=SplitPrintStat(mat,matSD,xlabel,17, Sets,0, DS, FALSE, "Dist. cenA-cenB","Distance  ( Å )","none",Multi,DS,LabScale)
      CorData[17,]=R17[[8]]
    }
    if(mode[18])
    {
      R18=SplitPrintStat(mat,matSD,xlabel,18, Sets,0, DS, FALSE, "Angle of heme and CFF1 planes","Angle  ","ANG",Multi,DS,LabScale)
      CorData[18,]=R18[[8]]
    }
    if(mode[19])
    {
      R19=SplitPrintStat(mat,matSD,xlabel,19, Sets,0, DS, FALSE, "Angle of heme and CFF2 planes","Angle  ","ANG",Multi,DS,LabScale)
      CorData[19,]=R19[[8]]
    }
    if(mode[20])
    {
      R20=SplitPrintStat(mat,matSD,xlabel,20, Sets,0, DS, FALSE, "Angle C2>CFF1 cen/cen>Proj. Fe/O","Angle  ","ANG",Multi,DS,LabScale)
      CorData[20,]=R20[[8]]
    }
    if(mode[21])
    {
      R21=SplitPrintStat(mat,matSD,xlabel,21, Sets,0, DS, FALSE, "Angle C2>CFF2 cen/cen>Proj. Fe/O","Angle  ","ANG",Multi,DS,LabScale)
      CorData[21,]=R21[[8]]
    }
    if(mode[22])
    {
      R22=SplitPrintStat(mat,matSD,xlabel,22, Sets,0, DS, FALSE, "Closest CFF1 atom to CFF2 cen","","ORI",Multi,DS,LabScale)
      CorData[22,]=R22[[8]]
    }
    if(mode[23])
    {
      R23=SplitPrintStat(mat,matSD,xlabel,23, Sets,0, DS, FALSE, "Closest CFF2 atom to CFF1 cen","","ORI",Multi,DS,LabScale)
      CorData[23,]=R23[[8]]
    }
    if(mode[24])
    {
      R24=SplitPrintStat(mat,matSD,xlabel,24, Sets,0, DS, FALSE, "Closest CFF* atom to Fe/O","","ORI",Multi,DS,LabScale)
      CorData[24,]=R24[[8]]
    }
    if(mode[25])
    {
      R25=SplitPrintStat(mat,matSD,xlabel,25, Sets,0, DS, FALSE, " "," ","none",Multi,DS,LabScale)
      CorData[25,]=R25[[8]]
    }
    mode=mode[1:nrow(CorData)]
    CorData=CorData[mode==1,]
    rn=rownames(CorData)
    temp=NULL
    for(i in  1:dim(CorData)[1])
    {
      if(sum(CorData[i,])==0) next
      if(is.null(temp))
      {
        temp=CorData[i,]
        temp1=rn[i]
      }else{
        temp=rbind(temp,CorData[i,])
        temp1=c(temp1,rn[i])
      } 
    }
    rownames(temp)=temp1
    CorData=temp
    TCorData=t(CorData)
    Mt=cor(TCorData,method="kendall")
    if(ShowCorrelation)
    {
    par(mfrow=c(1,1))
    test=cor.mtest(TCorData, method="kendall",conf.level = 0.95)
   # cp=corrplot(Mt, p.mat = test$p, sig.level = 0.01,addCoef.col ='white',type = 'lower',number.cex=0.5,insig = "blank") # p-value 0.01 for 95% confidence
    cp=corrplot(Mt, p.mat = test$p, sig.level = 0.01,type = 'lower',number.cex=0.5,insig = "blank",col = COL2('PiYG')) # p-value 0.01 for 95% confidence
     text(cp$corrPos$x, cp$corrPos$y, round(cp$corrPos$corr, 2),cex=0.5,font=2,col="black")
    }
    return(CorData)
  }
  
  # Function to evaluate parameter contribution and explained variance
  evaluate_parameters <- function(data, response, param1, param2, additional_hypothesis = FALSE)
    {
    # Fit a non-linear model with both parameters
    model <- nls(response ~ param1+ param2, data)
    
    # Calculate coefficient of determination (R-squared) for the model
    rsq <- summary(model)$rsq
    
    # Analyze individual parameter importance
    # Option 1: Standardized beta coefficients (relative weights)
    # beta_weights <- coefficient(model) / sd(response)
    
    # Option 2: Profile likelihood analysis - requires 'nlme' package
    # if (require(nlme, quietly = TRUE)) {
    #   profile_lik <- p.l(model, param1, param2)
    #   beta_weights <- profile_lik %>% 
    #     dplyr::mutate(importance = 1 - deviance / summary(model)$deviance)
    # }
    
    # Evaluate explained variance and potential for a third parameter
    if (!additional_hypothesis) {
      return(list(rsq = rsq, beta_weights = beta_weights))
    }
    
    # Fit a model with both parameters and the additional hypothesis
    additional_model <- nls(response ~ f(param1, param2, additional_hypothesis), data)
    
    # Compare R-squared values
    additional_rsq <- summary(additional_model)$rsq
    improvement <- additional_rsq - rsq
    
    # Conclusion
    conclusion <- 
      if (improvement > 0.05) "Significant improvement suggests a potential third parameter" else
        "No significant improvement; current model sufficient"
    
    return(list(rsq = rsq, beta_weights = beta_weights, improvement = improvement,
                conclusion = conclusion))
  }
  

  TruncKin=function(DataSet,SelCFF,Ofset=-1) #Extract subset of data from NormDat.txt file
  # SelCFF: selected CFF
  # Ofset: Number of timesteps to suppress at start of each OMMs (10 to suppress initial state) -1 fro no ofset
  {
    NDL=paste(DataSet,'/NormData.txt',sep="") 
    SIDX=c(1)
    Tab=read.table(NDL)
    Tab=Tab[Tab[,"CFF_Id"]==SelCFF,]
    Doc=Tab[,1:4]
    DOCX=Doc[1,]
    Idx=array(0,dim=c(dim(Tab)[1],1))
    Tab=cbind(Idx,Tab)
    a=Doc[1,1]
    cte=0
    n=0
    Inited=FALSE
    idx=1
    for(i in 1:dim(Tab)[1])
    {
      if(Doc[i,1]!=a)
      {
        n=0
        a=Doc[i,1]
        idx=idx+1
        DOCX=rbind(DOCX,Doc[i,])
        SIDX=c(SIDX,dim(Tabf)[1]+1)
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
    SIDX=c(SIDX,dim(Tabf)[1]+1)
    return(list(Tabf=Tabf,SIDX=SIDX,DOCX=DOCX))
    #Tabf: ed file SIDX: index of start for each OMMM
  }
  
  StatKin=function(TB,Parm1,Parm2,mask=NULL,ShowCorPlot=TRUE) #evaluate correlations between mask items
  {
    Tab1=TB$Tabf
    TLS=TB$SIDX
    x=Tab1[mask,Parm1]
    y=Tab1[mask,Parm2]
    mx=median(x)
    my=median(y)
    Tab2=Tab1[,mask]
    Tab2=Tab2[,apply(Tab2,2,sd)!=0]
    Mt=cor(Tab2,method="pearson")
    test=cor.mtest(Tab2, method="pearson",conf.level = 0.95)$p
    testL=test[test==0]=1e-50
    dd=dim(test)[1]
    re=-array(sapply(test,log10),dim=c(dd,dd))
    if(ShowCorPlot)
    {
      par(mfrow=c(1,1))
      cp=corrplot(Mt, p.mat = test, sig.level = 0.01,addCoef.col ='white',type = 'lower',number.cex=0.5,insig = "blank") # p-value 0.01 for 95% confidence
      text(cp$corrPos$x, cp$corrPos$y, round(cp$corrPos$corr, 2),cex=0.5,font=2,col="black")
      cat("Clic on graph to continue")
      locator(1)
    }
    return (list(medParm1=mx,medParm2=my,Cortast=test,Parm1=Tab1[,Parm1],Parm2=Tab1[,Parm2],TLS=TLS))
  }
  
  valf <- function(RC, D2, D4, A2)
  {
    A2=A2*pi/180
    fu=sqrt(RC*(RC-2*D4*cos(A2))+D4^2+D2^2)
    return(sqrt(RC*(RC-2*D4*cos(A2))+D4^2+D2^2))
  } 
  
  NormDataPlot=function(DataSet,TB,Parm1,Parm2,mask, line=TRUE, ShowCorPlot=TRUE)
  {
  RC=3.233081 # CFF cen to C2 distance
  NDL=paste(DataSet,'/NormData.txt',sep="") 
  rr=StatKin(TB,Parm1,Parm2,mask,ShowCorPlot) #correlation plot
  Border=rr$TLS
  col=c()
  ncol=rainbow(length(Border))
  par(mfrow=c(1,1),mar = c(5, 5, 1, 1))
  for(i in 1:(length(Border)-1))
  {
    for(j in Border[i]:Border[i+1]) col[j]=ncol[i]
  }
  nn=sort(rr$Parm1,decreasing=FALSE,index.return=TRUE)$ix
  x=rr$Parm1[nn]
  y=rr$Parm2[nn]
  c=col[nn]
  plot(x,y, col=c, xlab=paste(Parm1,"( )"),ylab=paste(Parm2,"( )"),cex.axis=2,cex.lab=2) # bipartite parm plot
  if (line) lines(lowess(x,y,f=0.2), col="black",lwd=2)
  #if (line) lines(smooth.spline(x, y), col="black",lwd=2)
  }
  
 t_col= function(color, percent = 50)
   {
   rgb.val=col2rgb(color)
   acol=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,
                alpha = (100 - percent) * 255 / 100)
   return(acol)
 }
 
 FilterNorm=function(DataSet,Ofset,SelCFF,Mask,Obs,Id,Type,FilterTarget)
 {
   Tabx=GetNewTab(DataSet,TRUE,Ofset,SelCFF)
   Tab=Tabx[[1]]
   Tax=Tabx[[2]]
   limits=Tabx[[3]]
   DOCX=Tabx[[4]]
   subId=c()
   for (i in 1:(length(limits)-1))
   {
     if(is.null(FilterTarget) || DOCX$Fe.O[i]==FilterTarget) subId=c(subId,limits[i]:(limits[i+1]-1))
   }
   res=GetData(Tab, Type,Mask,Obs) 
   n=dim(res[[1]])[1]
   Env=res[[1]][subId,]
   Enx=Tax[subId,]
   if(!is.numeric(res[[2]])) Resp=res[[2]][subId,] else Resp=res[[2]][subId]
   
   FullDat=as.matrix(res[[3]])
   return(list(Env=Env,Resp=Resp,Enx=Enx,Select=subId,FullDat=FullDat,limits=limits,DOCX=DOCX))
 }
 
 
 BuildMSDProj=function(Env, filename)
 {
   f=scale(Env)
   g=dist(f)
   Cmd=cmdscale(g, k = 2)
   write.table(Cmd, file=filename)
   return(Cmd)
 }
 
 MDSProj=function(FileName,TB,IminRg,DminRg,Filter=NULL, legP=1,barL="topleft") #Projection MDS reduced coordinate
 {
   Cmd=read.table(FileName,header=TRUE) 
   Tab=TB$FullDat[TB$Select,]
   Tax=TB$Enx
   if(!is.null(Filter))
   {
     Tab=Tab[Tax[,3]==Filter,]
     Cmd=Cmd[Tax[,3]==Filter,]
   }
   rb=c("blue","green","red","darkgrey","purple","cyan","magenta","orange","brown")#,"grey","darkgreen","darkblue","darkred","darkcyan","darkmagenta","darkorange","darkgrey","darkviolet","darkturquoise")
   id=c('C2','C2N4','N4','C8','C3','C3O2','O2','C1','O1')
   col1=rb[1:9][IminRg]
   id1=id[IminRg]
   rb=sapply(rb,t_col,0)
   crit=Tab[,"imin"]
   dmin=Tab[,"dmin"]
   col=rb[crit]
   ThDmin=DminRg[1]
   ThDmax=DminRg[2]
   
  titre=paste("MDS(A0, A1, A3, D2, D4) dmin= ",ThDmin,' to ',ThDmax, ' A range"',sep="")
  sel1=crit%in%IminRg & dmin<=ThDmax & dmin>ThDmin
  plot(Cmd[sel1,1],Cmd[sel1,2],col=col[sel1],xlim=c(-5,5),ylim=c(-3,3),xlab="AX1",ylab="AX2",main=titre,cex.lab=1.5,cex.axis=1.5)
  legend(legP[1],legP[2],legend=as.character(id1),col=col1,pch=20)
  cat("Clic on graph to continue")
  locator(1)
  
  cte=array(0,dim=c(9,7))
  DRange=array(c(2,3,3,3.5,3.5,4,4,5,5,7,7,9,9,15),dim=c(2,7))
  DText=c("2-3 A","3-3.5 A","3.5-4 A","4-5 A","5-7 A","7-9 A","9-15 A")
 
  for(i in 1:7) 
  {
    for(j in 1:9) cte[j,i]=sum(crit==j & dmin<=DRange[2,i] & dmin>DRange[1,i]) 
    Ylim=c(0,max(cte))
  }
  pop=apply(cte,2,sum)
  cten=cte
  for(i in 1:7)
  {
    for(j in 1:9) cten[j,i]= cte[j,i]/pop[i]*100
  }
  cte1=cten[IminRg,]
  cte2=cte[IminRg,]
  rb1=rb[IminRg]
  id1=id[IminRg]
  barplot(cte1,ylim=c(0,100),col=rb1,names.arg=DText,ylab="Imin frame counts",legend.text=as.character(id1),
          args.legend=list(x=barL,bty = "n"),beside=F,main="Fraction (%) of caffeine orientations within dmin range",cex.lab=1.5)
  cat("Clic on graph to continue")
  locator(1)
  barplot(cte2,ylim=c(0,max(pop)),col=rb1,names.arg=DText,ylab="Imin frame counts",legend.text=as.character(id1),
          args.legend=list(x=barL,bty = "n"),beside=F,main="Count of caffeine orientation within dmin range",cex.lab=1.5)
 }
 
 
  
  FilterMDS=function(DataSet,TB,Type,OfSet,ThD=NULL,SelOMM=NULL,IminRange=NULL,ClustNum=0)
  {
    Cmd=read.table(paste(DataSet,"/",Type,sep=""),header=TRUE)
    Tab=TB$Tabf
    TLS=c(TB$SIDX,dim(Tab)[1])
    if(!is.null(SelOMM))
    { 
      Zone=c(TB$SIDX,dim(Cmd)[1])
      Frame=array(FALSE,dim=dim(Cmd)[1])
      range=Zone[SelOMM[1]]:Zone[SelOMM[2]]
      Frame[range]=TRUE
    }else Frame=array(TRUE,dim=dim(Cmd)[1])
    if(is.null(IminRange)) IminRange=c(1,9)
    if(is.null(ThD)) ThD=c(0,100)
    crit=Tab$imin
    dmin=Tab$dmin
    col=rgb(1,0,0,0.30)#"red" #
    sel1=crit%in%IminRange & dmin<=ThD[2] & dmin>ThD[1] & Frame
    plot(Cmd[sel1,1],Cmd[sel1,2],col=col,xlim=c(-5,5),ylim=c(-4,3),xlab="AX1",ylab="AX2",main="MDS coordinates",cex.lab=2,cex.axis=2)
    if(ClustNum==0) nc=1 else nc=ClustNum
    IDC=data.frame()
    resKM=kmeans(Cmd[sel1,],centers=nc,nstart=20) 
    ncc=resKM$centers
    
    for(j in 1:nc)
    {
      Tpt=resKM$cluster==j
      pt=Cmd[sel1,][Tpt,]
      Count=sum(Tpt)
      u=pt-ncc[j,]
      SDD=sqrt(sum(u[,1]^2+u[,2]^2)/dim(u)[1])
      #draw.circle(x=ncc[j,1], y=ncc[j,2], radius=SDD)
      symbols(ncc[j,1],ncc[j,2],circles=SDD,inches=FALSE,add=TRUE) 
      y=Cmd[sel1,]-ncc[j,]
      z=sqrt(y[,1]^2+y[,2]^2)
      pos=which.min(z)
      for(i in 1: (length(TLS)-1))
      {
        if(pos>=TLS[i] & pos<TLS[i+1])
        {
          uid=i
          break
        }
      }
      FTab=Tab[sel1,]
      pp=pos-TLS[uid]+1
      cte=sum(sel1)
      TP=data.frame(uid=uid,Idx=FTab[pos,]$Idx,frame=pp, Count=cte)
      TR=FTab[pos,]
      TM=apply(FTab[Tpt,6:18],2,mean)
      TS=apply(FTab[Tpt,6:18],2,sd)
      if(j==1)IDD=TP else IDD=rbind(IDD,TP)
      if(j==1) ITM=TM else ITM=rbind(ITM,TM)
      if(j==1) ITS=TS else ITS=rbind(ITS,TS)
      if(j==1) ITR=TR else ITR=rbind(ITR,TR)
      if(j==1) CTE=Count else CTE=rbind(CTE,Count)
    }
    return(list(data.frame(uid=uid,frame=pp, Count=cte),FTab[pos,],ITM=ITM,ITS=ITS,ITR=ITR,CTE=CTE))
  }
  
  GetNewTab=function(DataSet,IsNew=TRUE,Ofset=25,SelCFF=1)
  {
    limits=c(1)
    for(i in 1:length(DataSet))
    {
      if(IsNew)
      {
        TB=TruncKin(DataSet[i],SelCFF[i],Ofset) # data extraction
        tabt=TB[[1]][,6:dim(TB[[1]])[2]]
        tabx=TB[[1]][,1:5]
        SI=TB[[2]]
        DX=TB[[3]]
      } else{
        tabt=read.table(DataSet[i], sep=' ',header=TRUE)
        tabt=tabt[,5:dim(tabt)[2]]
        tabx=tabt[,1:4]
        SI=c(1,dim(tabt)[1])
        DX=array(0,dim(4))
      } 
      if(i==1)
      {
        tab=tabt
        tax=tabx
        SIDX=SI
        DOCX=DX
      }else {
        tab=rbind(tab,tabt)
        tax=rbind(tax,tabx)
        DOCX=rbind(DOCX,DX)
        SIDX=c(SIDX,SI)
      }
    }
    return(list(tab,tax,SIDX,DOCX))
  }
  GetData=function(Tab, Type=NULL,Mask=NULL,Obs=1) 
  {
    #Type: 0=raw data, 1=normalized data, 2=angular data for each column of the file
    #Mask: columns to include  in environment
    #Obs: column to predict
    tab=Tab[,c('A0','A1','A2','A3','D2','D4','dmin','imin')]
    dd=dim(tab)
    if(is.null(Type)) Type=rep(0,dd[2])
    if(is.null(Mask)) Mask=rep(1,dd[2])
    ds=sum(Mask==1)
    init=TRUE
    # Build extended parameter set frame
    lc=c()
    lt=c()
    id=1
    for(i in 1:dd[2])
    {
      if(Type[i]==0 || Type[i]==1) 
      {
        if(i %in% Mask) lc=c(lc,id)
        if(i==Obs) lt=id
        id=id+1
      }
      if(Type[i]==2)
      {
        if(i %in% Mask) lc=c(lc,id,id+1)
        if(i==Obs) lt=c(lt,id,id+1)
        id=id+2
      }
    }
    
    nc=sum(Type==0)+sum(Type==1)+2*sum(Type==2)
    tmat=array(0,dim=c(dd[1],nc) ) #extended parameter temp matrix 
    
    index=1 #build extended parameter marix
    nn=c()
    cn=colnames(tab)
    for(i in 1:dd[2])
    {
      if(Type[i]==0) #raw data
      {
        tmat[,index]=tab[,i]
        nn=c(nn,cn[i])
        index=index+1
      }
      if(Type[i]==1) #normalized data not angular
      {
        temp=tab[,i]
        nn=c(nn,cn[i])
        tmat[,index]=(temp-mean(temp))/sd(temp)
        index=index+1
      }
      if(Type[i]==2) #angular data
      {
        temp=tab[,i]
        tmat[,index]=cos(temp*pi/180)
        tmat[,index+1]=sin(temp*pi/180)
        nn=c(nn,paste(cn[i],"_s",sep=""),paste(cn[i],"_c",sep=""))
        index=index+2
      }
    }
    colnames(tmat)=nn
    
    Dmat=tmat[,lc]
    Tmat=tmat[,lt]
    Fmat=tmat
    return(list(Dmat,Tmat,Fmat))
  }
  
 ModNorm=function(ycos,ysin) 
 {
   v1=rep(0,length(ycos))
   for(i in 1:length(ycos))
   {
     mod=sqrt(ycos[i]^2+ysin[i]^2)
     v1[i]=(1-min(abs(mod-1),1)) 
   }
   return (v1)
 }
 
  
  


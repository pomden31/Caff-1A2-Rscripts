AnalyseCFF=function(DataSet,Target, CFFtyp, Mask, FileName,twoCFF,VarListFile,SubVarListFile,IndexFile,TimeSeriesFile,PDBRef,
                    MaskType=NULL,StatMode=NULL,ReadRefPDB=FALSE,DoPrint=TRUE, report=NULL, AltMode=FALSE)
{
  AngRec<--c(0,0,0,0,0,0,0,0,0)
  mode=NULL#"CEN"
  config=2
  subtype=2
  Prefix="data.dcd"
  res=SetUp1(DataSet,Prefix, IndexFile,TimeSeriesFile,PDBRef, ReadRefPDB)
  xyz=res$xyz1
  if(is.null(MaskType))
  {
    temp=as.numeric(unlist(strsplit(Mask,":")))
    MaskData=temp[1]:temp[2]
  }
  else
  {
    sz=dim(xyz)[1]
    if(MaskType==1)  MaskData=1:trunc(0.3*sz)
    if(MaskType==2)  MaskData=trunc(0.3*sz):trunc(0.7*sz)
    if(MaskType==3)  MaskData=trunc(0.7*sz):99999
    if(MaskType==0)  MaskData=1:99999
  }
  if(!ReadRefPDB) MaskData=MaskData[MaskData<=dim(xyz)[1]]
  else
    {
      MaskData=MaskData[MaskData<dim(xyz)[1]]-10
      MaskData=MaskData[MaskData>0]
    }
  pdb=res$CPX1PDB
  parm=res$Parm1
  d1m=array(0,dim=c(7,dim(xyz)[1]))
  data1=data.frame()
  data1T=array()
  data1P=data.frame()
  data2=data.frame()
  data3=data.frame()
  dnorm=data.frame()
  XP=SetNormParm(pdb, parm, CFFtyp,mode, Target)
  XP2=SetNormParm2(pdb, parm, mode, Target)
  pg="-"
  res=NULL
  resp=NULL
  errr=array(0,dim=c(dim(xyz)[1],5))
  AngRec=c(0,0,0,0,0,0,0,0,0)
  sz=dim(xyz)[1]
  for(i in 1:sz)
  {
    cat(".")
    if(i%%50==0)cat("\n")
    if(AngleCorMode==0) AngCor=0
    if(AngleCorMode==1){if(i==1) AngCor=1 else AngCor=2}
    if(AngleCorMode==2){if(NewBatch && i==1) AngCor=1 else AngCor=2}
    res=NormIntersect(i,xyz, parm, CFFtyp,mode, Target,XP,res,AngCor)
    errr[i,]=unlist(res$err)
    AngRec=AngRec+res$AngRec
    data1=rbind(data1,res$df)
    data1T=rbind(data1T,res$Idl)
    resp=NormIntersect2(i,xyz, parm,mode, Target,XP2,AngCor)
    data1P=rbind(data1P,resp)
    #M2L=res$DP$M2
    Way=res$Way
    if(ActivateNorm>0)
    {
      ndata=SetNormData(i,xyz, Target,XP,Way)
      ndata=cbind(ndata,resp[c("IDENA","IDENB","IDENC")])
      ndata=cbind(ndata,res$df[c("dmin","imin","imin1")])
      dnorm=rbind(dnorm,ndata) 
    }
  }
  print(errr[1:20,])
  AngRec=AngRec/dim(xyz)[1]
  if(ActivateNorm>0)
  {
    mm=c("A3")
    for(y in 1 :length(mm))
    {
     tp=dnorm[,mm[y]]
     dnorm[,mm[y]]=sapply(tp,NormTo180)
    }
  }
  
  data1T=data1T[2:dim(data1T)[1],]
  XP3=SetNormParm3(pdb, parm, CFFtyp)
  for(i in 1:dim(xyz)[1])
  {
    cat("*")
    if(AngleCorMode==0) AngCor1=0
    if(AngleCorMode==1){if(i==1) AngCor1=1 else AngCor1=2}
    if(AngleCorMode==2){if(NewBatch && i==1) AngCor1=1 else AngCor1=2}
    if(i%%50==0)cat("\n")
    sta=StackAnalysis(pdb, i,xyz, parm, CFFtyp,NULL,XP3,AngCor1)
    data2=rbind(data2,sta[[1]])
    data3=rbind(data3,sta[[2]])
  }
  
  if(!is.null(ForceTarget) )Target1=ForceTarget else Target1=Target
  if(Target=="FE")
  {
    M1="FE to CFFx plane distance (D2)"
    M2="Shorter dist. from \n C1238 CFF atoms to FE"
  }else
  {
    M1="OXY to CFFx plane distance (D2)"
    M2="Shorter distance from \n C1238 CFF atoms to OXY"
  }
  RR=PrintSet4(data1P,xyz,MaskData,mode,DoPrint)
  data1P=RR[[1]]
  res4=RR[[2]]
  
  rv=PrintSet1(data1,data2,data3,xyz,MaskData,mode,DoPrint)
  data2=rv[[1]][[1]]
  data3=rv[[1]][[2]]
  res1=rv[[2]]
  res1b=rv[[3]]
 
  PS=PrintSet2(data1,M1,MaskData,mode,DoPrint)
  data1=as.data.frame(PS[[1]])
  res2=PS[[2]]
  
  PS3=PrintSet3(MaskData,data1,data2,M2,data1T,mode,DoPrint)
  stat=PS3$PS$stat
  data1=PS3$PS$data1
  res3=PS3$PSH
  
  SRM=data.frame(SCM=res1[[1]],CPPA=res1[[2]],CCTPP=res1[[3]],DPCPC=res1[[4]],
                 CHPA=res2[[1]],OCPD=res2[[2]],ACHAC=res2[[3]],ANG=res2[[4]],
                 DPCCO=res3[[1]],COC=res3[[2]],ANGC2=res3[[3]], DCATO=res3[[4]],
                 SCMb=res1b[[1]],CPPAb=res1b[[2]],CCTPPb=res1b[[3]],DPCPCb=res1b[[4]]
                 ,CTC=data1$Angc2[1])
  SRM1=data.frame()
  ahcb=res4[[3]]
  ahca=res4[[2]]

  SRM1=data.frame(DCAB=res4[[1]],AHCA=res4[[2]],AHCB=res4[[3]],VCENA=res4[[4]],VCENB=res4[[5]], IDENA=res4[[6]], IDENB=res4[[7]],IDENC=res4[[8]],AHAB=res4[[9]],CTC=data1$Angc2[1])

  T1=read.table(VarListFile,header=TRUE)
  T2=read.table(file=SubVarListFile,header=TRUE)
  rr=t(array(c("Type",T1$Code[T1$SubRes!=0]),dim=c(sum(T1$SubRes!=0)+1,1)))
  PrintLen=sum(T2$Print==1)
  rr1=array(c(T2$Legend[T2$Print==1]),dim=c(PrintLen,1))
  lres=as.data.frame(T2$Legend)
  for(i in 1:dim(T1)[1])
  {
    if(T1$SubRes[i]==0) next
    j=T1$Res[i]
    jl=0
    if(j=="PS1-2") jl=res1[[T1$SubRes[i]]]
    if(j=="PS1-3") jl=res1b[[T1$SubRes[i]]]
    if(j=="PS2") jl=res2[[T1$SubRes[i]]]
    if(j=="PS3") jl=res3[[T1$SubRes[i]]]
    if(j=="PS4") jl=res4[[T1$SubRes[i]]]
    rr2=array(0,dim=c(sum(T2$Print==1),1))
    ind=1
    for(k in 1:35)
    {
      if(T2$Print[k]==0) next
      if(length(dim(jl))==2) rr2[ind]=jl[T2$RS[k],T2$Sub[k]]
      else rr2[ind]=0
      ind=ind+1
    }
    rr1=cbind(rr1,rr2 )
  }
  rr=rbind(rr,rr1)
  if(!is.null(report)) write.table(rr, file=report)
  #if(!is.null(report)) write(knitr::kable(rr, format = "markdown"), file=report)
  NewBatch<<-FALSE
  return (list(SRM,data1,data2,data3,data1P,SRM1,rr,dnorm, AngRec=AngRec)) # SRM see matrix definition
}


  SetNormParm=function(Apdb, Parm, CFFTyp=1,mode="CEN", Target="FE")
  {
    if(!is.null(ForceTarget) )Target=ForceTarget
    Na="N1"
    Nb="N2"
    Nc="N3"
    Nd="N4"
    RHM=Parm$HEME
    if(CFFTyp==1)
    {
      RCF=Parm$CFF
    }else
      RCF=Parm$CFF1
    C3="C3"
    C2="C2"
    C1="C1"
    C8="C8"
    N4="N4"
    C7="C7"
    C3="C3"
    C4="C4"
    N4="N4"
    O2="O2"
    O1="O1"
    AtSG=AtomSel(Apdb,434,"SG")
    AtNA=AtomSel(Apdb,RHM,Na)
    AtNB=AtomSel(Apdb,RHM,Nb)
    AtNC=AtomSel(Apdb,RHM,Nc)
    AtND=AtomSel(Apdb,RHM,Nd)
    AtC1=AtomSel(Apdb,RCF,C1)
    AtC4=AtomSel(Apdb,RCF,C4)
    AtC2=AtomSel(Apdb,RCF,C2)
    AtN4=AtomSel(Apdb,RCF,N4)
    AtC8=AtomSel(Apdb,RCF,C8)
    AtC3=AtomSel(Apdb,RCF,C3)
    AtC7=AtomSel(Apdb,RCF,C7)
    AtN4=AtomSel(Apdb,RCF,N4)
    AtO2=AtomSel(Apdb,RCF,O2)
    AtO1=AtomSel(Apdb,RCF,O1)
    RFE=Parm$Fe
    AtFE=AtomSel(Apdb,RFE,"FE")
    RFO=Parm$OXY
    AtO=AtomSel(Apdb,RFO,"O1")

    return (list(AtC1=AtC1,AtC2=AtC2,AtC3=AtC3,AtC4=AtC4,AtC7=AtC7,AtC8=AtC8,AtN4=AtN4,AtO2=AtO2,AtO1=AtO1,AtFE=AtFE,AtO=AtFE,
                       AtNA=AtNA,AtNB=AtNB,AtNC=AtNC,AtND=AtND,AtSG=AtSG ))
  }


  SetAngle=function(V1,V2,OldNorm=NULL, LastAng=NULL)
  {
    Ang=acos(V1%*%V2)*180/pi
    Norm=cross(V1,V2)
    if(!is.null(LastAng)) Sgn=sign(LastAng)
    else Sgn=1
    if(!is.null(OldNorm))
    {
      Way=Norm%*%OldNorm
      if(sign(way<0)) Sgn=-Sgn
    }
    else Sgn=1
    Ang=Sgn*Ang
    return(list(Ang=Ang,Norm=Norm))
  }

#' Title
#'
#' @param NewAng
#' @param LastAng
#'
#' @return
#' @export
#'
#' @examples
 AdjustAngle=function(NewAng,LastAng=NULL)
 {
   if(!is.null(LastAng))
   {
     if(abs(NewAngle+360-LastAngle) < abs(NewAngle-LastAngle)) NewAngle=NewAngle+360
     if(abs(NewAngle-360-LastAngle) < abs(NewAngle-LastAngle)) NewAngle=NewAngle-360
   }
   return(NewAng)
 }


  NormIntersect=function(index,xyz, Parm, CFFTyp=1,mode="CEN", Target="FE",AT,res0, Init=0) #1: set 2:use
  {
    if(!is.null(ForceTarget) )Target=ForceTarget
    calibrate=FALSE
    c1=xyz[index,AT$AtC1$xyz]
    c4=xyz[index,AT$AtC4$xyz]
    c2=xyz[index,AT$AtC2$xyz]
    n4=xyz[index,AT$AtN4$xyz]
    c8=xyz[index,AT$AtC8$xyz]
    c3=xyz[index,AT$AtC3$xyz]
    c7=xyz[index,AT$AtC7$xyz]
    na=xyz[index,AT$AtNA$xyz]
    nc=xyz[index,AT$AtNC$xyz]
    nb=xyz[index,AT$AtNB$xyz]
    nd=xyz[index,AT$AtND$xyz]
    sg=xyz[index,AT$AtSG$xyz]
    fe=xyz[index,AT$AtFE$xyz]
    o1=xyz[index,AT$AtO1$xyz]
    o2=xyz[index,AT$AtO2$xyz]
    tt=as.vector(unlist(GCenter(c1,c2,c3,TRUE)$center))
    if(calibrate) AngRec=CalibrateCFFAngle(c1,o1,c2,n4,c8,c3,o2) else AngRec=c(0,0,0,0,0,0,0,0,0)
    DHO=1.66 #distance Fe-O
    RCFF=3.233081
    if(Target=="FE") mode1=2
    if(Target=="OXY") mode1=3
    NP1=SetNP1(na,nb,nc,nd)$NP1
    F2=GetM1(na,nb,nc,nd, fe,NP1,DHO, mode1)
    M1=F2$Parm$M1
    O1=F2$Parm$O1

    F0=GetM3(c1,c2,c3,1) #mode: 1:center of circle, 2:mean
    M3=F0$cen #CFF cen
    RCFF=F0$RCFF #distance CFF center of C1,C2,C3 circle to C1 (C2,C3)
    TP=SetNP1(na,nb,nc,nd)
    NP1=TP$NP1
    NP2=SetNP2(c1,c2,c3)$NP2
    if(NP2%*%(nc-na)/norm_vec(na-nc)>0)
    {
      Face='L'
    }else{
      Face='R'
      NP2=-NP2
    } 
    dM1P=abs(NP2%*%(M1-M3))
    M21=M1-NP2*dM1P
    M22=M1+NP2*dM1P
    #M2 determination
    if(abs(NP2%*%(M21-M3))<abs(NP2%*%(M22-M3)))
    {
      M2=M21
      Way1=1
    }else{
      M2=M22
      Way1=-1
    } 
    
    F3=GetSysAngleDist(M1,M2,M3,c2,NP1,NP2,nb,nd,Way1)
    MP1=F3$MP1
    LI=F3$LI
    nbd=(nd-nb)/norm_vec(nb-nd)
    dMP1M1=norm_vec(MP1-M1)
    dMP1M2=norm_vec(MP1-M2)
    #MP1 validation
    if(abs(NP2%*%(MP1-M3))>1e-10 || abs(NP1%*%(MP1-M1))>1e-10
       || abs(abs((MP1-M1)%*%(MP1-M2)/dMP1M1/dMP1M2)-abs(NP1%*%NP2))>1e-5)
    {
      stop("Error UtilB 320")
    }
    
    A0=F3$Parm$A0 # Angle of N4->N2 vector with intersect of CFF and heme plane
    A1=F3$Parm$A1 # #Angle of CFF and heme planes
    A2=F3$Parm$A2 # Angle of C2->cen with M2->cen. M2 (Proj. Fe/Oxy on CFF plane)
    A3=F3$Parm$A3 # Angle M2->Mp1 with M2->cen. MP1 (Proj. of Fe/Oxy on intersect of Heme and CFF planes)
    A4=F3$Parm$A4 # Angle M2->Mp1 with cen->Mp1. MP1 (Proj. of Fe/Oxy on intersect of Heme and CFF planes)
   
    #formally sign of way indicate if M2 is on the same side of heme than M3
    # sign of A1 indicate if CFF is on face 'L' or 'R'
    
    if(Target=="FE") P=fe
    if(Target=="OXY") P=O1
    F4=IDCffProject(c1,c2,c3,c8, P,Init) #Fe/Oxy proj. on CFF plane used as reference (see target)
    dmin=F4$dmin # Distance from P to closest CFF Oxyd. atom (P: Fe or Oxy see target)
    imin1=F4$imin1 # Id of CFF  Oxyd. atom (P: Fe or Oxy see target) closest to P: Fe or Oxy see target
    imin=F4$imin # orientation code index
    Idl=F4$ASTAT # orientation matrix
    
    if(Target=="FE") mode=2
    if(Target=="OXY") mode=3
    D4=norm_vec(M2-M3)
    D2=Way1*norm_vec(M2-M1)
   
   
    CCFC=GetBaseDescFromCoord(fe,na,nb,nc,nd,D2,D4,A0,A1,A2,A3,Way1,mode, Face)
    m1=CCFC$Coord$M1
    m2=CCFC$Coord$M2
    m3=CCFC$Coord$M3
    m4=CCFC$Coord$M4
    mp1=CCFC$Coord$MP1
    np2=CCFC$Coord$NP2
   
    err=max(norm_vec(m1-M1),norm_vec(m2-M2),norm_vec(m3-M3),norm_vec(m4-c2))
    DP=data.frame(M1=M1,M2=M2,M3=M3,M4=c2,MP1=MP1)
    df=data.frame(AHC=A1,DFP=D2, AAC=A3,ADB=A4,DPHPC=D4,dmin=dmin,imin=imin,imin1=imin1, Ang=A0,Dist=dmin, Angc2=A2)
    errf=data.frame(err,err1=norm_vec(m1-M1),err2=norm_vec(m2-M2),err3=norm_vec(m3-M3),err4=norm_vec(m4-c2))
    if(err>0.2)
    {
      stop("Invalid err UtilB 373")
    }
    return (list(Way=Way1,df=df,Idl=Idl,DP=DP,err=errf))
  } 
 

  DefAng=function(RefAtom, Center, Rayon, CFA,UseProj=TRUE)
  {
    CFA=t(CFA)
    MatRC=array(c(CFA[,1],(CFA[,1]+CFA[,2])/2,CFA[,2],(CFA[,2]+CFA[,3])/2,CFA[,3],(CFA[,3]+CFA[,4])/2,
    CFA[,4],(CFA[,4]+CFA[,5])/2, CFA[,5],(CFA[,5]+CFA[,6])/2, CFA[,6],(CFA[,6]+CFA[,7])/2,CFA[,7],(CFA[,7]+CFA[,1])/2), dim=c(3,14))
    dmin=1000
    imin=0
    imin1=0
    for (i in 1:14)
    {
      if(UseProj)proj1=GProject(Center,Rayon, MatRC[,i])
      else proj1=MatRC[,i]
      temp=norm_vec(proj1-RefAtom)
      if(temp<dmin)
      {
        dmin=temp
        imin=i
      }
    }
    id=c("C1","C1C4","C4","C4C2","C2","C2N4","N4","N4C8","C8","C8C3","C3","C3C7","C7","C7C1")
    return (list(imin=imin,dmin=dmin,id=id[imin]))
  }


  DefAng3=function(RefAtom, Center, Rayon, CFA,UseProj=TRUE)
  {
    CFA=t(CFA)
    MatRC=array(c(CFA[,1],(CFA[,1]+CFA[,2])/2,CFA[,2],(CFA[,2]+CFA[,3])/2,CFA[,3],(CFA[,3]+CFA[,4])/2,
                  CFA[,4],(CFA[,4]+CFA[,5])/2, CFA[,5],(CFA[,5]+CFA[,6])/2, CFA[,6],(CFA[,6]+CFA[,7])/2,CFA[,7],(CFA[,7]+CFA[,1])/2), dim=c(3,14))
    dmin=1000
    imin=0
    for (i in c(1,5,9,11))
    {
      if(UseProj)proj1=GProject(Center,Rayon, MatRC[,i])
      else proj1=MatRC[,i]
      temp=norm_vec(proj1-RefAtom)
      if(temp<dmin)
      {
        dmin=temp
        imin=i
      }
    }
    id=c("C1","C1C4","C4","C4C2","C2","C2N4","N4","N4C8","C8","C8C3","C3","C3C7","C7","C7C1")
    return (list(imin=imin,dmin=dmin,id=id[imin]))
  }

  DefAng2=function(RefAtom, Center, Rayon, CFA)
  {
    CFA=t(CFA)
    MatRC=array(c(CFA[,1],CFA[,2], CFA[,3], CFA[,4]), dim=c(3,4))
    dmin=1000
    imin=0
    imin1=0
    for (i in 1:4)
    {
      proj1=GProject(Center,Rayon, MatRC[,i])
      temp=norm_vec(proj1-RefAtom)
      if(temp<dmin)
      {
        dmin=temp
        imin=i
      }
    }
    id=c("C1","C2","C8","C3")
    return (list(imin=imin,dmin=dmin,id=id[imin]))
  }

  SetNormParm2=function(Apdb, Parm,mode="CEN", Target="FE")
  {
    if(!is.null(ForceTarget) )Target=ForceTarget
    Na="N1"
    Nb="N2"
    Nc="N3"
    Nd="N4"
    RHM=Parm$HEME
    RCF=Parm$CFF
    RCF1=Parm$CFF1
    AtSG=AtomSel(Apdb,434,"SG")
    C3="C3"
    C2="C2"
    C1="C1"
    C8="C8"
    N4="N4"
    C7="C7"
    C4="C4"
    AtNA=AtomSel(Apdb,RHM,Na)
    AtNB=AtomSel(Apdb,RHM,Nb)
    AtNC=AtomSel(Apdb,RHM,Nc)
    AtND=AtomSel(Apdb,RHM,Nd)

    if(RCF!=0)
    {
      AtC1a=AtomSel(Apdb,RCF,C1)
      AtC4a=AtomSel(Apdb,RCF,C4)
      AtC2a=AtomSel(Apdb,RCF,C2)
      AtN4a=AtomSel(Apdb,RCF,N4)
      AtC8a=AtomSel(Apdb,RCF,C8)
      AtC3a=AtomSel(Apdb,RCF,C3)
      AtC7a=AtomSel(Apdb,RCF,C7)
    } else{
      AtC1a=0
      AtC4a=0
      AtC2a=0
      AtN4a=0
      AtC8a=0
      AtC3a=0
      AtC7a=0
    }
    if(RCF1!=0)
    {
      AtC1b=AtomSel(Apdb,RCF1,C1)
      AtC4b=AtomSel(Apdb,RCF1,C4)
      AtC2b=AtomSel(Apdb,RCF1,C2)
      AtN4b=AtomSel(Apdb,RCF1,N4)
      AtC8b=AtomSel(Apdb,RCF1,C8)
      AtC3b=AtomSel(Apdb,RCF1,C3)
      AtC7b=AtomSel(Apdb,RCF1,C7)
    }else{
      AtC1b=0
      AtC4b=0
      AtC2b=0
      AtN4b=0
      AtC8b=0
      AtC3b=0
      AtC7b=0
    }
    RFE=Parm$Fe
    RFO=Parm$OXY
    AtO=AtomSel(Apdb,RFO,"O1")
    AtFE=AtomSel(Apdb,RFE,"FE")
    return (list(AtC1a=AtC1a,AtC2a=AtC2a,AtC3a=AtC3a,AtC4a=AtC4a,AtC7a=AtC7a,AtC8a=AtC8a,AtN4a=AtN4a,
                 AtC1b=AtC1b,AtC2b=AtC2b,AtC3b=AtC3b,AtC4b=AtC4b,AtC7b=AtC7b,AtC8b=AtC8b,AtN4b=AtN4b
                       ,AtFE=AtFE,AtNA=AtNA,AtNB=AtNB,AtNC=AtNC,AtND=AtND,AtSG=AtSG))
  }
 
  GetParmab=function(Target,fe,na,nb,nc,nd,c1,c2,c3)
  {
    if(!is.null(ForceTarget) )Target=ForceTarget
    DHO=1.66 #distance Fe-O
    if(Target=="FE") mode1=2
    if(Target=="OXY") mode1=3
    NP1=SetNP1(na,nb,nc,nd)$NP1
    F2=GetM1(na,nb,nc,nd, fe,NP1,DHO, mode1)
    M1=F2$Parm$M1
    M3=GetM3(c1,c2,c3,1)$cen
    NP2=SetNP2(c1,c2,c3)$NP2
    dM1P=abs(NP2%*%(M1-M3))
    M21=M1-NP2*dM1P
    M22=M1+NP2*dM1P
    if(abs(NP2%*%(M21-M3))<abs(NP2%*%(M22-M3)))
    {
      M2=M21
      Way=1
    }else
    {
      M2=M22
      Way=-1
    }
    NP2=NP2*Way
    AngRec[19,]=SetLI(NP1,NP2,nd,nb)
    AA=GetAngle(NP1,NP2,19,AngRec)
    if(AA<0) AA=180+AA
    return (AA)
  }
 
  NormIntersect2=function(index,xyz, Parm,mode="CEN", Target="FE",AT, Init=0) # 1init 2:use
  {
    if(!is.null(ForceTarget) )Target=ForceTarget
    DHO=1.66 #FE to OXY distance
    RHM=Parm$HEME
    RCF=Parm$CFF
    RCF1=Parm$CFF1
    fe=xyz[index,AT$AtFE$xyz]
    na=xyz[index,AT$AtNA$xyz]
    nb=xyz[index,AT$AtNB$xyz]
    nc=xyz[index,AT$AtNC$xyz]
    nd=xyz[index,AT$AtND$xyz]
    NP1=SetNP1(na,nb,nc,nd)$NP1
    if(Target=="FE") M1=GetM1(na,nb,nc,nd, fe,NP1,DHO, 2)$Parm$M1
    if(Target=="OXY") M1=GetM1(na,nb,nc,nd, fe,NP1,DHO,3)$Parm$M1
    
    if(RCF!=0)
    {
    c1a=xyz[index,AT$AtC1a$xyz]
    c2a=xyz[index,AT$AtC2a$xyz]
    c8a=xyz[index,AT$AtC8a$xyz]
    c3a=xyz[index,AT$AtC3a$xyz]
    NP2a=SetNP2(c1a,c2a,c3a)$NP2
    M3a=GetM3(c1a,c2a,c3a,1)$cen
    M2a=dot(M3a-M1, NP2a)*NP2a+M1
    }
    if(RCF1!=0)
    {
    c1b=xyz[index,AT$AtC1b$xyz]
    c2b=xyz[index,AT$AtC2b$xyz]
    c8b=xyz[index,AT$AtC8b$xyz]
    c3b=xyz[index,AT$AtC3b$xyz]
    NP2b=SetNP2(c1b,c2b,c3b)$NP2
    M3b=GetM3(c1b,c2b,c3b,1)$cen
    M2b=dot(M3b-M1, NP2b)*NP2b+M1 
    }
    if(RCF!=0 && RCF1!=0) DCAB=norm_vec(M3a-M3b) else DCAB=0  # Distance des centromeres caffeine A et B ****
    
    if(RCF!=0)#effecteur
    {
    AHCa=GetParmab(Target,fe,na,nb,nc,nd,c1a,c2a,c3a)
    ResHa=IDCffProject(c1a,c2a,c3a,c8a,M1,Init)
    anga=ResHa$Ang
    ResH=ResHa
    }else{
      AHCa=0
      anga=0
      IDENA=0
      IDENB=0
    }

    if(RCF1!=0)#producteur
    {
      AHCb=GetParmab(Target,fe,na,nb,nc,nd,c1b,c2b,c3b)
      ResHb=IDCffProject(c1b,c2b,c3b,c8b,M1,Init)
      angb=ResHb$Ang
      IDENC=ResHb$imin
      ResH=ResHb
    }else{
      AHCb=0
      angb=0
      IDENA=0
      IDENB=0
      IDENC=0
    }

    if(RCF!=0 && RCF1!=0)
    {
    VNA=cross(c1a-c2a,c1a-c3a)
    VNB=cross(c1b-c2b,c1b-c3b)
    
    if(Init==0)AHab=GetAngle(VNA,VNB,0)
    if(Init==1) AHab=GetAngle(VNA,VNB,-11,AngRec,TRUE)
    if(Init==2) AHab=GetAngle(VNA,VNB,11,AngRec,TRUE)
    
    ResH=ResHb
    ResCa=IDCffProject(c1a,c2a,c3a,c8a,M3b,Init)
    ResCb=IDCffProject(c1b,c2b,c3b,c8b,M3a,Init)
    IDENA=ResCa$imin
    IDENB=ResCb$imin
    }else{
      AHab=0
    }
    IDENC=ResH$imin
    if(RCF!=0 || RCF1!=0)Indc=ResH$imin
    return (data.frame(DCAB=DCAB,AHCA=AHCa,AHCB=AHCb,AHAB=AHab,VCENA=anga,VCENB=angb, IDENA=IDENA, IDENB=IDENB,IDENC=Indc))

    # DCAB Distance des centromeres caffeine A et B
    #AHAB angle plan cafb plan cafa
    #AHCA angle plan heme plan CFFa
    #AHCB angle plan heme et CFFb
    #VCENA angle (cena-Proj Fe/plan CFFa) avec (cena-c2a CFF)
    #VCENB angle (cena-Proj Fe/plan CFFa) avec (cena-c2a CFF)
    #IDENA index du residue de la CFFa le plus proche du cen de la CFF b
    #IDENB index du residue de la CFFb le plus proche du cen de la CFF a
    #IDENC index du residue de la CFF2 le plus proche du O/FE
}


  SetNormParm3=function(Apdb, Parm, CFFTyp=1, target=NULL)
  {
    TyrRes=202
    TyrRes1=236
    if(CFFTyp==1)
    {
      RCF=Parm$CFF
    }else
      RCF=Parm$CFF1
    AtSG=AtomSel(Apdb,434,"SG")
    if(is.null(target))
    {
      AtC1=AtomSel(Apdb,RCF,"C1")
      AtC2=AtomSel(Apdb,RCF,"C2")
      AtC3=AtomSel(Apdb,RCF,"C3")
      AtP1=AtomSel(Apdb,TyrRes,"CD1")
      AtP2=AtomSel(Apdb,TyrRes,"CD2")
      AtP3=AtomSel(Apdb,TyrRes,"CZ")
      AtP1b=AtomSel(Apdb,TyrRes1,"CD1")
      AtP2b=AtomSel(Apdb,TyrRes1,"CD2")
      AtP3b=AtomSel(Apdb,TyrRes1,"CZ")
    }
    else{
      AtC1=AtomSel(Apdb,target$p1,target$c11)
      AtC2=AtomSel(Apdb,target$p1,target$c21)
      AtC3=AtomSel(Apdb,target$p1,target$c31)
      AtP1=AtomSel(Apdb,target$p2,target$c12)
      AtP2=AtomSel(Apdb,target$p2,target$c22)
      AtP3=AtomSel(Apdb,target$p2,target$c32)
      AtP1b=0
      AtP2b=0
      AtP3b=0
    }
   return(list(AtC1=AtC1,AtC2=AtC2,AtC3=AtC3,AtP1=AtP1,AtP2=AtP2,AtP3=AtP3
                     ,AtP1b=AtP1b,AtP2b=AtP2b,AtP3b=AtP3b,AtSG=AtSG))
  }



  StackAnalysis=function(Apdb, index,xyz, Parm, CFFTyp=1, target=NULL,AT,AngCor)
  {
    c1=xyz[index,AT$AtC1$xyz]
    c2=xyz[index,AT$AtC2$xyz]
    c3=xyz[index,AT$AtC3$xyz]
    p1=xyz[index,AT$AtP1$xyz]
    p2=xyz[index,AT$AtP2$xyz]
    p3=xyz[index,AT$AtP3$xyz]
    p1b=xyz[index,AT$AtP1b$xyz]
    p2b=xyz[index,AT$AtP2b$xyz]
    p3b=xyz[index,AT$AtP3b$xyz]

    GC_CFF=GCenter(c1,c2,c3,TRUE)
    GC_PHE=GCenter(p1,p2,p3,TRUE)
    GC_PHEb=GCenter(p1b,p2b,p3b,TRUE)

    Vdir1=cross(unlist(c1-c2),unlist(c3-c2))
    Vdir1=Vdir1/norm_vec(Vdir1)
    Vdir2=cross(unlist(p1-p2),unlist(p3-p2))
    Vdir2b=cross(unlist(p1b-p2b),unlist(p3b-p2b))
    Vdir2=Vdir2/norm_vec(Vdir2)
    Vdir2b=Vdir2b/norm_vec(Vdir2b)
    cv1=as.numeric(GC_CFF$center)
    cv2=as.numeric(GC_PHE$center)
    cv2b=as.numeric(GC_PHEb$center)

    res=SurfInterf(cv1,GC_CFF$Rayon,Vdir1,cv2,GC_PHE$Rayon,Vdir2, 0.1,AngCor, 1)
    resb=SurfInterf(cv1,GC_CFF$Rayon,Vdir1,cv2b,GC_PHEb$Rayon,Vdir2b, 0.1, AngCor,2)
    return (list(res,resb))
  }


  SurfInterf=function(Cen1,Ray1,vdir1,Cen2,Ray2,vdir2, step=0.1, AngCor,Typ)
  {
    vdir1=vdir1/norm_vec(vdir1)
    vdir2=vdir2/norm_vec(vdir2)
    proj12=Cen1+vdir2%*%(Cen2-Cen1)*vdir2 #coordonn?e projection ortho Cent_CFF sur plan Phe
    proj21=Cen2+vdir1%*%(Cen1-Cen2)*vdir1 #coordonn?e projection ortho Cent_Phe sur plan CFF
    VP=cross(vdir1,vdir2) # // vect droite intersect plan Phe plan CFF
    NVP=norm_vec(VP)
    if(NVP!=0)
    {
      VP=VP/NVP
      VN=cross(VP, vdir2) # normal droite dans plan Phe (VN,VP forme une base locale du plan phe)
      VN=VN/norm_vec(VN)
      Tmat=(matrix(c(VP,VN, vdir2), nrow=3,ncol=3)) #VN,VP,Vdir2 forment une base locale du plan phe)
      Tmat1=solve(Tmat)
      CC1=Tmat1 %*% t(t(proj12-Cen2)) #Proj Cent_phe_cffp en coordonn?es locales phe
      ry=Ray1*(vdir1%*%vdir2) # rayons de l'ellipse projet?e (rx, ry)
    }else{
      CC1=Cen1-Cen2
      CC1[3]=0
      ry=Ray1 # rayons de l'ellipse projet?e (rx, ry)
    }

    CC2=c(0,0,0)  #Cent2 sur plan2 en coordonn?es locales P2
    rx=Ray1
    rx2=rx*rx
    ry2=ry*ry
    ct=0
    cte=0

    for(i in seq(-Ray2,+Ray2, step))
    {

      for(j in seq(-Ray2,Ray2, step))
      {
        if(i*i+j*j>Ray2*Ray2) next #not in c2
        cte=cte+1
        CXX=i-CC1[1]
        CYY=j-CC1[2]
        if(CXX*CXX/rx2+CYY*CYY /ry2>1) next
        ct=ct+1
      }
    }

    dist=norm_vec(Cen1-proj12) #distance Cent_cff to plan phe
    if(length(vdir2)!=0) dist2=norm_vec(CC1)
    

      if(AngCor==0) angl=GetAngle(vdir1,vdir2,0) #angle Plan1 Plan2
      if(AngCor==1) angl=GetAngle(vdir1,vdir2,-12,AngRec,TRUE)
      if(AngCor==2) angl=GetAngle(vdir1,vdir2,12,AngRec,TRUE)

    
    Rec=ct/cte*100 # %recouvrement project cercle 1 sur plan2 with cercle 2
    AngRec<-- AngRec
    return (data.frame(dist=dist,dist2=dist2, angl=angl,rec=Rec))
  }


  #return coordonnées xx,xyxz du point du cercle dand l'axe vect(Pt-center)

  GProject=function(Center,Rayon, Pt){
    V=(Pt-Center)/norm_vec(as.matrix(Pt-Center))
    res=as.matrix(Center)+V*Rayon
    return (data.frame(x=res[1],y=res[2],z=res[3]))
  }


  GRProject=function(C1,C2,C3, At){
  Center=unlist(GCenter(C1,C2,C3,TRUE)$center)
  catAvar1=c("C2","C2N4","N4","C8","C3", "C3O2","O2","O2C1","C1","C1C4","C4","C2C2")
  VN=cross((C2-Center),(C3-Center))
  VN=VN/norm_vec(VN)
  resp=array(0,dim=c(3,12))
  for(i in 1:12)
  {
  a=30*pi/180*(i-1)
  c=cos(a)
  s=sin(a)
  ux=VN[1]
  uy=VN[2]
  uz=VN[3]
  RotMat=array(0,dim=c(3,3))
  RotMat[1,1]=ux*ux*(1-c)+c
  RotMat[1,2]=ux*uy*(1-c)-uz*s
  RotMat[1,3]=ux*uz*(1-c)+uy*s
  RotMat[2,1]=ux*uy*(1-c)+uz*s
  RotMat[2,2]=uy*uy*(1-c)+c
  RotMat[2,3]=uy*uz*(1-c)-ux*s
  RotMat[3,1]=ux*uz*(1-c)-uy*s
  RotMat[3,2]=uy*uz*(1-c)+ux*s
  RotMat[3,3]=uz*uz*(1-c)+c
  X=RotMat%*%(C2-Center)+Center
  resp[,i]=c(X[1],X[2],X[3])
  }
  imin=0
  dmin=1000
  for(i in 1:12)
  {
    d=norm_vec(At-resp[,i])
    if(d<dmin)
    {
      dmin=d
      imin=i
    }
  }
  imin1=imin
  return (list(resp=resp,imin=imin,dmin=dmin,Id=catAvar1[imin1]))
  }


  NullsetcafAng=function()
  {
    C2=c(57.727, 50.313, 31.249)
    C5=c(55.699, 48.834, 30.815)
    N4=c(55.797, 48.855, 29.459)
    C8=c(54.772, 48.126, 29.044)
    C3=c(53.012, 46.559, 30.012)
    C6=c(54.780, 47.935, 31.202)
    O2=c(53.723, 47.028, 33.095)
    C1=c(55.140, 48.524, 34.874)
    O1=c(56.893, 50.137, 33.878)
    C2N4=(C2+N4)/2
    C3O2=(C3+O2)/2
    GC=GCenter(C1,C2,C3,TRUE)
    Center=as.vector(unlist(GC$center))
    Rayon=unlist(GC$Rayon)
    VN=cross((C2-Center),(C3-Center))
    VN=VN/norm_vec(VN)
    ss=sign(VN%*%(cross(C2-Center,C3O2-Center)))
    AC2=acos((C2-Center)%*%(C2-Center)/norm_vec(C2-Center)/norm_vec(C2-Center))*180/pi
    AC24=acos((C2N4-Center)%*%(C2-Center)/norm_vec(C2N4-Center)/norm_vec(C2-Center))*180/pi
    AN4=acos((N4-Center)%*%(C2-Center)/norm_vec(N4-Center)/norm_vec(C2-Center))*180/pi
    AC8=acos((C8-Center)%*%(C2-Center)/norm_vec(C8-Center)/norm_vec(C2-Center))*180/pi
    AC3=acos((C3-Center)%*%(C2-Center)/norm_vec(C3-Center)/norm_vec(C2-Center))*180/pi
    AC32=acos((C3O2-Center)%*%(C2-Center)/norm_vec(C3O2-Center)/norm_vec(C2-Center))*180/pi
    AO2=360-acos((O2-Center)%*%(C2-Center)/norm_vec(O2-Center)/norm_vec(C2-Center))*180/pi
    AC1=360-acos((C1-Center)%*%(C2-Center)/norm_vec(C1-Center)/norm_vec(C2-Center))*180/pi
    AO1=360-acos((O1-Center)%*%(C2-Center)/norm_vec(O1-Center)/norm_vec(C2-Center))*180/pi
    Id=c("C2","C2N4","N4","C8","C3","C3O2","O2","C1","O1")
  }



  GProject=function(Center,Rayon, Pt){
    V=(Pt-Center)/norm_vec(as.matrix(Pt-Center))
    res=as.matrix(Center)+V*Rayon
    return (data.frame(x=res[1],y=res[2],z=res[3]))
  }

  PrintRes=function(RES, id)
  {
    rss1=data.frame()
    rss2=data.frame()
    Res=RES[[1]]
    Res1=RES[[6]]
    CTC=Res[1,113]

    for(i in c(0:15))
    {
      df=data.frame(Res[1,1+7*i], Res[2,1+7*i],Res[3,1+7*i], Res[4,1+7*i],Res[5,1+7*i])
      df1=data.frame(Res[1,2+7*i], Res[2,2+7*i],Res[3,2+7*i], Res[4,2+7*i],Res[5,2+7*i])
      df=cbind(df,df1)
      df2=data.frame(Res[1,3+7*i], Res[2,3+7*i],Res[3,3+7*i], Res[4,3+7*i],Res[5,3+7*i])
      df=cbind(df,df2)
      df3=data.frame(Res[1,4+7*i], Res[2,4+7*i],Res[3,4+7*i], Res[4,4+7*i],Res[5,4+7*i])
      df=cbind(df,df3)
      df4=data.frame(Res[1,5+7*i], Res[2,5+7*i],Res[3,5+7*i], Res[4,5+7*i],Res[5,5+7*i])
      df=cbind(df,df4)
      df5=data.frame(Res[1,6+7*i], Res[2,6+7*i],Res[3,6+7*i], Res[4,6+7*i],Res[5,6+7*i])
      df=cbind(df,df5)
      df6=data.frame(Res[1,7+7*i], Res[2,7+7*i],Res[3,7+7*i], Res[4,7+7*i],Res[5,7+7*i])
      df=cbind(df,df6)
      if(i==0) rss1=df
      else rss1=rbind(rss1,df)
    }
    Res=Res1
    for(i in c(0:8))
    {
      df=data.frame(Res[1,1+7*i], Res[2,1+7*i],Res[3,1+7*i], Res[4,1+7*i],Res[5,1+7*i])
      df1=data.frame(Res[1,2+7*i], Res[2,2+7*i],Res[3,2+7*i], Res[4,2+7*i],Res[5,2+7*i])
      df=cbind(df,df1)
      df2=data.frame(Res[1,3+7*i], Res[2,3+7*i],Res[3,3+7*i], Res[4,3+7*i],Res[5,3+7*i])
      df=cbind(df,df2)
      df3=data.frame(Res[1,4+7*i], Res[2,4+7*i],Res[3,4+7*i], Res[4,4+7*i],Res[5,4+7*i])
      df=cbind(df,df3)
      df4=data.frame(Res[1,5+7*i], Res[2,5+7*i],Res[3,5+7*i], Res[4,5+7*i],Res[5,5+7*i])
      df=cbind(df,df4)
      df5=data.frame(Res[1,6+7*i], Res[2,6+7*i],Res[3,6+7*i], Res[4,6+7*i],Res[5,6+7*i])
      df=cbind(df,df5)
      df6=data.frame(Res[1,7+7*i], Res[2,7+7*i],Res[3,7+7*i], Res[4,7+7*i],Res[5,7+7*i])
      df=cbind(df,df6)
      if(i==0) rss2=df
      else rss2=rbind(rss2,df)
    }
    RR=rbind(rss1,rss2)
    CTC=data.frame(array(CTC,dim=c(1,dim(RR)[2])))
    colnames(CTC)=colnames(RR)
    RR=rbind(RR,CTC)
    Rname=c("SCM1","CPPA1","CCTPP1","DPCPC1","CHPA","OCPD","ACHAC","ANG","DPCCO","COC","ANGC2","DCMO","SCM2","CPPA2","CCTPP2","DPCPC2",
            "DCAB","AHCA","AHCB","VCENA","VCENB","IDENA","IDENB","IDENC","AHAB","CTC")
    Cname=c("MaxE","SdEx","DifQtl","GCen11","StartVal", "GSd11","DifExFt", "NormTest","Sum1","Sum2","GCen12","GCen22","GCen13","GCen23","GCen33",
            "GSd12","GSd22","GSd13","GSd23","GSd33","MaxE12","MaxE22","MaxE13","MaxE23","MaxE33","SCC12","SCC22","SCC13","SCC23","SCC33",
            "CNT12","CNT22","CNT13","CNT23","CNT33")
    rownames(RR)=Rname
    colnames(RR)=Cname
    return (RR)
  }


  PrintSet1=function(data1,data2,data3,xyz,Mask,mode,DoPrint)
  {
    #par(mfrow=c(1,2), mar=c(2.5,5,3,0.5))
    par(mfrow=c(1,2),mar = c(2, 4, 2, 3))
    f1=sum(data2$rec==100)/dim(xyz)[1]*100
    f2=sum(data2$rec>80)/dim(xyz)[1]*100
    f3=sum(data2$rec>50)/dim(xyz)[1]*100
    f1b=sum(data3$rec==100)/dim(xyz)[1]*100
    f2b=sum(data3$rec>80)/dim(xyz)[1]*100
    f3b=sum(data3$rec>50)/dim(xyz)[1]*100
    col1=rep("grey",length(data2$rec))
    col2="blue"
    if(ReadRefPDB)
    {
      col1[1:10]=rep("red",10)
    }
    text=c(paste("Frequency=100: ",sprintf("%.2f",f1)),paste("Frequency >80 : ",sprintf("%.2f",f2)),paste("Frequency >50 : ",sprintf("%.2f",f3)))
    
    MT="Stacking coverage CFFx/Phe1(226)"
    plot(data2$rec,ylim=c(min(50,min(data2$rec)),100),col=col1,ylab="REC1 surf. overl. (%)",main=MT,cex.main=0.8)
    lines(lowess(data2$rec,f=0.02), col=col2,lwd=2)
    a=PlotFreq(data2$rec,Mask, mtext=MT, mode,DoPrint)

    MT="Stacking coverage (%) CFFx/Phe2(260)"
    plot(data3$rec,ylim=c(min(50,min(data3$rec)),100),col=col1,ylab="REC2 Surf. overl. (%)",main=MT,cex.main=0.8)
    lines(lowess(data3$rec,f=0.02), col=col2,lwd=2)
    #legend(300, 40,cex=0.75,title="Overlapping surface  (%)", legend=text)
    a1=PlotFreq(data3$rec,Mask, mtext=MT, mode,DoPrint)

    MT="Angle CFFx | Phe1(226) planes"
    data2$angl=CentAngl(NormAngl(data2$angl))
    ymin=min(data2$angl)-10
    ymax=max(data2$angl)+10
    plot(data2$angl,ylim=c(ymin,ymax),col=col1, main =MT,ylab="ANGD1 Ang.",cex.main=0.8)
    lines(lowess(data2$angl,f=0.02), col=col2,lwd=2)
    b=PlotFreq(data2$angl,Mask, mtext=MT, mode,DoPrint)

    MT="Angle CFFx | Phe2(260) planes"
    data3$angl=CentAngl(NormAngl(data3$angl))
    ymin=min(data3$angl)-10
    ymax=max(data3$angl)+10
    plot(data3$angl,ylim=c(ymin,ymax),col=col1, main =MT,ylab="ANGD2 Ang.",cex.main=0.8)
    lines(lowess(data3$angl,f=0.02), col=col2,lwd=2)
    b1=PlotFreq(data3$angl,Mask, mtext=MT, mode,DoPrint)
    
    MT="CENx to Phe1(226) plane dist. (A)"
    x=1:dim(xyz)[1]
    y=data2$dist
    l=lowess(y,f=0.1)
    ymin=min(y)-5
    ymax=max(y)+5
    plot(y,xlim=c(0,dim(xyz)[1]),ylim=c(0,ymax),xlab="Frames",ylab="DISP1 dist.(Å)",col=col1, main =MT,cex.main=0.8)
    lines(l, col=col2,lwd=2)
    c=PlotFreq(data2$dist,Mask, mtext=MT, mode,DoPrint)
    
    MT="Dist. CENx proj/Phe1\n to CEN_Phe1 (226)"
    x=1:dim(xyz)[1]
    y=data2$dist2
    l=lowess(y,f=0.1)
    ymin=min(y)-5
    ymax=max(y)+5
    plot(y,xlim=c(0,dim(xyz)[1]),ylim=c(0,ymax),xlab="Frames",ylab="DISJ1 dist.(A)",col=col1, main =MT,cex.main=0.8)
    lines(l, col=col2,lwd=2)
    d=PlotFreq(data2$dist2,Mask, mtext=MT, mode,DoPrint)
    
   
    MT="CENx to Phe2(260) plane dist.(A)"
    x=1:dim(xyz)[1]
    y=data3$dist
    l=lowess(y,f=0.1)
    ymin=min(y)-5
    ymax=max(y)+5
    plot(y,xlim=c(0,dim(xyz)[1]),ylim=c(0,ymax),col=col1,xlab="Frames",ylab="DISP2 dist.(A)", main =MT,cex.main=0.8)
    lines(l, col=col2,lwd=2)
    c1=PlotFreq(data3$dist,Mask, mtext=MT, mode,DoPrint)

    MT="Dist. CENx proj/Phe2\n to CEN_Phe2 (260)"
    x=1:dim(xyz)[1]
    y=data3$dist2
    l=lowess(y,f=0.1)
    ymin=min(y)-5
    ymax=max(y)+5
    plot(y,xlim=c(0,dim(xyz)[1]),ylim=c(0,ymax),col=col1,xlab="Index", main =MT, ylab="DISJ2 dist.(A)",cex.main=0.8)
    lines(l, col=col2,lwd=2)
    d1=PlotFreq(data3$dist2,Mask,mtext=MT, mode,DoPrint)
    
    return(list(DT23=list(data2,data3), AD1=list(REC=a,ANGL=b,DIST=c,DIST2=d),AD2=list(REC=a1,ANGL=b1,DIST=c1,DIST2=d1)))
    
  }
  #-----------------

  PrintClass=function(amat,MainLabel,YLabel,col1,col2)
  {
    data=amat
    dataf=c(data[1])
    codemask=rep(0,length(data))
    for(i in 2:length(data))
    {
      u=data[i]
      if(abs(data[i]+9-dataf[i-1])<abs(data[i]-dataf[i-1]))
        {
        u=data[i]+9
        codemask[i]=9
        }
      if(abs(data[i]-9-dataf[i-1])<abs(data[i]-dataf[i-1]))
        {
        u=data[i]-9
        codemask[i]=-9
        }
      dataf=c(dataf,u)
    }
    ymin=min(dataf)-1
    ymax=max(dataf)+1
    code=c("C2","C2N4","N4","C8","C3","C3O2","O2","C1","O1")
    code=c("C2(C12N3)","C2N4","N4(N9)","C8(C8)","C3(C14N7)","C3O2","O2(O13)","C1(C10N1)","O1(O11)")

    labels=c()
    for( i in ymin:ymax) labels=c(labels,code[(i+8)%%9+1])
    plot(dataf,ylim=c(ymin,ymax),col=col1, main =MainLabel,ylab=YLabel,cex.main=1,axes=FALSE)
    u=10.0/length(dataf)
    ff=max(0.02,u)
    lines(lowess(dataf,f=ff), col=col2,lwd=2)
    axis(2,las = 2,at=ymin:ymax,labels=labels, cex.axis=0.6)
    axis(1)
    return (list(dataf, code, codemask))
  }


  PrintSet4=function(data1P,xyz,Mask,mode,DoPrint)
  {
    data1P_bk=data1P
    par(mfrow=c(1,2))
    col1=rep("grey",length(data1P$DCAB))
    col2="blue"
    if(ReadRefPDB)
    {
      col1[1:10]=rep("red",10)
    }
    
    a=b=c=d=e=f=g=h=i=array(0,dim=c(5,7))

    data=data1P$DCAB
    if(!sum(data)==0)
    {
    ymin=min(data)-10
    ymax=max(data)+10
    MainText="DCAB:Dist.CENa-CENb (A)"
    plot(data,ylim=c(ymin,ymax),col=col1, main =MainText,ylab="DCAB Dist.(A)",cex.main=1)
    lines(lowess(data,f=0.02), col=col2,lwd=2)
    a=PlotFreq(data1P$DCAB,Mask, mtext=MainText, mode,DoPrint)
    
    }
   
    data1P$AHCA=CentAngl(data1P$AHCA)
    data=data1P$AHCA
    if(!sum(data)==0)
    {
    ymin=min(data)-10
    ymax=max(data)+10
    MainText="Angle of heme and CFFa planes (A1a)"
    plot(data,ylim=c(ymin,ymax),col=col1, main =MainText,ylab="Ang.A1a [AHCA]",cex.main=1)
    lines(lowess(data,f=0.02), col=col2,lwd=2)
    b=PlotFreq(data1P$AHCA,Mask, mtext=MainText, mode,DoPrint)
    }
   
    
    data1P$AHCB=CentAngl(data1P$AHCB)
    data=data1P$AHCB
    if(!sum(data)==0)
    {
    ymin=min(data)-10
    ymax=max(data)+10
    MainText="Angle of heme and CFFb planes (A1b)"
    plot(data,ylim=c(ymin,ymax),col=col1, main =MainText,ylab="Angle A1b  [AHCB]",cex.main=1)
    lines(lowess(data,f=0.02), col=col2,lwd=2)
    c=PlotFreq(data1P$AHCB,Mask, mtext=MainText, mode,DoPrint)
    }
    
    
    data1P$AHAB=CentAngl(data1P$AHAB)
    data=data1P$AHAB
    if(!sum(data)==0)
    {
    ymin=min(data)-10
    ymax=max(data)+10
    MainText="Angle of CFFa and CFFb planes"
    plot(data,ylim=c(ymin,ymax),col=col1, main =MainText,ylab="Angle [AHAB]",cex.main=1)
    lines(lowess(data,f=0.02), col=col2,lwd=2)
    i=PlotFreq(data1P$AHAB,Mask, mtext=MainText, mode,DoPrint)
    }
    
  
    data1P$VCENA=CentAngl(NormAngl(data1P$VCENA))
    data=data1P$VCENA
    if(!sum(data)==0)
    {
    ymin=min(data)-10
    ymax=max(data)+10
    MainText="Angle C2-CENa | CENa-M2 (A2a)"
    plot(data,ylim=c(ymin,ymax),col=col1, main =MainText,ylab="Angle A2a [VCENA]",cex.main=1)
    lines(lowess(data,f=0.02), col=col2,lwd=2)
    d=PlotFreq(data1P$VCENA,Mask, mtext=MainText, mode,DoPrint)
    }
   
    
    data1P$VCENB=CentAngl(NormAngl(data1P$VCENB))
    data=data1P$VCENB
    if(!sum(data)==0)
    {
    ymin=min(data)-10
    ymax=max(data)+10
    MainText="Angle C2-CENb | CENb-M2 (A2b)"
    plot(data,ylim=c(ymin,ymax),col=col1, main =MainText,ylab="Angle A2b [VCENB]",cex.main=1)
    lines(lowess(data,f=0.02), col=col2,lwd=2)
    e=PlotFreq(data1P$VCENB,Mask, mtext=MainText, mode,DoPrint)
    }
    
    
    #par(mfrow=c(1,2), mar=c(2.5,4.5,3,0.5))
    if(!sum(data1P$IDENA)==0)
    {
      val=data1P$IDENA
      val=((val-1)%%9)+1
      MainText="CFFa C2 orient.|CENaCENb (IDENA)"
      re=PrintClass(val, MainText,"",col1,col2)
      f=PlotFreq(val,Mask, mtext= MainText, mode,FALSE,FALSE,AltMode)
    }
    
    if(!sum(data1P$IDENB)==0)
    {
      val=data1P$IDENB
      val=((val-1)%%9)+1
      MainText="CFFb C2 orient.|CENbCENa (IDENB)"
      re=PrintClass(val,MainText,"",col1,col2)
      g=PlotFreq(val,Mask, mtext="Closest CFFb atom to CENa", mode,FALSE,FALSE,AltMode)
    }
    
    if(!sum(data1P$IDENC)==0)
    {
      val=data1P$IDENC
      val=((val-1)%%9)+1
      MainText="CFFb C2 orient.|FE/OXY proj. (IDENC)"
      re=PrintClass(val,MainText,"",col1,col2)
      h=PlotFreq(val,Mask, mtext=MainText, mode,FALSE,FALSE, AltMode)
    }
    
    return (list(D1P=data1P, PS4H=list(DCAB=a,AHCA=b,AHCB=c,VCENA=d,VCENB=e,IDENA=f,IDENB=g,IDENC=h,AHAB=i)))
  }



  
  #----------------------

  PrintSet2=function(data1,M1,Mask,mode,DoPrint)
  {
    par(mfrow=c(1,2),mar = c(2, 4, 2, 3))
    col1=rep("grey",length(data1$AHC))
    col2="blue"
    if(ReadRefPDB)
    {
      col1[1:10]=rep("red",10)
    }
    
    # #Angle of CFF and heme planes
    MT="Angle CFFx and heme planes (A1)"
    data1$AHC=CentAngl(data1$AHC)
    avar=data1$AHC
    mx=10*ceiling(ceiling(1.2*max(avar)+1)/10)
    plot(avar, ylim=c(min(avar)-10,1.2*max(avar)), col=col1,ylab="Ang.A1 [ACH] ",main=MT,cex.main=1)
    lines(lowess(avar,f=0.2)$y, col=col2)
    a= PlotFreq(data1$AHC,Mask, mtext=MT, mode,DoPrint)
    
    
    #distance iron -CFF plane
    MT="Dist.D2:FE/OXY-CFFx plane (A)"
    avar=data1$DFP
    plot(avar, ylim=c(min(avar)-10,max(avar)+10), col=col1,ylab=MT,xlab="MD frames",main=M1,cex.main=1)
    lines(lowess(avar,f=0.2)$y, col=col2)
    b= PlotFreq(data1$DFP,Mask, mtext=M1, mode,DoPrint)

    # Angle M2->Mp1 with M2->cen. MP1 (Proj. of Fe/Oxy on intersect of Heme and CFF planes)
    MT="Angle M2-MP1 | M2-CENx (A3)"
    #data1$AAC=NormAngl(NoNegAngl(data1$AAC))
    data1$AAC=NormTo180(data1$AAC)
    avar=data1$AAC
    plot(avar, ylim=c(min(min(avar)-10,min(1.2*avar)),1.2*max(avar)), col=col1,ylab="Ang.A3 [AAC] ",main= MT,cex.main=1)
    lines(lowess(avar,f=0.1)$y, col=col2)
    #mx=10*ceiling(ceiling(1.2*max(avar)+1)/10), yaxp=c(0,mx,10)
    c= PlotFreq(data1$AAC,Mask, mtext= MT, mode,DoPrint)
    
    # Angle of N4->N2 vector with intersect of CFF and heme plane
    MT="Angle N4-N2 axis with \n inters. CFFx-heme planes (A0)"
    data1$Ang=CentAngl(NormAngl(data1$Ang))
    avar=data1$Ang
    mx=10*ceiling(ceiling(1.2*max(avar)+1)/10)
    plot(avar, ylim=c(min(min(avar)-10,min(1.2*avar)),1.2*max(avar)), col=col1,xlab="MD frames",ylab="Ang.A0 [ANG] ",main=MT,cex.main=1)
    lines(lowess(avar,f=0.2)$y, col=col2) #, yaxp=c(0,mx,10),
    d= PlotFreq(data1$Ang,Mask,  mtext=MT, mode,DoPrint)
    
    # Angle M2->Mp1 with cen->Mp1. MP1 (Proj. of Fe/Oxy on intersect of Heme and CFF planes)
    MT="Angle M2-MP1 | CENx-MP1 (A4)"
   # data1$ADB=NormAngl(NoNegAngl(data1$ADB))
    data1$ADB=NormTo180(data1$ADB)
    avar=data1$ADB
    plot(avar, ylim=c(min(avar)-10,1.2*max(avar)), col=col1,xlab="MD frames",ylab="Ang.A4 [ADB] ",main=MT,cex.main=1)
    lines(lowess(avar,f=0.1)$y, col=col2)#mx=10*ceiling(ceiling(1.2*max(avar)+1)/10), yaxp=c(0,mx,10)
    e= PlotFreq(data1$ADB,Mask, mtext=MT, mode,DoPrint)
    
    return (list(D1=list(data1),DH=list(AHC=a,DFP=b,AAC=c,ANG=d,ADB=e)))
  }

  #----------

 
  #----------------------

  PrintSet3=function(mask,data1,data2,M2,data1T,mode,DoPrint)
  {
    par(mfrow=c(1,2),mar = c(2, 4, 2, 3))
    col1=rep("grey",length(data1$DPHPC))
    col2="blue"
    if(ReadRefPDB)
    {
      col1[1:10]=rep("red",10)
    }

    #Distance from CFF cen to proj. Fe/Oxy on CFF plane (see target)
    MT="Dist.D4: CENx-FE/OXY proj.(A)"
    avar=data1$DPHPC
    plot(avar, ylim=c(0,1.2*max(avar)), col=col1,main=MT,ylab="D4 dist.[DPHPC]",cex.main=1)
    lines(lowess(avar,f=0.1)$y, col=col2)
    a=PlotFreq(data1$DPHPC,mask, mtext=MT, mode,DoPrint)


    avar1=data1$imin1
    avar=data1$imin
    Dist1=data1$Dist
    cc=c(8,1,5,4)
    A=cc[avar1]
    avar1=A
    
    par(mfrow=c(1,1), mar=c(2.5,5,3,0.5))
    MT="CFFx orientation code"
    re=PrintClass(data1$imin,MT,"IMIN index",col1,col2) #Orientation
    b=PlotFreq(data1$imin,mask, mtext=MT, mode,FALSE,FALSE)
    
    par(mfrow=c(1,1))
    catAvar=c("C1","C2","C3","C8")
    catAvar1=c("C2(C12N3)","C2N4","N4(N9)","C8(C8)","C3(C14N7)","C3O2","O2(O13)","C1(C10N1)","O1(O11)")
    stat=array(0,dim=c(3,9))
    stat2=array(0,dim=c(1,4))
    for(i in 1:length(avar)){
      stat[1,avar[i]]=stat[1,avar[i]]+1
      stat[2,avar1[i]]=stat[2,avar1[i]]+1
      if(data1$Dist[i]<4) stat[3,avar1[i]]=stat[3,avar1[i]]+1
      if(i<=4) stat2[1,]=stat2[1,]+data1T[i,]
    }
    stat2=100*stat2/length(avar)
    colnames(stat)=catAvar1
    colnames(stat2)=catAvar
    rowname=c("Orientation","Closer_CFF_CH3","Less_Than_4A")
   # try
    {
    barplot(height=stat,xlab="Caffeine atoms",ylab="Frequency", ylim=c(0,1.5*max(stat)), beside=T,col=c("black","red","green"),
            legend.text = rowname,args.legend = list(x = "topright",inset = c( 0, 0)),main="CFFx orientation and expected metabolite",cex.main=1
            ,name.arg=catAvar1,cex.names=0.7)
   }
    par(mfrow=c(1,2),mar = c(2, 4, 2, 3))
     # Angle of C2->cen with M2->cen. M2 (Proj. Fe/Oxy on CFF plane)
    MT="Ang. C2-CENx | M2-CENx (A2)"
    data1$Angc2=CentAngl(NormAngl(data1$Angc2))
    data=data1$Angc2
    plot(data, ylim=c(min(0,1.2*min(data)),max(data)+20), col=col1,xlab="MD frames",ylab="Ang.A2 [ANGC2] ", main=MT,cex.main=1)
    lines(lowess(data,f=0.01)$y, col=col2)
    c=PlotFreq(data1$Angc2,mask, mtext=MT, mode,DoPrint)

    # Distance from P to closest CFF Oxyd. atom (P: Fe or Oxy see target)
    MT="Dist. FE/OXY-C1238 CFF atom"
    avar=data1$Dist
    plot(avar, ylim=c(0,1.2*max(avar)), col=col1,xlab="MD frames",ylab="Distance (A)", main=MT,cex.main=1)
    lines(lowess(avar,f=0.01)$y, col=col2)
    d=PlotFreq(data1$Dist,mask,mtext=MT, mode,DoPrint)
    
    return (list(PS=list(stat=stat,data1=data1),PSH=list(DPHPC=a,IMIN=b,ANGC2=c,DIST=d)))
  }
  #-----------

  SetUp2=function(Id1,Id2,prefix,IndexFile,TimeSeriesFile){
    table=read.table(IndexFile, sep='\t')
    colnames(table)=table[1,]
    table=table[2:dim(table)[1],]
    Parm1=ExtractParms(Id1, table)
    Parm2=ExtractParms(Id2 , table)
    DataSet1=c(paste("Set",Id1, sep='_'),table[Id1,2])
    DataSet2=c(paste("Set",Id2, sep='_'),table[Id2,2])
    ff=ReadFile(prefix,DataSet1, DataSet2,Parm1,Parm2,"P450d")
    CPX1PDB=ff$CPX1PDB
    CPX2PDB=ff$CPX2PDB
    xyz1=ff$xyz1
    xyz2=ff$xyz2
    TS=ReadTimeSeries(TimeSeriesFile)
    X1=GetFrameTime(DataSet1[1], DataSet1[2],FALSE, TS)
    X2=GetFrameTime(DataSet2[1], DataSet2[2],FALSE, TS)
    return (list(config=config, Parm1=Parm1,Parm2=Parm2,DataSet1=DataSet1,DataSet2=DataSet2, CPX1PDB=CPX1PDB,CPX2PDB=CPX2PDB, X1=X1,X2=X2,xyz1=xyz1,xyz2=xyz2))
  }

  SetUp1=function(Id1,prefix,IndexFile,TimeSeriesFile,PDBRef,ReadRefPDB=FALSE){
    table=read.table(IndexFile, sep='\t')
    colnames(table)=table[1,]
    table=table[2:dim(table)[1],]
    Parm=ExtractParms(Id1, table)
    DataSet=c(paste("Set",Id1, sep='_'),table[Id1,2])
    ff=ReadFile1(prefix,DataSet,Parm,"P450d",PDBRef, ReadRefPDB)
    APDB=ff$CPX1PDB
    xyz=ff$xyz1
    TS=ReadTimeSeries(TimeSeriesFile)
    X=GetFrameTime(DataSet[1], DataSet[2],FALSE, TS)
    return (list(Parm1=Parm,DataSet1=DataSet, CPX1PDB=APDB,X1=X,xyz1=xyz))
  }

  ExtractParms=function(Id, table){
    parm=list(P450d=0,FMNd=0,HEME=0,Fe=0,OXY=0,CytLig=0,CFF=0,FMN=0,CFF1=0, FMNT='',HEMET='', frame=0)
    code=table[Id,3]
    cc=c('P','R','H','E','O','F','C','c','X')
    cci=c(1,2,3,4,5,6,7,8,9)
    Mat=array('', dim=9)
    Ofset=0
    for (i in 1:nchar(code)){
      t=cci[cc==substr(code,i,i)]
      if (t==1)
      {
        parm$P450d=table[Id,5]:table[Id,6]
        Ofset=as.numeric(table[Id,6])
        parm$CytLig=434
      }
      if (t==2)
      {
        parm$FMNd=table[Id,7]:table[Id,8]
        Ofset=as.numeric(table[Id,8])
      }
      if (t>2) Ofset=Ofset+1
      if (t==3) parm$HEME=Ofset
      if (t==4) parm$Fe=Ofset
      if (t==5) parm$OXY=Ofset
      if (t==6) parm$FMN=Ofset
      if (t==7) parm$CFF=Ofset
      if (t==8) parm$CFF1=Ofset
    }
    parm$FMNT=table[Id,9]
    parm$HEMET=table[Id,10]
    parm$frame=as.numeric(table[Id,4])
    return (parm)
  }


  ReadFile1=function(Prefix,ASet1, Parm1,RefDomain, PDBRef,StartPDB=FALSE){
    CPX1DCD <- read.dcd(paste(ASet1[2],"/",Prefix,sep=""))
    CPX1PDB  <- read.pdb(paste(ASet1[2],"/data.pdb", sep=""))
    if(RefDomain=="FMNd") RefDomain1=Parm1$FMNd
    if(RefDomain=="P450d")RefDomain1=Parm1$P450d
    if (StartPDB)
    {
      RefCPX=read.pdb(paste(ASet1[2],"/RefPDB.pdb", sep=""))
      xyz0=RefCPX$xyz
      xyz1=ReadTraj(CPX1PDB,CPX1DCD,RefDomain1,PDBRef)
      d=dim(xyz0)[2]
      u=matrix(xyz0[1,],ncol=d)
      xyz2=u
      for(i in 1:9)
      {
        xyz2=rbind(xyz2,u)
      }
      if(file.exists(paste(ASet1[2],"/RefPDB_M.pdb", sep="")))
      {
        RefCPX1=read.pdb(paste(ASet1[2],"/RefPDB_M.pdb", sep=""))
        xyz01=RefCPX1$xyz
        d=dim(xyz01)[2]
        u=matrix(xyz01[1,],ncol=d)
        for(i in 1:9)xyz2=rbind(xyz2,u)
      }
      xyz2=rbind(xyz2,xyz1)
      xyz2=as.xyz(xyz2)
    }else{
     xyz2=ReadTraj(CPX1PDB,CPX1DCD,RefDomain1,PDBRef)
    }
    return (list(CPX1PDB=CPX1PDB,xyz1=xyz2))
  }


  ReadTimeSeries=function(TimeSeriesFile){
    table=read.table(TimeSeriesFile, sep='\t')
    t1=table[,1]
    set=table[grep("%",t1),1]
    set=substr(set,start=2,stop=nchar(set[]))
    Nframe=as.integer(table[grep("%",t1),3])
    Nframe=Nframe[2:length(Nframe)]
    RState=table[grep("%",t1),2]
    RState=RState[2:length(RState)]
    NL=(1:dim(table)[1])[grep("%",t1)]
    NL=NL[NL[]>1]+1
    Frame=data.frame()
    for (i in 1:(length(set)-1))
    {
      fid=1
      sps=0
      for(j in 1:Nframe[i])
      {
        id=NL[i]+j-1
        ad=data.frame(set[i+1],table[id,2],t1[id],RState[i],table[id,2],fid,table[id,3], sps , table[id,4])
        Frame=rbind(Frame,ad)
        fid=fid+as.integer(table[id,2])
        sps=sps+as.integer(table[id,2])*as.numeric(table[id,3])
      }
    }
    names(Frame)[]=c("SetId","FrameId", "Type", "Redox", "NFrame","StartF", "TStep", "Spsec", "Valid")
    return(Frame)
  }


  GetFrameTime=function(SetName, Redox,IncludeEqu, FB)
  {
    F0=(FB[,9]=="1")
    F1=(FB[,1]==SetName)
    F2=(tolower(FB[,4])==tolower(Redox))
    F4=(FB[,3]=="Prod")
    F5=(FB[,3]=="Equi")
    if( ! IncludeEqu) FS= F0 & F1 & F2 & F4 else
      FS=F0 & F1 & F2 & (F4 | F5)
    SubTab=FB[FS,]
    sz=sum(as.integer(SubTab[,5]))
    tstep=as.numeric(SubTab[1,7])
    Tmat=array(0, dim=sz)
    Tmat[1]=tstep
    FBloc=1
    ns=as.integer(SubTab[1,5])-1
    for(i in 2:sz)
    {
      if(ns<=0)
      {
        FBloc=FBloc+1
        tstep=as.integer(SubTab[FBloc,7])
        ns=as.integer(SubTab[FBloc,5])
      }
      Tmat[i]=Tmat[i-1]+tstep
      ns=ns-1
    }
    return (Tmat)
  }


  ReadTraj=function(Apdb,Adcd, Adomain,RefPDB){
    ca1.inds <- atom.select(Apdb, elety="CA",verbose=FALSE)
    ca2.inds <- atom.select(Apdb, resno=Adomain,verbose=FALSE)
    cal.inds=combine.select(ca1.inds,ca2.inds)
    xyz <- fit.xyz(fixed=RefPDB$xyz, mobile=Adcd,
                   fixed.inds=cal.inds$xyz,
                   mobile.inds=cal.inds$xyz)
    dim(xyz) = dim(Adcd)
    return (xyz)
  }


  AtomSel=function(apdb,resnum,Atom){
    res.inds= atom.select(apdb, resno=resnum,verbose=FALSE)
    ato.inds=atom.select(apdb, elety=Atom,verbose=FALSE)
    res.inds=combine.select(res.inds,ato.inds)
    return (res.inds)
  }

  # retourne le centre du cercle circonscrit de 3 points
  GCenter=function(P1,P2,P3,IsVector=FALSE){
    if(IsVector){
      P1=data.frame(x=P1[1],y=P1[2],z=P1[3])
      P2=data.frame(x=P2[1],y=P2[2],z=P2[3])
      P3=data.frame(x=P3[1],y=P3[2],z=P3[3])
    }
    V12=(P2-P1)/norm_vec(as.matrix(P2-P1))
    V23=(P3-P2)/norm_vec(as.matrix(P3-P2))
    VM12=(P1+P2)/2
    VM23=(P3+P2)/2
    d12=-(V12$x*VM12$x+V12$y*VM12$y+V12$z*VM12$z) # plan perpendiculaire passant milieu de AB
    d23=-(V23$x*VM23$x+V23$y*VM23$y+V23$z*VM23$z) # plan perpendiculaire passant milieu de Bc
    temp=cross(as.matrix(V12),as.matrix(V23))
    V123=temp/norm_vec(as.matrix(temp))
    d123=-(V123[1]*P1$x+V123[2]*P1$y+V123[3]*P1$z) # plan  passant parABC
    amat=array(c(V12$x,V12$y,V12$z,V23$x,V23$y,V23$z,V123[1],V123[2],V123[3]), dim=c(3,3))
    amat=t(amat)
    coe=array(c(-d12,-d23,-d123),dim=c(3,1))
    iamat=solve(amat)
    gg=iamat%*%coe
    gc=data.frame(x=gg[1],y=gg[2],z=gg[3])
    ry=norm_vec(as.matrix(gc-P1))
    return (list(center=gc, Rayon=ry))
  }

  norm_vec <- function(x) sqrt(sum(x^2))
  #DC7, DC9 dist droite passant par Fe normale heme to flavin C7, C9



  PlotFreq=function(data,mask, mtext=NULL, mode=NULL, DoPlot=TRUE, UseFit=TRUE, AltMode=FALSE)
  {
  
    if(ReadRefPDB) mask=c(1:10,mask+10)
    data=data[mask]
    if(!is.null(mode)) PrintPlus=TRUE
    else PrintPlus=FALSE
    EM11=0
    ES11=0
    EN11=0
    m11=0
    m21=0
    m22=0
    m31=0
    m32=0
    m33=0
    ss11=0
    ss21=0
    ss22=0
    ss31=0
    ss32=0
    ss33=0
    ft11=c(0,0)
    ft21=c(0,0)
    ft22=c(0,0)
    ft31=c(0,0)
    ft32=c(0,0)
    ft33=c(0,0)
    nn11=0
    nn21=0
    nn22=0
    nn31=0
    nn32=0
    nn33=0
    EM11=0
    ES=0
    EN11=0
    SS11=0
    M11=0
    SD1=0
    SD2=0
    SD3=0
    SD4=0
    SD5=0
    cat("\n1427")
    nb=trunc(max(15,length(data)/30))
    fr=1:length(data)
    colo=seq(1,nb)
    res=hist(data, breaks=nb, plot=FALSE)

    fr=1:length(data)
    xxval=paste("Value ranges for frames ",min(fr),"-",max(fr))
    xxval=""
   
 #filter
    exc=which(res$counts<max(res$counts)/20)
    cl=c(-1000,res$mids,1000)
    for(i in 1:length(exc))
    {
      #data[data>]
    }
    mm=res$mids #res$mids[res$counts>=max(res$counts)/20]

 #end filter

    avg=array(0,dim=c(5,7))
    cdz=sum(res$count!=0)
    SD5=mean(data[1:10])
    if(cdz<=3) # Less than 4 different values
    {
      set=as.numeric(levels(as.factor(as.character(data))))
      vc=c(0,0,0)
      for(j in 1:length(set))
      {
        vc[j]=sum(data==set[j])
      }
      ss=sort(vc, decreasing=TRUE, index=TRUE)

      set=set[ss$ix]
      EM11=set[1]
      EN11=ss$x[1]
      m21=set[1]
      m31=set[1]
      ft21[1]=set[1]
      ft21[2]=0
      ft31[1]=set[1]
      ft31[2]=0
      nn21=ss$x[1]
      nn31=ss$x[1]
      if(cdz>=2)
      {
        m22=set[2]
        ft22[1]=set[2]
        ft22[2]=0
        m32=set[2]
        ft32[1]=set[2]
        ft32[2]=0
        nn22=ss$x[2]
        nn32=ss$x[2]
        if(cdz==3)
        {
          m33=set[3]
          ft33[1]=set[3]
          ft33[2]=0
          nn33=ss$x[3]
        }
      }
    }
# end less that 4 different values
 
    if(cdz>3){
       try({
         # Gauss 0
         hh=hist(data,freq=TRUE,breaks=nb,,plot=FALSE,warn.unused=FALSE)
         EM11=mean(res$mids[res$counts==max(res$counts)]) #not the mean of all but the mean of egual maximum
         ES11=sd(data)
         EN11=length(data)
         ft=unlist(fitdistr(data,"normal"))
         M11=ft[1]  # center of gauss curve fit
         ss11=ft[2] # SD of gauss curve fit
         if(ft[2]!=0) SD2=abs(mean(res$mids[res$counts==max(res$counts)])-ft[1])/ft[2] # ?cart relatif du max exp et du fit center/SD fit
         vnorm=dnorm(res$mids,mean=ft[1],sd=ft[2])
         vnorm=vnorm*length(data)/sum(vnorm)
         mm=array(hh$count, dim=c(length(hh$count),1))
         mm=cbind(mm,vnorm)
         rr=apply(mm,1,min)
         SD1=sum(rr)/sum(hh$counts)
#end Gauss0

         if(FALSE) #  TRUE:guess the number of components
         {
           nq=trunc(length(data)/30)
           nq=max(3,trunc(nq/2)+1)
           NQ=min(9,nq)
           y = makemultdata(data, cuts = quantile(data, (1:nq)/(nq+1)))$y
           a = boot.comp(y = y, max.comp = 3,B=200, arbmean=TRUE, mix.type ="normalmix", epsilon = 1e-3,hist=FALSE)
           ng=length(a$p.values)
         }else ng=3
#end guess components
#triple Gauss
          hh=hist(data,freq=TRUE,breaks=nb,xlab=xxval,ylab="Counts", main=mtext,cex.main=1,plot=DoPlot,warn.unused=DoPlot)
          out=normalmixEM(data,arbvar=FALSE,arbmean=TRUE,k=ng,epsilon = 1e-03)
          ng=length(out$lambda)
#double gauss
          if(max(out$mu)-min(out$mu)<ss11)
          {
            ng=2
            out=normalmixEM(data,arbvar=FALSE,arbmean=TRUE,k=ng,epsilon = 1e-03)
          }
#single gau
          if(max(out$mu)-min(out$mu)<ss11)
          {
            ng=1
            ng1=1
            sres=array(0,dim=3)
            sres[1]=sum(hh$count)
            sres[2]=M11
            sres[3]=ss11
            y=dnorm(hh$mids,mean=M11,sd=ss11)
            yy=sum(hh$count)/sum(y)*y
            lines(hh$mids,yy, col="green")
            valT=paste("Gauss: m:",format(M11,2),"sd:",format(ss11,2),"\n", sep=" ")
            valT=paste(valT,"Exp: m:",format(mean(data),2),"sd:",format(sd(data),2),"\n", sep=" ")
            mtext(valT, side=4, cex=0.5,line=0)
          }


          if(ng>1)
          {
            rmat=array(0,dim=c(ng,3))
            rmat[,1]=out$lambda
            rmat[,2]=out$mu
            rmat[,3]=out$sigma
            st=sort(rmat[,2], index=TRUE)
            srmat=rmat[st$ix,]
            tmat=srmat[1,]

            sres=NULL
            for (j in 2:ng)
            {
              if(abs(srmat[j,2]-tmat[2])<abs(srmat[j,3]+tmat[3])/4)
              {
                mm=(srmat[j,2]*srmat[j,1]+tmat[1]*tmat[2])/(srmat[j,1]+tmat[1])
                sm=tmat[3]
                nm=srmat[j,1]+tmat[1]
                tmat=array(c(nm,mm,sm), dim=c(1,3))
              }else{
                if(is.null(sres)) sres=tmat else sres=rbind(sres,tmat)
                tmat=srmat[j,]
              }
            }
            if(is.null(sres)) sres=tmat else sres=rbind(sres,tmat)
            sres=rbind(sres,array(0,dim=c(1,3)))

            s=sum(hh$counts)
            if(length(dim(sres))>1) sres[,1]=s*sres[,1] else sres[1]=s*sres[1]
            ng1=dim(sres)[1]-1
            fit=array(0,dim=length(hh$mids))
            sdd=array(0,dim=ng1)
            valT=""
            for(i in 1:ng1)
            {
              y=dnorm(hh$mids,mean=sres[i,2],sd=sres[i,3])
              yy=sres[i,1]/sum(y)*y
              if(DoPlot) lines(hh$mids,yy,col="blue")
              fit=fit+yy
              mm=array(hh$count, dim=c(length(hh$count),1))
              mm=cbind(mm,vnorm)
              rr=apply(mm,1,min)
              sdd[i]=sum(rr)/sum(hh$counts)
              ssn=100*sres[i,1]/length(data)
              valT=paste(valT,
              format(ssn,0),"% m:"
              ,format(sres[i,2],2),
              "sd:",format(sres[i,3],2),"\n")
            }
            if(DoPlot)
            {
              lines(hh$mids,fit, col="red")
              options( "digits"=2, "scipen"=0)
              mtext(valT, side=4, cex=0.5,line=0)
            }

        }

          if(length(dim(sres))>1)
          {
            sse=sort(sres[,1], decreasing =TRUE, index=TRUE)
            tp=sres[sse$ix,]
            sres=tp
            m11=sres[1,2]
            EN11=sum(sres[,1])
          }else{
            m11=sres[2]
            nn21=EN11
          }

          if(ng1>=2)
          {
           ft21[1]=sres[1,2]
           ft22[1]=sres[2,2]
           ft21[2]=sres[1,3]
            ft22[2]=sres[2,3]
            ss21=sres[1,3]
            ss22=sres[2,3]
            nn21=sres[1,1]
            nn22=sres[2,1]
            SD3=sdd[1]
            SD4=sdd[2]
          }
          if(ng1>=3)
          {
            if(length(dim(sres))>1) ft31[1]=sres[1,2] else ft31[1]=sres[2]
            ft32[1]=sres[2,2]
            ft33[1]=sres[3,2]
            ft31[2]=sres[1,3]
            ft32[2]=sres[2,3]
            ft33[2]=sres[3,3]
            ss31=sres[1,3]
            ss32=sres[2,3]
            ss33=sres[3,3]
            nn31=sres[1,1]
            nn32=sres[2,1]
            nn33=sres[3,1]
          }

         }
      )
    }
    cat("\n1651")

  if(cdz>3 && FALSE){

    hh=hist(data,freq=TRUE,breaks=nb,xlim=c(min(mm),max(mm)),xlab=xxval,ylab="Counts", main=mtext,cex.main=1,plot=DoPlot,warn.unused=DoPlot)
    EM11=mean(res$mids[res$counts==max(res$counts)]) #not the mean of all but the mean of egual maximum
    ES11=sd(data)
    EN11=length(data)

    if(UseFit)
    {
      ft=unlist(fitdistr(data,"normal"))
      vnorm=dnorm(res$mids,mean=ft[1],sd=ft[2])
      vnorm=vnorm*length(data)/sum(vnorm)
      lines(hh$mids,vnorm, col="black")

      M11=ft[1] # center of gauss curve fit
      sS11=ft[2] # SD of gauss curve fit
      if(ft[2]!=0) SD2=abs(mean(res$mids[res$counts==max(res$counts)])-ft[1])/ft[2]  # ?cart relatif du max exp et du fit center/SD fit
      if(length(hh$mids)>25)
      {
        filt=loess(hh$count~hh$mids)
        pred=predict(filt,res$mids)
        lines(res$mids,pred, col="black") # fit liss? count des classes distrib exp
      }
      mm=array(hh$count, dim=c(length(hh$count),1))
      mm=cbind(mm,vnorm)
      rr=apply(mm,1,min)
      SD1=sum(rr)/sum(hh$count)
    }else{
        M= EM11
        SS11= ES11
    }

  # start clustering
      dist_mat = dist(data, method = 'euclidean')
      arbre = hclust(dist_mat, method ="average")
      classe2 = cutree(arbre, k = 2)
      classe3 = cutree(arbre, k = 3)
      datam=data
      s21=datam[classe2==1]
      s22=datam[classe2==2]
      if(length(s22)>length(s21))
      {
        tp=s21
        s21=s22
        s22=tp
      }
      s31=datam[classe3==1]
      s32=datam[classe3==2]
      s33=datam[classe3==3]

      tt=c(length(s31),length(s32),length(s33))# sort clustering
      re=sort(tt,index=TRUE, decreasing=TRUE)
      si=list(s31,s32,s33)
      so=si[re$ix]
      s31=so[[1]]
      s32=so[[2]]
      s33=so[[3]]
      nn21=length(s21)
      nn22=length(s22)
      nn31=length(s31)
      nn32=length(s32)
      nn33=length(s33)

  #start decomposition

      if(UseFit)
      {
        try(
          {
          ms=mean(s21)
          ec=sd(s21)
          sq21=s21[abs(s21-ms)<2*ec]
          ft21=unlist(fitdistr(sq21,"normal"))
          h21=hist(s21,freq=TRUE,breaks=nb,xlim=c(min(s21),max(s21)),xlab=xxval,ylab="Counts", main=paste("2C1b_",mtext),cex.main=1,plot=DoPlot,warn.unused=DoPlot)
          vnorm=dnorm(h21$mids,mean=ft21[1],sd=ft21[2])
          vnorm=vnorm*length(sq21)/sum(vnorm)
          mm=array(hh21$count, dim=c(length(hh21$count),1))
          mm=cbind(mm,vnorm)
          rr=apply(mm,1,min)
          SD3=sum(abs(rr))/sum(h21$counts)
          try(if((PrintPlus && DoPlot) || (PrintPlus && !DoPlot && (SD3<SD1)))lines(h21$mids,vnorm, col="blue"), silent=TRUE)
          m21=mean(s21[res$counts==max(res$counts)])
          ss21=sd(s21)
          }
          ,silent=TRUE
        )
      }


      if(UseFit)
      {
        try(
          {
            if(length(s22)/EN11>=0.05 )
            {
              ms=mean(s22)
              ec=sd(s22)
              sq22=s22[abs(s22-ms)<2*ec]
              ft22=unlist(fitdistr(sq22,"normal"))
              h22=hist(s22,freq=TRUE,breaks=nb,xlim=c(min(s22),max(s22)),xlab=xxval,ylab="Counts", main=paste("2C2_",mtext),cex.main=1,plot=FALSE,warn.unused=FALSE)
              vnorm=dnorm(h22$mids,mean=ft22[1],sd=ft22[2])
              vnorm=vnorm*length(sq22)/sum(vnorm)
              mm=array(hh22$count, dim=c(length(hh22$count),1))
              mm=cbind(mm,vnorm)
              rr=apply(mm,1,min)
              SD4=sum(rr)/sum(h22$counts)
              if(SD4<0.3)
              {
                h22=hist(s22,freq=TRUE,breaks=nb,xlim=c(min(s22),max(s22)),xlab=xxval,ylab="Counts", main=paste("2C2_",mtext),cex.main=1,plot=DoPlot,warn.unused=DoPlot)
                if((PrintPlus && DoPlot) || (PrintPlus && !DoPlot && (norm2<norm1))) lines(h22$mids,vnorm, col="blue")
                m22=mean(s22[res$counts==max(res$counts)])
                ss22=sd(s22)
              }
            }
          }
          ,silent=TRUE
        )
      }


      if(UseFit)
      {
        try(
          {
            ms=mean(s31)
            ec=sd(s31)
            sq31=s31[abs(31-ms)<2*ec]
            ft31=unlist(fitdistr(sq31,"normal"))
            h31=hist(s31,freq=TRUE,breaks=nb,xlim=c(min(s31),max(s31)),xlab=xxval,ylab="Counts", main=paste("3C1_",mtext),cex.main=1,plot=DoPlot,warn.unused=DoPlot)
            vnorm=dnorm(h31$mids,mean=ft31[1],sd=ft31[2])
            vnorm=vnorm*length(sq31)/sum(vnorm)
            norm3=sum(abs(vnorm-h31$counts))/sum(h31$counts)
            if((PrintPlus && DoPlot) || (PrintPlus && !DoPlot && (norm3<norm2) && (norm2<norm1))) lines(h31$mids,vnorm, col="green")
            m31=mean(s31[res$counts==max(res$counts)])
            ss31=sd(s31)
          }
          ,silent=TRUE
        )
      }

      if(UseFit)
      {
        try(
          {
            if(length(s32)/EN11>=0.05)
            {
              ms=mean(s32)
              ec=sd(s32)
              sq32=s32[abs(s32-ms)<2*ec]
              ft32=unlist(fitdistr(sq32,"normal"))
              h32=hist(s32,freq=TRUE,breaks=nb,xlim=c(min(s32),max(s32)),xlab=xxval,ylab="Counts", main=paste("3C2_",mtext),cex.main=1,plot=FALSE,warn.unused=FALSE)
              vnorm=dnorm(h32$mids,mean=ft32[1],sd=ft32[2])
              vnorm=vnorm*length(sq32)/sum(vnorm)
              err=sum(abs(vnorm-h32$counts))/sum(h32$counts)
              if(err<0.3)
              {
              h32=hist(s32,freq=TRUE,breaks=nb,xlim=c(min(s32),max(s32)),xlab=xxval,ylab="Counts", main=paste("3C2_",mtext),cex.main=1,plot=DoPlot,warn.unused=DoPlot)
              if((PrintPlus && DoPlot) || (PrintPlus && !DoPlot && (norm3<norm2) && (norm2<norm1)))  lines(h32$mids,vnorm, col="green")
              m32=mean(s32[res$counts==max(res$counts)])
              ss32=sd(s32)
              }
            }
          }
          ,silent=TRUE
        )
      }

        if(UseFit)
        {
          try(
            {
              if(length(s33)/EN11>=0.05 )
              {
                ms=mean(s33)
                ec=sd(s33)
                sq33=s23[abs(s33-ms)<2*ec]
                ft33=unlist(fitdistr(sq33,"normal"))
                h33=hist(s33,freq=TRUE,breaks=nb,xlim=c(min(s33),max(s33)),xlab=xxval,ylab="Counts", main=paste("3C3_",mtext),cex.main=1,plot=FALSE,warn.unused=FALSE)
                vnorm=dnorm(h33$mids,mean=ft33[1],sd=ft33[2])
                vnorm=vnorm*length(sq33)/sum(vnorm)
                err=sum(abs(vnorm-h33$counts))/sum(h33$counts)
                if(err<0.3)
                {
                h33=hist(s33,freq=TRUE,breaks=nb,xlim=c(min(s33),max(s33)),xlab=xxval,ylab="Counts", main=paste("3C3_",mtext),cex.main=1,plot=DoPlot,warn.unused=DoPlot)
                if((PrintPlus && DoPlot) || (PrintPlus && !DoPlot && (norm3<norm2) && (norm2<norm1)))  lines(h33$mids,vnorm, col="green")
                m33=mean(s33[res$counts==max(res$counts)])
                ss33=sd(s33)
                }
              }
            }
            ,silent=TRUE
          )
        }
  }
  if(AltMode)
  {
    data1=(data-1)%%9+1
    ft21[1]=sum(data1==1)
    ft21[2]=sum(data1==2)
    ft22[1]=sum(data1==3)
    ft22[2]=sum(data1==4)
    ft31[1]=sum(data1==5)
    ft31[2]=sum(data1==6)
    ft32[1]=sum(data1==7)
    ft32[2]=sum(data1==8)
    ft33[1]=sum(data1==9)
    ft33[2]=sum(data1==10)
  }
    avg[1,1]=EM11 # Max Exp
    avg[2,1]=ES11 # SD Exp
    avg[3,1]=EN11 # Nb Exp
    avg[3,2]=SD1
    avg[2,2]=SD2
    avg[4,1]=M11  # mean gauss 1/1
    avg[1,2]=ss11 # sd gauss 1/1
    avg[4,2]=SD3
    avg[5,2]=SD4
    avg[5,1]=SD5
    avg[1,3]=ft21[1] # mean class 1/2 Fit Gauss
    avg[2,3]=ft22[1] # mean class 2/2 Fit Gauss
    avg[3,3]=ft31[1] # mean class 1/3 Fit Gauss
    avg[4,3]=ft32[1] # mean class 2/3 Fit Gauss
    avg[5,3]=ft33[1] # mean class 3/3 Fit Gauss
    avg[1,4]=ft21[2] # SD class 1/2 Fit Gauss
    avg[2,4]=ft22[2] # SD class 2/2 Fit Gauss
    avg[3,4]=ft31[2] # SD class 1/3 Fit Gauss
    avg[4,4]=ft32[2] # SD class 2/3 Fit Gauss
    avg[5,4]=ft33[2] # SD class 3/3 Fit Gauss
    avg[1,5]=m21 # Max class 1/2
    avg[2,5]=m22 # MAX class 2/2
    avg[3,5]=m31 # MAX class 1/3
    avg[4,5]=m32 # MAX class 2/3
    avg[5,5]=m33 # MAX class 3/3
    avg[1,6]=ss21 # SD class 1/2
    avg[2,6]=ss22 # SD class 2/2
    avg[3,6]=ss31 # SD class 1/3
    avg[4,6]=ss32 # SD class 2/3
    avg[5,6]=ss33 # SD class 3/3
    avg[1,7]=nn21 # count 1/2
    avg[2,7]=nn22 # count 2/2
    avg[3,7]=nn31 # count 1/3
    avg[4,7]=nn32 # count 2/3
    avg[5,7]=nn33 # count 3/3
    par(mar=c(5.1,4.1,4.1,2.1))
    cat("-1397")
    return(avg)
  }


  SumSS=function(ss, fname, AltParm=nULL)
  {
    mode="new"
    ss[is.na(ss)]=0
    mres=data.frame()
    rrnames=c()
    rr=array(0,dim=c(6,4))
    for(i in 1:(dim(ss)[1]-1))
    {
      mm=c(0,0,0)
      if(mode=="old")
      {
      nn=ss[i,"CNT12"]+ss[i,"CNT22"]
      V1=c(ss[i,"GCen12"]!=0,ss[i,"GCen22"]!=0)
      V2=c(ss[i,"GCen13"]!=0,ss[i,"GCen23"]!=0,ss[i,"GCen33"]!=0)
      S1=abs(ss[i,"GCen12"]-ss[i,"GCen22"])>=0.5*(ss[i,"GSd12"]+ss[i,"GSd22"])
      S2=c(abs(ss[i,"GCen13"]-ss[i,"GCen23"])>=0.5*(ss[i,"GSd13"]+ss[i,"GSd23"]),abs(ss[i,"GCen33"]-ss[i,"GCen23"])>=0.5*(ss[i,"GSd33"]+ss[i,"GSd23"]))

      if((sum(V2)==3) && (sum(S2)==2) ) #case 1: triple Gauss decomposition  valid
      {
        rr[3:5,1]=c(ss[i,"GCen13"],ss[i,"GCen23"],ss[i,"GCen33"])
        rr[3:5,2]=c(ss[i,"GSd13"],ss[i,"GSd23"],ss[i,"GSd33"])
        rr[3:5,3]=c(ss[i,"CNT13"],ss[i,"CNT23"],ss[i,"CNT33"])
        rr[6,1:3]=c(1,0,0)
        mm[3]=3
      }
      if(mm[3]==0 && sum(V2)==2 && S2[1]) #case 2: only 2 of the triple Gauss decomposition  valid
      {
        rr[3:4,1]=c(ss[i,"GCen13"],ss[i,"GCen23"])
        rr[3:4,2]=c(ss[i,"GSd13"],ss[i,"GSd23"])
        rr[3:4,3]=c(ss[i,"CNT13"],ss[i,"CNT23"])
        rr[5,1:3]=c(0,0,0)
        rr[6,1:3]=c(2,0,0)
        mm[3]=2
      }

      if(mm[3]==0 && sum(V1)==2 && S1[1]) #case 3:  the double Gauss decomposition  valid but not the triple

      {
        rr[3:4,1]=c(ss[i,"GCen12"],ss[i,"GCen22"])
        rr[3:4,2]=c(ss[i,"GSd12"],ss[i,"GSd22"])
        rr[3:4,3]=c(ss[i,"CNT12"],ss[i,"CNT22"])
        rr[5,1:3]=c(0,0,0)
        rr[6,1:3]=c(3,0,0)
        mm[2]=2
      }

      if(mm[2]==0 && S2[1]) #case 4:  the double Gauss decomposition  not valid
      {
        rr[3,1]=ss[i,"GCen12"]
        rr[3,2]=ss[i,"GSd12"]
        rr[3,3]=ss[i,"CNT12"]
        rr[4,1:3]=c(0,0,0)
        rr[5,1:3]=c(0,0,0)
        rr[6,1:3]=c(4,0,0)
        mm[2]=1
      }
      if(ss[i,"GExCover"]<0.6 || ss[i,"NormTest"]>0.8)
      {
        rr[1,1:3]=c(ss[i,"MaxE"],ss[i,"SdEx"],ss[i,"CNT12"]+ss[i,"CNT22"])
        rr[3,1:3]=c(0,0,0)
        rr[4,1:3]=c(0,0,0)
        rr[5,1:3]=c(0,0,0)
        rr[6,1:3]=c(5,0,0)
        mm[2]=1
      }
      rr[1,1]=ss[i,"MaxE"]
      rr[1,2]=ss[i,"SdEx"]
      rr[1,3]=ss[i,"CNT12"]+ss[i,"CNT22"]

      rr[2,1]=ss[i,"GCen11"]
      rr[2,2]=ss[i,"GSd11"]
      rr[2,3]=ss[i,"CNT12"]+ss[i,"CNT22"]

      rr[1,4]=ss[i,"DifExFt"]
      rr[2,4]=ss[i,"NormTest"]
      rr[3,4]=ss[i,"DifQtl"]
      }else{
      if(ss[i,"GCen33"]!=0)
      {
        rr[3:5,1]=c(ss[i,"GCen13"],ss[i,"GCen23"],ss[i,"GCen33"])
        rr[3:5,2]=c(ss[i,"GSd13"],ss[i,"GSd23"],ss[i,"GSd33"])
        rr[3:5,3]=c(ss[i,"CNT13"],ss[i,"CNT23"],ss[i,"CNT33"])
        rr[6,1:3]=c(1,0,0)
        mm[3]=3
      }
      if(ss[i,"GCen22"]!=0 && ss[i,"GCen33"]==0){
        rr[3:4,1]=c(ss[i,"GCen12"],ss[i,"GCen22"])
        rr[3:4,2]=c(ss[i,"GSd12"],ss[i,"GSd22"])
        rr[3:4,3]=c(ss[i,"CNT12"],ss[i,"CNT22"])
        rr[5,1:3]=c(0,0,0)
        rr[6,1:3]=c(3,0,0)
        mm[2]=2
      }
      if(ss[i,"GCen22"]==0 ){
        rr[3,1]=ss[i,"GCen12"]
        rr[3,2]=ss[i,"GSd12"]
        rr[3,3]=ss[i,"CNT12"]
        rr[4,1:3]=c(0,0,0)
        rr[5,1:3]=c(0,0,0)
        rr[6,1:3]=c(4,0,0)
        mm[2]=1
      }
      rr[1,1]=ss[i,"MaxE"]
      rr[1,2]=ss[i,"SdEx"]
      rr[1,3]=ss[i,"CNT12"]+ss[i,"CNT22"]
      rr[1,4]=ss[i,"DifExFt"]

      rr[2,1]=ss[i,"GCen11"]
      rr[2,2]=ss[i,"GSd11"]
      rr[2,3]=ss[i,"CNT12"]+ss[i,"CNT22"]
      rr[2,4]=ss[i,"NormTest"]

      rr[3,4]=ss[i,"DifQtl"]
      rr[4,4]=ss[i,"StartVal" ]
      if(!is.null(AltParm) && (i%in%AltParm) )
      {
        rr[3:5,1]=c(ss[i,"GCen12"],ss[i,"GSd12"],ss[i,"GCen22"])
        rr[3:5,2]=c(ss[i,"GSd22"],ss[i,"GCen13"],ss[i,"GSd13"])
        rr[3:5,3]=c(ss[i,"GCen23"],ss[i,"GSd23"],ss[i,"GCen33"])
        rr[6,1:3]=c(ss[i,"GSd33"],0,0)
      }
    }

      rrf=data.frame(array(rr, dim=c(6,4)))
      nms=rownames(ss)[i]
      names=c(paste(nms,"Exp", sep="_"),paste(nms,"Gauss0", sep="_"),paste(nms,"Gauss1", sep="_"),paste(nms,"Gauss2", sep="_")
                                              ,paste(nms,"Gauss3", sep="_"),paste(nms,"Desc", sep="_"))
      if(i==1)
      {
        mres=rrf
        ccnames=names
      }
      else
      {
          mres=rbind(mres,rrf)
          ccnames=c(ccnames,names)
      }
    }
    mres=rbind(mres,c(ss[dim(ss)[1],1],0,0,0))
    cnames=c(paste(fname,"Center",sep="_"),paste(fname,"StdDev",sep="_"),paste(fname,"Count",sep="_"),paste(fname,"Stat",sep="_"))
    rnames=c("Exp","Gauss0","Gauss1","Gauss2","Gauss3","Desc")
    colnames(mres)=cnames
    rownames(mres)=c(ccnames,"CTC")
    return (mres)
  }


  FuzeResultSet2=function(ASet,Desc,expname, ResFile=NULL, select=2,AltProt=NULL)
  {
    nrec=dim(ASet)[2]/4
    nparm=(dim(ASet)[1]-1)/6

    Mat=rep(seq(1,3*length(select)+2),nrec)
    for(i in 1:nrec)
    {
      m1=array(as.numeric(Desc[i,c(1:2,4)]),dim=c(3,1))
      if(i==1)
      {
        mt=m1
        for(j in 1:(3*length(select)+1)) mt=cbind(mt,m1)
      }
      else
        {
          for(j in 1:(3*length(select)+2)) mt=cbind(mt,m1)
        }
    }
    Mat=rbind(Mat,mt)
    for(i in 1:nrec)
    {
      pc=1+(i-1)*4
      for(j in 1:nparm)
      {
        pr0=1+(j-1)*6
        dbloc=ASet[pr0:(pr0+5),pc:(pc+3)]
        InitVal=dbloc[4,4]
        if(j==23)
        {
          u=1
          u=2
        }
        if(is.null(AltProt) || !(j%in%AltProt))
        {
        if(dbloc[2,1]==0)dbloc[2,]=dbloc[1,]
        if(dbloc[3,1]==0)dbloc[3,]=dbloc[2,]
        if(dbloc[3,1]==0) dbloc[3,]=dbloc[1,]
        }
        t1=array(c(dbloc[select,1],dbloc[select,2],dbloc[select,3],dbloc[1,1],InitVal), dim=c(1,3*length(select)+2))
        if (j==1) t2=t1 else t2=rbind(t2,t1)
      }
      if(i==1) t3=t2
      else t3=cbind(t3,t2)
    }
    rr=rbind(Mat,t3)

    for(j in 1:length(select))
    {
      if(j==1)
      {
        m1=paste("Moy",j,sep="")
        m2=paste("Std",j,sep="")
        m3=paste("Cnt",j,sep="")
      }else{
      m1=c(m1,paste("Moy",j,sep=""))
      m2=c(m2,paste("Std",j,sep=""))
      m3=c(m3,paste("Cnt",j,sep=""))
      
      }
    }
    cn=c(m1,m2,m3,"MoyExp","InitVal")
    cn=rep(cn,nrec)
    colnames(rr)=cn
    fcc=c()
    st=seq(1,length(rownames(ASet)),6)
    for(i in st)
    {
      tt=as.array(unlist(str_split(rownames(ASet)[i],pattern="_")))
      fcc=c(fcc,tt[1])
    }
    uu=array(unlist(str_split(rownames(ASet),"_")),dim=c(2,dim(ASet)[1]))[1,]
    f=seq(1,dim(ASet)[1]-1,by=6)
    uuu=uu[f]
    bn=c("Key","ID","OMM1","CFTyp1",uuu)
    rownames(rr)=bn
    if(!is.null(ResFile)) write.table(rr, file=ResFile)
    return (rr)
  }
  
  FuzeResultSet3=function(ASet,Desc,expname, ResFile=NULL, select=2)
  {
    nrec=dim(ASet)[2]/4
    nparm=(dim(ASet)[1]-1)/6
    
    Mat=rep(seq(1,3*length(select)+2),nrec)
    for(i in 1:nrec)
    {
      m1=array(as.numeric(Desc[i,c(1:2,4)]),dim=c(3,1))
      if(i==1)
      {
        mt=m1
        for(i in 1:(3*length(select)+1)) mt=cbind(mt,m1)
      }
      else
      {
        for(i in 1:(3*length(select)+2)) mt=cbind(mt,m1)
      }
    }
    Mat=rbind(Mat,mt)
    
    for(i in 1:nrec)
    {
      pc=1+(i-1)*4
      for(j in 1:nparm)
      {
        pr0=1+(j-1)*6
        dbloc=ASet[pr0:(pr0+5),pc:(pc+3)]
        InitVal=dbloc[4,4]
        if(dbloc[2,1]==0)dbloc[2,]=dbloc[1,]
        if(dbloc[3,1]==0)dbloc[3,]=dbloc[2,]
        if(dbloc[3,1]==0) dbloc[3,]=dbloc[1,]
        t1=array(c(dbloc[select,1],dbloc[select,2],dbloc[select,3],dbloc[1,1],InitVal), dim=c(1,3*length(select)+2))
        if (j==1) t2=t1 else t2=rbind(t2,t1)
      }
      if(i==1) t3=t2 else t3=cbind(t3,t2)
    }
    rr=rbind(Mat,t3)
    
    for(j in 1:length(select))
    {
      if(j==1)
      {
        m1=paste("Moy",j,sep="")
        m2=paste("Std",j,sep="")
        m3=paste("Cnt",j,sep="")
      }else{
        m1=c(m1,paste("Moy",j,sep=""))
        m2=c(m2,paste("Std",j,sep=""))
        m3=c(m3,paste("Cnt",j,sep=""))
      }
    }
    cn=c(m1,m2,m3,"MoyExp","InitVal")
    cn=rep(cn,nrec)
    colnames(rr)=cn
    fcc=c()
    st=seq(1,length(rownames(ASet)),6)
    for(i in st)
    {
      tt=as.array(unlist(str_split(rownames(ASet)[i],pattern="_")))
      fcc=c(fcc,tt[1])
    }
    uu=array(unlist(str_split(rownames(ASet),"_")),dim=c(2,dim(ASet)[1]))[1,]
    f=seq(1,dim(ASet)[1]-1,by=6)
    uuu=uu[f]
    bn=c("Key","ID","OMM1","CFTyp1",uuu)
    rownames(rr)=bn
    if(!is.null(ResFile)) write.table(rr, file=ResFile)
    return (rr)
  }


 # tp=sapply(rownames(ASet),str_split,pattern="_")
 # sapply(nn,unlist)

  selectRec=function(PdfSet, WorkSet, NbCFF=1, SelCFF=1, InitPos=NULL, InitOxDist=NULL, MeanType, ErrorType, CountType,InitType,OrientType=NULL)
  {
    #type=11
    #MeanType=1 # 1: use mode, 3: use center of gauss fit
    tab=read.table(PdfSet)
    tab=as.matrix(tab)
    Ws=read.table(WorkSet, sep='\t')
    sub=as.matrix(Ws[,c(2,4)])
    filt2=which(tab[4,]==2)  #  CFF 2present
    filt1=which(tab[4,]==1)  # CFF 1 present
    idd2=as.numeric(levels(factor(tab[3,filt2]))) #OMM with CFf 2
    idd1=as.numeric(levels(factor(tab[3,filt1]))) #OMM with CFf 1
    if(length(idd2)>=length(idd1))
      {
      idd=idd2[idd2%in%idd1]
      }else{
      idd=idd1[idd1%in%idd2]
      }
    tp=array(FALSE,dim=dim(tab)[2])
    if(NbCFF==2) tp=array(FALSE,dim=dim(tab)[2]) 
    if(NbCFF==1) tp=array(TRUE,dim=dim(tab)[2])
    if(length(idd)>0)
    {
      for(k in 1:length(idd))
      {
        tpi=which(tab[3,]==idd[k])
        if(NbCFF==2) tp[tpi]=TRUE
        if(NbCFF==1) tp[tpi]=FALSE
      }
    }
    if(!is.null(InitPos))
    {
      filt=which((tab[15,]== InitPos)  & (tab[1,]==3) & (tab[4,]==SelCFF))
      idd=as.numeric(levels(factor(tab[3,filt])))
      tp1=array(FALSE,dim=dim(tab)[2])
      for(i in 1:length(idd))
      {
        tpi=which(tab[3,]==idd[i])
        tp1[tpi]=TRUE
      }
    }
    if(!is.null(InitOxDist))
    {
      filt=which((tab[15,]<= InitOxDist) & (tab[1,]==2)& (tab[4,]==SelCFF))
      idd=as.numeric(levels(factor(tab[3,filt])))
      tp2=array(FALSE,dim=dim(tab)[2])
      for(i in 1:length(idd))
      {
        tpi=which((tab[3,]==idd[i]))
        tp2[tpi]=TRUE
      }
    }
    sel=tp & (tab[1,] %in% MeanType) & (tab[4,]==SelCFF) #   Select mode or mean of fit
    sel1=tp & (tab[1,] %in% ErrorType) & (tab[4,]==SelCFF) #   select error
    sel2=tp & (tab[1,] %in% CountType) & (tab[4,]==SelCFF) #   select Count
    sel3=tp & (tab[1,] %in% InitType) & (tab[4,]==SelCFF) #   select InitVal
    if(!is.null(OrientType)) sel4=tp & (tab[1,] %in% OrientType) & (tab[4,]==SelCFF) #   select InitVal
    moresel=rep(TRUE,dim(tab)[2])
    if(!is.null(InitPos)) moresel=(moresel & tp1)
    if(!is.null(InitOxDist)) moresel=(moresel & tp2)
    tt=tab[3,which(moresel)]
    setsel=array(FALSE,dim=dim(tab)[2])
    for(i in 1:length(tt))
    {
      tp=which(tab[3,]==tt[i])
      setsel[tp]=TRUE
    }
    sel=sel & setsel
    sel1=sel1 & setsel
    sel2=sel2 & setsel
    if(!is.null(OrientType)) sel4=sel4 & setsel else rep(FALSE,dim(tab)[2])
    return(data.frame(sel=sel,sel1=sel1,sel2=sel2,sel3=sel3,sel4=sel4))
  }


  selectRec1=function(PdfSet, WorkSet, OMMId, SelCFF=1)
  {
    MeanType=1 # 1: use mode, 3: use center of gauss fit
    tab=read.table(PdfSet)
    tab=as.matrix(tab)
    Ws=read.table(WorkSet, sep='\t')
    Id=Ws[Ws[,3]%in% OMMId,2]
    id1=tab[3,]%in%as.numeric(Id)
    sel=id1 & (tab[1,]==MeanType) & (tab[4,]==SelCFF) #   Select mode or mean of fit
    sel1=id1 & (tab[1,]==2) & (tab[4,]==SelCFF) #   select error
    return(data.frame(sel=sel,sel1=sel1))
  }


  SelectOMMSet=function()
  {
  OMMSel2=TRUE
  OMMSel2=FALSE

  if (OMMSel1)
  {
    OMMA=c("OMM02","OMM14","OMM15","OMM18","OMM19","OMM39")
    OMMB=c("OMM03","OMM13","OMM23","OMM25","OMM24","OMM26B","OMM30","OMM47","OMM56")
    OMMC=c("OMM05","OMM09","OMM11","OMM46")
    OMMD=c("OMM06","OMM08","OMM37","OMM37B","OMM37C")
    OMMF=c("OMM22","OMM27","OMM21")
    OMMG=c("OMM28","OMM26","OMM32","OMM32B","OMM53B")
    OMMH=c("OMM31","OMM33","OMM34","OMM54","OMM55")
    OMMJ=c("OMM36","OMM38","OMM43","OMM44")
    OMMK=c("OMM45","OMM52","OMM57")
    OMM2AE_effA=c("OMM1B","OMM41","OMM07","OMM01","OMM7B")
    OMM2AE_effB=c("OMM49","OMM50","OMM51","OMM04","OMM42","OMM48")
    OMM2AP_prodA=c("OMM50","OMM42","OMM49")
    OMM2AP_prodB=c("OMM01","OMM1B")
    OMM2AP_prodC=c("OMM07","OMM7B","OMM41","OMM04","OMM4B")
    OMMId=OMM2AP_prodC
    asel=selectRec1(PdfSet, WorkSet, OMMId, SelCFF)
  }
  if (OMMSel2)
  {
    OMMA=c("OMM02","OMM14","OMM15","OMM18","OMM19","OMM39")
    OMMB=c("OMM03","OMM13","OMM23","OMM25","OMM24","OMM26","OMM26B","OMM28","OMM30","OMM32"
           ,"OMM32B","OMM47","OMM53","OMM56")
    OMMC=c("OMM05","OMM09","OMM11","OMM46")
    OMMD1=c("OMM06","OMM08","OMM37","OMM37B","OMM37C")
    OMMD2=c("OMM43","OMM44")
    OMME=c("OMM22","OMM27","OMM21")
    OMMF=c("OMM31","OMM33","OMM34","OMM54","OMM55")
    OMMG=c("OMM35")
    OMMJ=c("OMM36","OMM38","OMM45","OMM52","OMM57")
  }
  OMMId=OMMB
  asel=selectRec1(PdfSet, WorkSet, OMMId, SelCFF)
  return (sel)
  }
  #endsel


  PrintRec=function(type, PdfSet, WorkSet, sel, sel1, Sort=TRUE, Bar=TRUE, Legend=NULL,Color=NULL)
  {
    MeanType=1 # 1: use mode, 3: use center of gauss fit
    tab=read.table(PdfSet)
    tab=as.matrix(tab)
    Ws=read.table(WorkSet, sep='\t')
    rownames(tab)=c("Type","Order","RunId","CFF id","overlap","CFF-PHE Ang","CFF-Phe dist","CFF-Phe ProjDist","Ang CFF-HEM","Fe-O CFFP dist",
                    "Angl CFFP-NaNc","Angl CFFP-NbNd","CFF-Fe-O Proj Dist", "orient","CTC","COC","Shorter O to ox CFF")
    mt=tab[,sel] #mean
    mte=tab[,sel1] # deviation
    d=1:length(sel)
    tx=tab[3,which(sel)]
    rr=vector()
    for (i in 1:length(tx))
    {
      u=which(Ws[,2]==tx[i])
      rr=c(rr,Ws[u[1],3])
    }
    colnames(mt)=rr
    if(length(Sort)==1)
    {
      if(Sort)
      {
        ss=sort(mt[type,], index.return=TRUE,decreasing=TRUE)
        dat=mt[type,ss$ix]
        date=mte[type,ss$ix]
      }else{
        ss=1:dim(mt)[2]
        dat=mt[type,ss]
        date=mte[type,ss]
      }
    }else{
      dat=mt[type,Sort]
      date=mte[type,Sort]
    }
    if(is.null(Color)) Color=rep("grey",length(dat))
    if(is.null(Legend)) mm=rownames(mt)[type]
    else mm=Legend
    if(Bar) bp=barplot(dat, main=mm,las=2, col=Color)
    else bp=plot(dat,main=mm,pch=19)
    arrows(x0 = bp, y0 = dat + date,y1 = dat- date,angle = 90,code = 3,length = 0.1)
    print(c("Mean: ",mean(mt[type,])))
    print(c("Stddev.",sd(mt[type],)))
  }


  DefOrientSet=function(mat, RefParm, Thresh,Multi=1)
  {
    ss=seq(1,ncol(mat),Multi)
    sel=which(mat[RefParm,ss]<Thresh) # select cff face
    sel1=which(mat[RefParm,ss]>=Thresh) # select cff face
    return(list(Set1=sel,Set2=sel1))
  }


  SplitPrintStat=function(mat,matSD,names,ParmId, OrientSets,CorMode=0, DoSort=NULL, Identify=FALSE, Title="", YAxis="", filter="none",Multi=3, DS=FALSE,LabScale=1, Ofset=0)
  {
    nc=ncol(mat)
    cs=nc/Multi #number of OMMs
    OrientSet1=OrientSets[[1]]
    OrientSet2=OrientSets[[2]]
    OrientSet1=OrientSet1[!is.na(OrientSet1)]
    OrientSet2=OrientSet2[!is.na(OrientSet2)]
    if(CorMode==1 || CorMode==3) mat=180-mat
    if(CorMode==4) mat[mat>90]=180-mat[mat>90]
    
    name=names$id
    xlab=names$xlab
    TypParm=ParmId
    if(length(OrientSet1)>0)
    {
      
      pos=(OrientSet1-1)*Multi+1
      M=mat[TypParm,pos]
      S=matSD[TypParm,pos]
      sn=name[OrientSet1]
      if(Multi>=2)
      {
        Mb=mat[TypParm,pos+1]
        Sb=matSD[TypParm,pos+1]
        colb=c(rep(rgb(0,0.5,0),length(Mb)))
      }
      if(Multi>=3)
      {
        Mc=mat[TypParm,pos+2]
        Sc=matSD[TypParm,pos+2]
        colc=c(rep("orange",length(Mc)))
      }
      I=1:length(M)
      if(length(DoSort)==1)
      {
        ss=sort(M, index=TRUE)
        ord=ss$ix
      }else ord=DoSort[OrientSet1]
      if(!is.null(DoSort))
      {
          nn=sn[ord]
          M=M[ord]
          Mb=Mb[ord]
          Mc=Mc[ord]
          Sb=Sb[ord]
          Sc=Sc[ord]
          S=S[ord]
      }else{
          nn=sn
      }
      so=M
      si=S
      sob=Mb
      sib=Sb
      soc=Mc
      sic=Sc
      lm1 = median(M)
      lm2 = 0
      sd1=sd(M)
      sd2=0
      col=c(rep("black",length(M)))
      ccol=col
      ccolb=colb
      ccolc=colc
      symb=rep(19,length(M))
      ssymb=symb
      nnn=nn
      ORD=ord
      }

    if(length(OrientSet2)>0)
    {
      pos=(OrientSet2-1)*Multi+1
      M1=mat[TypParm,pos]
      S1=matSD[TypParm,pos]
      sn1=name[OrientSet2]
      if(Multi>=2)
      {
        M1b=mat[TypParm,pos+1]
        S1b=matSD[TypParm,pos+1]
        col1b=c(rep(rgb(0,0.5,0),length(M1b)))
      }
      if(Multi>=3)
      {
        M1c=mat[TypParm,pos+2]
        S1c=matSD[TypParm,pos+2]
        col1c=c(rep("orange",length(M1c)))
      }
      I=1:length(M1)
      cordata=M1
      cordatab=M1b
      cordatac=M1c
      if(length(DoSort)==1)
      {
        ss=sort(M1, index=TRUE)
        ord1=ss$ix
      }else ord1=DoSort[OrientSet2]
      if(!is.null(DoSort))
      {
        nn1=sn1[ord1]
        M1=M1[ord1]
        M1b=M1b[ord1]
        M1c=M1c[ord1]
        S1b=S1b[ord1]
        S1c=S1c[ord1]
        S1=S1[ord1]
      }else{
        nn1=sn1
      }
      so=M1
      si=S1
      sob=M1b
      sib=S1b
      soc=M1c
      sic=S1c
      lm1 = median(M1)
      lm2 = 0
      sd1=sd(M1)
      sd2=0
      col1=c(rep("black",length(M1)))
      ccol=col1
      ccolb=col1b
      ccolc=col1c
      symb1=rep(20,length(M1))
      ssymb=symb1
      nnn=nn1
      ORD=ord1
    }
    resultat =0
    if (length(OrientSet1)>0 && length(OrientSet2)>0)
    {
      cordata=c(M,M1)
      cordatab=c(Mb,Mb1)
      cordatac=c(Mc,Mc1)
      if(!is.null(DoSort))
      {
        nn=c(s[ss$ix],s1[ss1$ix])
      }else{
        nn=c(s,s1)
      }
      
      so=c(M,M1)
      si=c(S,S1)
      ccol=c(col,col1)
      ssymb=c(symb,symb1)
      ORD=c(ord,ord1)
      nnn=c(nn,nn1)
      if(Multi>=2)
      {
        sob=c(Mb,M1b)
        sib=c(Sb,S1b)
        ccolb=c(colb,colb1)
      }
      if(Multi>=3)
      {
        soc=c(Mc,M1c)
        sic=c(Sc,S1c)
        ccolc=c(colc,colc1)
      }
      cordata=so
      cordatab=sob
      cordatac=soc
      lm1 = median(M)
      lm2 = median(M1)
      sd1=sd(M)
      sd2=sd(M1)
      if(length(M)>1 && length(M1)>1) resultat = t.test(M, M1, var.equal=F) else resultat=NULL
    }
    ii=1:length(so)
    SLY="n"
    if(filter=="ANG")
    {
      so=so+Ofset
      sob[sob!=0]=sob[sob!=0]+Ofset
      soc[soc!=0]=soc[soc!=0]+Ofset
      so[1]=so[1]%%360
      for(i in 1:length(so))
      {
        if(i>1)
        {
        a=seq(so[i]-1080,so[i]+1080,360)
        so[i]=a[abs(a-so[i-1])==min(abs(a-so[i-1]))][1]
        }
        if(Multi>=2)
        {
        a=seq(sob[i]-1080,sob[i]+1080,360)
        sob[i]=a[abs(a-so[i])==min(abs(a-so[i]))][1]
        }
        if(Multi>=3)
        {
        a=seq(soc[i]-1080,soc[i]+1080,360)
        soc[i]=a[abs(a-so[i])==min(abs(a-so[i]))][1]
        }
      }
     
    }
      if(filter=="ORI")
      {
       SLY="y"
        lab=c("C2(C12N3)","C2N4","N4(N9)","C8(C8)","C3(C14N7)","C3O2","O2(O13)","C1(C10N1)","O1(O11)")
        so[1]=so[1]%%9
        for(i in 1:length(so))
        {
          if(i>1)
          {
          a=seq(so[i]-27,so[i]+27,9)
          so[i]=a[abs(a-so[i-1])==min(abs(a-so[i-1]))][1]
          }
          if(Multi>=2)
          {
          a=seq(sob[i]-27,sob[i]+27,9)
          sob[i]=a[abs(a-so[i])==min(abs(a-so[i]))][1]
          }
          if(Multi>=3)
          {
          a=seq(soc[i]-27,soc[i]+27,9)
          soc[i]=a[abs(a-so[i])==min(abs(a-so[i]))][1]
          }
        }
        mm=trunc(min(c(so,sob,soc)))
        if(mm<1)
        {
          so=so+18
          sob=sob+18
          soc=soc+18
        }
      }
    cordata=so
    cordatab=sob
    cordatac=soc
    
    if(length(DS)>0) lg="Sorted OMMs" else "Unsorted OMMs"
    if(SLY=='y')of=9 else of=0
    if(Multi==3)
    {
    amin=trunc(min(c(so-si,sob-sib,soc-sic))-1)+of
    amax=trunc(max(c(so+si,sob+sib,soc+sic))+1)+of
    }
    if(Multi==2)
    {
      amin=trunc(min(c(so-si,sob-sib))-1)+of
      amax=trunc(max(c(so+si,sob+sib))+1)+of
    }
    if(Multi==1)
    {
      amin=trunc(min(c(so-si))-1)+of
      amax=trunc(max(c(so+si))+1)+of
    }
      
    if(LabScale>0)
    {
    if(SLY=='y') plot(ii,type='n',ylim=c(amin,amax),pch=ssymb, main=Title, ylab=YAxis,xlab=xlab,cex.main=1,xaxt="n",yaxt="n")
    if(SLY=='n') plot(ii,type='n',ylim=c(amin,amax),pch=ssymb, main=Title, ylab=YAxis,xlab=xlab,cex.main=1,xaxt="n")
    #axis(1,las = 1,at=min(ii):max(ii),labels=nnn, cex.axis=1)
    if(LabScale==1 ) axis(1,las = 2,at=min(ii):max(ii),labels=nnn, cex.axis=0.5)
    if(LabScale==2 || LabScale==3 ) axis(1,las = 1,at=min(ii):max(ii),labels=nnn, cex.axis=1)
    
    if(SLY=='y')
    {
      lab=c(lab,lab,lab,lab,lab)[amin:amax] 
      axis(2,las=2,at=amin:amax,labels=lab, cex.axis=0.8)
    }
    if(Multi>=3 )
    {
      val=(soc!=0 | sic!=0) & (soc!=360)
      ii1=ii[val]
      Mc1=soc[val]+of
      Sc1=sic[val]
      colc1=ccolc[val]
      arrows(ii1,Mc1-Sc1,ii1,Mc1+Sc1, code=0,col=colc1,lwd=1)
      points(ii1,Mc1,col=colc1,pch=ssymb)
    }
    if(Multi>=2)
    {
      val=(sob!=0 | sib!=0) & (sob!=360)
      ii1=ii[val]
      Mb1=sob[val]+of
      Sb1=sib[val]
      colb1=ccolb[val]
      arrows(ii1,Mb1-Sb1,ii1,Mb1+Sb1, code=0,col=colb1,lwd=1)
      points(ii1,Mb1,col=colb1,pch=ssymb,ty='b',lty=2)
    }
    so=so+of
    arrows(ii,so-si,ii,so+si, code=0,col=ccol,lwd=1)
    points(ii,so,col=ccol,pch=ssymb,ty='b')
    
    if (Identify) rt=identify(ii,so,nnn)
    }
    return (list(ttest=resultat, median=so,StDev=si,name=nnn,SortIdx=ORD,Median1=lm1,Median2=lm2,corData=cordata,corDatab=cordatab,corcatac=cordatac))
  }


  ccenter=function(c1,c2,c3)
  {
    a=mean(c(c1[1],c2[1],c3[1]))
    b=mean(c(c1[2],c2[2],c3[2]))
    c=mean(c(c1[3],c2[3],c3[3]))
    return(c(a,b,c))
  }

  rotatePDB=function(id,frame,Residue,angle,file,IndexFile,TimeSeriesFile, mode,mirror=0) # mirror axis code: 0:none 1:c1,2:c2,3c3 to cen or gcen
  {
    target=id
    RCF=Residue
    angle=angle/180*pi
    WorkSetFile="WorkSet.txt"
    ws=read.table("WorkSet.txt",header=TRUE,sep='\t')
    Prefix="data.dcd"
    DataSet=c(paste("Set",Id1, sep='_'),table[Id1,2])
    # !!! modified dataset
    CPX1PDB  <- read.pdb(paste(DataSet[2],"/data.pdb", sep=""))
    res=SetUp1(ws$Id[target],Prefix,IndexFile,TimeSeriesFile,CPX1PDB)
    xyz=res$xyz1
    Apdb=res$CPX1PDB
    ato.inds=atom.select(Apdb, resno=RCF)
    CFFxyz=xyz[frame,ato.inds$xyz]
    Mxyz=t(matrix(CFFxyz,nrow=3,ncol=24))
    AtC1=AtomSel(Apdb,RCF,"C1")
    AtC2=AtomSel(Apdb,RCF,"C2")
    AtC3=AtomSel(Apdb,RCF,"C3")
    AtC4=AtomSel(Apdb,RCF,"C4")
    AtC6=AtomSel(Apdb,RCF,"C6")
    c1=xyz[target,AtC1$xyz]
    c2=xyz[target,AtC2$xyz]
    c3=xyz[target,AtC3$xyz]
    c4=xyz[target,AtC4$xyz]
    c6=xyz[target,AtC6$xyz]
    if(mode==2)cen=as.matrix(GCenter(c1,c2,c3,TRUE)$center)
    if(mode==1)cen=as.matrix(ccenter(c1,c2,c3))
    V1=as.vector(c1-cen)
    V2=as.vector(c2-cen)
    V3=as.vector(c3-cen)
    V46=as.vector(c4-c6)
    if(mirror==0)VNC=as.vector(cross(V1,V2))
    else{
      if(mirror==1) VNC=V1
      if(mirror==2) VNC=V2
      if(mirror==3) VNC=V3
      if(mirror==4) VNC=V46
      angle=pi
    }
    VNC=VNC/norm_vec(VNC)
    sina=sin(angle)
    cosa=cos(angle)
    if(mirror!=4)
    {
      cenm=t(matrix(rep(cen,24),nrow=3,ncol=24))
      Mxyz=Mxyz-cenm
    }
    else{
    cenm=t(matrix(rep(c4,24),nrow=3,ncol=24))
     Mxyz=Mxyz-cenm

    }
    xx=VNC[1]*VNC[1]
    xy=VNC[1]*VNC[2]
    xz=VNC[1]*VNC[3]
    yy=VNC[2]*VNC[2]
    yz=VNC[2]*VNC[3]
    zz=VNC[3]*VNC[3]
    rmat1=matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
    rmat2=matrix(c(xx,xy,xz,xy,yy,yz,xz,yz,zz),nrow=3,ncol=3)
    rmat3=matrix(c(0,VNC[3],-VNC[2],-VNC[3],0,VNC[1],VNC[2],-VNC[1],0),nrow=3,ncol=3)
    rmat=cosa*rmat1+(1-cosa)*rmat2+sina*rmat3
    rr=rmat%*%t(Mxyz)
    rr=as.vector(rr+t(cenm))
    xyz[1,ato.inds$xyz]=rr
    Apdb$atom$resid[ato.inds$atom]=rep("CFF",24) #eventually rename target CFF
    write.pdb(Apdb,file,xyz[target,])
    cc1=xyz[target,AtC1$xyz]
    cc2=xyz[target,AtC2$xyz]
    cc3=xyz[target,AtC3$xyz]
    if(mode==2)cen1=as.matrix(GCenter(cc1,cc2,cc3,TRUE)$center)
    if(mode==1)cen1=as.matrix(ccenter(cc1,cc2,cc3))
    return(list(cen=cen,cen1=cen1,VNC=VNC))
  }


  rotateFPDB=function(PDBfile,Residue, angle,file, mode,mirror=0,source=NULL) # mirror axis code: 0:none 1:c1,2:c2,3c3 to cen or gcen
  {
    RCF=Residue
    target=1
    angle=angle/180*pi
    if(is.null(source)) source="PDBfile"
    Apdb=read.pdb(source)
    xyz=Apdb$xyz
    ato.inds=atom.select(Apdb, resno=RCF)
    CFFxyz=xyz[frame,ato.inds$xyz]
    Mxyz=t(matrix(CFFxyz,nrow=3,ncol=24))
    AtC1=AtomSel(Apdb,RCF,"C1")
    AtC2=AtomSel(Apdb,RCF,"C2")
    AtC3=AtomSel(Apdb,RCF,"C3")
    AtC4=AtomSel(Apdb,RCF,"C4")
    AtC6=AtomSel(Apdb,RCF,"C6")
    c1=xyz[target,AtC1$xyz]
    c2=xyz[target,AtC2$xyz]
    c3=xyz[target,AtC3$xyz]
    c4=xyz[target,AtC4$xyz]
    c6=xyz[target,AtC6$xyz]
    if(mode==2)cen=as.matrix(GCenter(c1,c2,c3,TRUE)$center)
    if(mode==1)cen=as.matrix(ccenter(c1,c2,c3))
    V1=as.vector(c1-cen)
    V2=as.vector(c2-cen)
    V3=as.vector(c3-cen)
    V46=as.vector(c4-c6)
    if(mirror==0)VNC=as.vector(cross(V1,V2))
    else{
      if(mirror==1) VNC=V1
      if(mirror==2) VNC=V2
      if(mirror==3) VNC=V3
      if(mirror==4) VNC=V46
      angle=pi
    }
    VNC=VNC/norm_vec(VNC)
    sina=sin(angle)
    cosa=cos(angle)
    if(mirror!=4)
    {
      cenm=t(matrix(rep(cen,24),nrow=3,ncol=24))
      Mxyz=Mxyz-cenm
    }
    else{
      cenm=t(matrix(rep(c4,24),nrow=3,ncol=24))
      Mxyz=Mxyz-cenm

    }
    xx=VNC[1]*VNC[1]
    xy=VNC[1]*VNC[2]
    xz=VNC[1]*VNC[3]
    yy=VNC[2]*VNC[2]
    yz=VNC[2]*VNC[3]
    zz=VNC[3]*VNC[3]
    rmat1=matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
    rmat2=matrix(c(xx,xy,xz,xy,yy,yz,xz,yz,zz),nrow=3,ncol=3)
    rmat3=matrix(c(0,VNC[3],-VNC[2],-VNC[3],0,VNC[1],VNC[2],-VNC[1],0),nrow=3,ncol=3)
    rmat=cosa*rmat1+(1-cosa)*rmat2+sina*rmat3
    rr=rmat%*%t(Mxyz)
    rr=as.vector(rr+t(cenm))
    xyz[1,ato.inds$xyz]=rr
    Apdb$atom$resid[ato.inds$atom]=rep("CFF",24) #eventually rename target CFF
    write.pdb(Apdb,file,xyz[target,])
    cc1=xyz[target,AtC1$xyz]
    cc2=xyz[target,AtC2$xyz]
    cc3=xyz[target,AtC3$xyz]
    if(mode==2)cen1=as.matrix(GCenter(cc1,cc2,cc3,TRUE)$center)
    if(mode==1)cen1=as.matrix(ccenter(cc1,cc2,cc3))
    return(list(cen=cen,cen1=cen1,VNC=VNC))
  }


  CorrCOC=function(data)
  {
    code=c("C2(C12N3)","C2N4","N4(N9)","C8(C8)","C3(C14N7)","C3O2","O2(O13)","C1(C10N1)","O1(O11)")
    CodeMask=array(0,dim=length(data))
    up=FALSE
    down=FALSE
    dmask=rep(0,dim(data)[2])
    for(j in 1:dim(data)[2])
    {
      if(is.na(data[1,j])) next
      for(i in 2:dim(data)[1])
      {
        if(is.na(data[i,j])) break
        if(abs(data[i,j]+9-data[i-1,j])<abs(data[i,j]-data[i-1,j]))
        {
          data[i,j]=data[i,j]+9
          up=TRUE
        }
        if(abs(data[i,j]-9-data[i-1,j])<abs(data[i,j]-data[i-1,j]))
        {
          data[i,j]=data[i,j]-9
          down=TRUE
        }
      }
    }
    code=c(code,code,code)
    if(down) data[,j]=data[,j]+9
    return(list(Data=data,Code=code))
  }


  CorrCOC1=function(data)
  {
    code=c("C2(C12N3)","C2N4","N4(N9)","C8(C8)","C3(C14N7)","C3O2","O2(O13)","C1(C10N1)","O1(O11)")
    up=FALSE
    down=FALSE
      for(i in 2:length(data))
      {
        if(abs(data[i]+9-data[i-1])<abs(data[i]-data[i-1]))
        {
          data[i]=data[i]+9
          up=TRUE
        }
        if(abs(data[i]-9-data[i-1])<abs(data[i]-data[i-1]))
        {
          data[i]=data[i]-9
          down=TRUE
        }
      }
    code=c(code,code,code)
    if(down) data=data+9
    return(list(Data=data,Code=code))
  }
  
Sround=function(x)
{
  y=((x-1)+0.5)%%9+0.5
  return(y)
}

  PrepStatplot=function(matCen,matSD,matNb, TargParm,MapSize=500,AsAng=FALSE,AsCOC=FALSE,AddOfset=NULL)
  {
    if(AsCOC)
    {
      ssel=grep("Moy1",colnames(matCen))
      matCen[TargParm,ssel]=sapply(matCen[TargParm,ssel],Sround)
      ssel=grep("Moy2",colnames(matCen))
      matCen[TargParm,ssel]=sapply(matCen[TargParm,ssel],Sround)
      ssel=grep("Moy3",colnames(matCen))
      matCen[TargParm,ssel]=sapply(matCen[TargParm,ssel],Sround)
    }
    if(!is.null(AddOfset))
    {
      ssel=grep("Moy1",colnames(matCen))
      matCen[TargParm,ssel]=matCen[TargParm,ssel]+AddOfset
      ssel=grep("Moy2",colnames(matCen))
      matCen[TargParm,ssel]=matCen[TargParm,ssel]+AddOfset
      ssel=grep("Moy3",colnames(matCen))
      matCen[TargParm,ssel]=matCen[TargParm,ssel]+AddOfset
    }
    dd=dim(matCen)[2]/3
    m=array(NA,dim=c(3,dd))
    for(j  in 1:dd)
    {
      u=(j-1)*3+1
      v1=matNb[TargParm,u]!=0
      v2=matNb[TargParm,u+1]!=0
      v3=matNb[TargParm,u+2]!=0
      if(v1) m[1,j]=matCen[TargParm,u]
      if(v2) m[2,j]=matCen[TargParm,u+1]
      if(v3) m[3,j]=matCen[TargParm,u+2]
    }
    if(AsAng)
    {
      rr=CorAngl(m[1,], m[2,],m[3,])
      m[1,]=rr[,"ang"]
      m[2,]=rr[,"ang1"]
      m[3,]=rr[,"ang2"]
    }
    
    nn=dim(matNb)[2]/3
    val=array(NA,dim=c(MapSize,3*nn))
    for(j  in 1:nn)
    {
      u=(j-1)*3+1
      tp=sum(matNb[TargParm,u:(j*3)])
      n1=trunc(matNb[TargParm,u]*MapSize/tp)
      n2=trunc(matNb[TargParm,u+1]*MapSize/tp)
      n3=trunc(matNb[TargParm,u+2]*MapSize/tp)
      val[1:n1,j]=rnorm(n1,m[1,j],matSD[TargParm,u])
      if(AsCOC)
      {
        val[1:n1,j]=sapply(val[1:n1,j],Sround)
        if(!is.null(AddOfset)) val[1:n1,j]=val[1:n1,j]+AddOfset[j]
      }
      if(n2!=0)
      {
        val[1:n2,j+nn]=rnorm(n2,m[2,j],matSD[TargParm,u+1])
        if(AsCOC)
        {
          val[1:n2,j+nn]=sapply(val[1:n2,j+nn],Sround)
          if(!is.null(AddOfset)) val[1:n2,j+nn]=val[1:n2,j+nn]+AddOfset[j]
        }
      }
      if(n3!=0)
      {
        val[1:n3,j+2*nn]=rnorm(n3,m[3,j],matSD[TargParm,u+2])
        if(AsCOC)
        {
          val[1:n3,j+2*nn]=sapply(val[1:n3,j+2*nn],Sround)
          if(!is.null(AddOfset)) val[1:n3,j+2*nn]=val[1:n3,j+2*nn]+AddOfset[j]
        }
      }
    }
    return(val)
  }

  CorAngl=function(ang, ang1,ang2)
  {
    
    for(i in 2:(length(ang)))
    {
      shift=0
      if(abs(ang[i]+360-ang[i-1])<abs(ang[i]-ang[i-1])) shift=360
      if(abs(ang[i]-360-ang[i-1])<abs(ang[i]-ang[i-1])) shift=-360
      for(j in i:length(ang))
      {
      ang[j]=ang[j]+shift
      ang1[j]=ang1[j]+shift
      ang2[j]=ang2[j]+shift
      }
    }
    repeat
    {
      shift=0
      if(max(ang)>=360 ) shift=-360 else{
      if(min(ang)<=-360) shift=360
      }
      if(shift==0) break
      for(j in 1:length(ang))
      {
        ang[j]=ang[j]+shift
        ang1[j]=ang1[j]+shift
        ang2[j]=ang2[j]+shift
      }
    }
    return(data.frame(ang=ang,ang1=ang1,ang2=ang2))
  }
  
  RainbowPlot=function(val,Limits=NA, Scale=NA,Xlab,Ylab,MainLab, Item,Tab, mark=NULL, Ofset=0,AddOfset=NULL)
  {
    AddOfset=AddOfset[dim(val)[2]/3]
    transparency=0.15
    if(is.na(Limits)[1]) yll=c(-1000,1000) else yll=Limits
    dd=dim(val)
    val=sapply(val,max,yll[1])
    val=array(sapply(val,min,yll[2]),dim=dd)
    if(is.na(Scale)[1]) ylim=c(max(Scale[1],min(val)), min(Scale[2],max(val))) else ylim=Scale
    nn=dim(val)[2]/3
    mm=array(0,dim=nn)
    for(i in 1:nn)
    {
      v1=val[!is.na(val[,i]),i]+Ofset
      v2=val[!is.na(val[,i+nn]),i+nn]+Ofset
      v3=val[!is.na(val[,i+2*nn]),i+2*nn]+Ofset
      mm[i]=mean(v1)
      x1=rep((i-1)*30,length(v1))
      x2=rep((i-1)*30,length(v2))
      x3=rep((i-1)*30,length(v3))
      if(i==1)
      {
        if(is.null(mark)) plot(x3,v3,xlim=c(0,330),ylim=ylim,pch=16,xlab=Xlab,ylab=Ylab,main=MainLab,col=rgb(0,0,1,transparency))
        else plot(x3,v3,xlim=c(0,330),ylim=ylim,pch=16,xlab=Xlab,ylab=Ylab,main=MainLab,col=rgb(0,0,1,transparency), yaxt="n")
        if(!is.null(mark))
        {
          ll=1:Scale[2]
          axis(2,at =ll,labels=mark[ll], las=2,cex.axis=0.6)
          #axis(1)
        }

      }else points(x1,v1, col=rgb(0,0,1,transparency),pch=16)
      grid()
      points(x2,v2, col=rgb(0,1,0,transparency),pch=1)
      points(x1,v1, col=rgb(1,0,0,transparency),pch=1)
    }
    InDat=Tab[Item,grep("InitVal",colnames(Tab))]
    InDat=CorrCOC1(InDat)$Data+Ofset+AddOfset
    lines(seq(0,30*(nn-1),by=30),InDat,lwd=2, col="orange",lty="dashed",type="b")
    lines(seq(0,30*(nn-1),by=30),mm,lwd=2,type='b')
  }
  
  RainbowPlotSP=function(val,Limits=NA, Scale=NA,Xlab,Ylab,MainLab, Item,Tab, mark=NULL, Ofset=0,LocOfset=NULL, mode=0)
  {
    par(mfrow=c(1,1))
    if(is.null(LocOfset)) LocOfset=rep(0,(dim(val)[2]/3))
    LocOfset=LocOfset[1:(dim(val)[2]/3)]
    transparency=1
    if(is.na(Limits)[1]) yll=c(-1000,1000) else yll=Limits
    dd=dim(val)
    val=sapply(val,max,yll[1])
    val=array(sapply(val,min,yll[2]),dim=dd)
    nn=dim(val)[2]/3
    mm=array(0,dim=nn)
    xx=1:nn
    Th=c(100,50,25,12,6,0)
    PP=c(16,20,1,2,4)
    InDat=Tab[Item,grep("InitVal",colnames(Tab))]
    InDat=CorrCOC1(InDat)$Data
    lcte=list()
    lpc=list()
    lcx=list()
    lcol=list()
    for(i in 1:nn)
    {
      dd=val[!is.na(val[,i]),i]
      tt=table(factor(dd, levels = c(1:9)))
      fr=as.numeric(tt)
      cte=as.numeric(names(tt))
      sfr=sort(fr,decreasing = T, index.return = T)
      fr=sfr$x
      cte=cte[sfr$ix]
      cte=cte[fr>0]
      fr=fr[fr>0]
      fr=fr/sum(fr)*100
      if(i==1) Ref=cte[1] else Ref=mm[i-1]
      Ref=InDat[i]
      col=c()
      pc=c()
      cx=c()
      
      for(j in 1:length(fr))
      {
       s=seq(cte[j],cte[j]+44, by=9)-27
       if(mode==0 || (i==1)) n= which(abs(s-Ref)==min(abs(s-Ref)))
       if(mode==1 && (i>1))
        {
         n=which(abs(s-LM)==min(abs(s-LM)))
       }
        cte[j]=s[n]
        if(j==1) Um=s[n]
        if(fr[j]>50)
        {
          col=c(col,rgb(0,0,0))
          pc=c(pc,16)
          cx=c(cx,3) 
        }
        if(fr[j]>25 && fr[j]<=50) 
        {
          col=c(col,rgb(0.5,0.5,0.5))
          pc=c(pc,16)
          cx=c(cx,3)
        }
        if(fr[j]>12 && fr[j]<=25) 
        {
          col=c(col,rgb(0.8,0.8,0.8))
          pc=c(pc,16)
          cx=c(cx,2)
        }
        if(fr[j]>6 && fr[j]<=12) 
        {
          col=c(col,rgb(0.8,0.8,0.8))
          pc=c(pc,16)
          cx=c(cx,1)
        }
        if(fr[j]>0 && fr[j]<=6) 
        {
          col=c(col,rgb(1,1,1))
          pc=c(pc,4)
          cx=c(cx,1)
        }
      }
      LM=Um
      mm[i]=cte[1]
      cte=cte+Ofset+LocOfset[i]
#      cte[cte<1]=cte[cte<1]+9
      InDat=InDat+Ofset+LocOfset
      mm=mm+Ofset+LocOfset
      lcte[[i]]=cte
      lpc[[i]]=pc
      lcx[[i]]=cx
      lcol[[i]]=col
    }
    
    mmax=max(c(unlist(lcte),mm,InDat))
    mmin=min(c(unlist(lcte),mm,InDat))
    ylim=c(mmin-2,mmax+2)
    if(!is.null(Scale))
    {
      ylim[1]=min(Scale[1],ylim[1])
      ylim[2]=max(Scale[2],ylim[2])
    }
    ll=ylim[1]:ylim[2]
    for(i in 1:nn)
    {
      cte=lcte[[i]]
      pc=lpc[[i]]
      cx=lcx[[i]]
      col=lcol[[i]]
      xx=rep(1+(i-1)*30,length(cte))
      par(mar=c(5.1, 5.5, 4.1, 2.1))
      if(i==1) 
      {
        plot(xx,cte,xlim=c(0,330),ylim=ylim,pch=pc,cex=cx,xlab=Xlab,ylab=Ylab,main=MainLab,col=col, yaxt="n",cex.axis=1.2,cex.lab=1.5)
      } else
        {
          points(xx,cte, col=col,pch=pc,cex=cx)
          points(xx,cte, col=rgb(0,0,0),pch=1,cex=cx)
        }
    }
    
    lab=(ll+Ofset-1)%%9+1
    axis(2,at =ll,labels=mark[lab], las=2,cex=cx,cex.axis=1)
    grid()
    lines(seq(0,30*(nn-1),by=30),InDat,lwd=2, col="red",lty="dashed",type="b")
    lines(seq(0,30*(nn-1),by=30),mm,lwd=2,type='l')
    par(mfrow=c(1,2))
  }


  ViolonPlot=function(val, Limits=NA,Scale=NA,Xlab,Ylab,MainLab, Item, Tab, mark=NULL, Ofset=0,AddOfset=NULL)
  {
    AddOfset=AddOfset[dim(val)[2]/3]
    transparency=0.15
    nn=dim(val)[2]/3
    if(is.na(Limits)[1]) ylim=c(-1000,1000) else ylim=Limits
    dd=dim(val)
    val=sapply(val,max,ylim[1])
    val=array(sapply(val,min,ylim[2]),dim=dd)
    u=as.vector(val)
    u=u[u!=0]
    if(is.na(Scale)[1])
    {
      ylim=c(min(u),max(u))
    }else ylim=Scale
    plot(c(0:11),ylim=ylim,col="white",type="n",axes=F,xlab=Xlab,ylab=Ylab,main=MainLab)
    grid()
    if(!is.null(mark))
    {
      ll=1:Scale[2]
      axis(2,at =ll,labels=mark[ll], las=2,cex.axis=0.6)
      #axis(1)
    }else axis(2)
   # mtext("Paramètre",side=2,col="black",line=2.5)
    grid()
    axis(1,at=1:12,labels=as.character(seq(0,330,by=30)))
    mm=array(0,dim=12)
    for(i in 1:nn)
    {
      v1=val[!is.na(val[,i]),i]+Ofset
      v2=val[!is.na(val[,i+nn]),i+nn]+Ofset
      v3=val[!is.na(val[,i+2*nn]),i+2*nn]+Ofset
      x1=rep(i,length(v1))
      x2=rep(i,length(v2))
      x3=rep(i,length(v3))
      if(length(v3)>0) vioplot(v3,at=i,horizontal=F,col=rgb(0.8,0.8,0.8,0.5),add=T)
      if(length(v2)>0)vioplot(v2,at=i,horizontal=F,col=rgb(0.8,0.8,0.8,0.5),add=T)
      if(length(v1)>0)vioplot(v1,at=i,horizontal=F,col=rgb(1,0.5,0.5,1),add=T)
      mm[i]=mean(v1)
    }
    InDat=Tab[Item,grep("InitVal",colnames(Tab))]
    InDat=CorrCOC1(InDat)$Data+Ofset+AddOfset
    lines(InDat,lwd=2, col="orange",lty="dashed",type="b")
    lines(mm,lwd=2,type='b')

  }
  
  GetPlotParms=function(ParmFile,Dataset,CffId,Parm)
  {
    PlotParms=read.table(ParmFile, sep='\t', header = TRUE)
    n1=which(PlotParms[,2]==Dataset & PlotParms[,3]==Parm & PlotParms[,4]==CffId)
    if(length(n1)>0) n2=which(PlotParms[n1,1]==max(PlotParms[n1,1])) else n2=NULL
    if(!is.null(n2) && length(n1[n2])>0) return(PlotParms[n1[n2],]) else{
      n1=which(PlotParms[,3]==Parm)
      if(length(n1)>0) n2=which(PlotParms[n1,1]==max(PlotParms[n1,1])) else n2=NULL
      if(!is.null(n2) && length(n2)>0) return(PlotParms[n1[n2],]) else return(NULL)
    }
  }
  
  SetPlotParms=function(ParmFile,Dataset,CffId,Parm,MyParms)
  {
    Plotparms=read.table(ParmFile, sep='\t', header = TRUE)
    n1=which(Plotparms[,2]==Dataset & Plotparms[,3]==Parm & Plotparms[,4]==CffId)
    if(length(n1)>0) n2=which(Plotparms[n1,1]==max(Plotparms[n1,1])) else n2=NULL
    
    if(!is.null(n2) && length(n1[n2])>0)
    {
      u=n1[n2][1]
      if(MyParms[1]==0) MyParms[1]=max(Plotparms[,1])+1
      Plotparms[u,]=MyParms
      write.table(Plotparms,file=ParmFile, sep='\t')
    }else{
    x=max(Plotparms[,1])+1
    MyParms[1]=x
    Plotparms=rbind(Plotparms,MyParms)
    }
    write.table(Plotparms,file=ParmFile, sep='\t')
  }
  
  

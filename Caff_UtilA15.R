#Accessory functions ----------------------------------------------------------------------------
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

# Rotates a point P in space around a line of director vector V crossing point R with angle theta
#to rotate a vector P set R to c(0,0,0)
rotate_in_space <- function(P, V, R, theta) {
    theta=theta/180*pi
    P_translated <- P - R  # Translate P to the origin
    VecN <- V / norm_vec(V)  #Normalize the rotation axis
    ST=sin(theta)
    CT=cos(theta)
    K <- matrix(c(
      CT + VecN[1]^2 * (1 - CT),
      VecN[1] * VecN[2] * (1 - CT) - VecN[3] * ST,
      VecN[1] * VecN[3] * (1 - CT) + VecN[2] * ST,
      VecN[1] * VecN[2] * (1 - CT) + VecN[3] * ST,
      CT + VecN[2]^2 * (1 - CT),
      VecN[2] * VecN[3] * (1 - CT) - VecN[1] * ST,
      VecN[1] * VecN[3] * (1 - CT) - VecN[2] * ST,
      VecN[2] * VecN[3] * (1 - CT) + VecN[1] * ST,
      CT + VecN[3]^2 * (1 - CT)
    ), nrow = 3, ncol = 3)
    P_rotated <- K %*% P_translated # Rotate the translated point
    P_final <- P_rotated + R #Translate the rotated point back to its original position
    return(P_final)
  }

#calculate the length of a vector  
Norm_vec=function(Vect){return(sqrt(dot(vect,vect)))}

#Calculate angle between 2 vectors (-180-+180 or 0-360)
GetAngle=function(v1,v2,Code,AngRec=NULL,Adjust=FALSE) #code negatif : set, code =0 read angle (as acos) other :apply
{
  v1N=v1/norm_vec(v1)
  v2N=v2/norm_vec(v2)
  CN=dot(v1N,v2N)
  ca=acos(CN)/pi*180
  PN=cross(v1N,v2N)
  if(Code==0) return(ca)
  if(Code<0)
  {
    CC=PN
    AngRec[-Code,]=PN
    Code=-Code
  } else CC=AngRec[Code,]
 
  WW=PN%*%CC
  AA=ca*sign(WW)
  if(Adjust)
  {
  PN1=CC+0.1*PN*sign(WW)
  AngRec[Code,]=PN1/norm_vec(PN1)
  }
  AngRec<<- AngRec
  return (AA)
}

GetAngle2=function(v1,v2,RNorm=NULL) 
{
  v1N=v1/norm_vec(v1)
  v2N=v2/norm_vec(v2)
  CN=dot(v1N,v2N)
  ca=acos(CN)/pi*180
  PN=cross(v1N,v2N)
  PN=PN/norm_vec(PN)
  PN1=RNorm$VV
  if(!is.null(RNorm))
  {
    CC=PN%*%RNorm$VV
    if(CC<0) ca=-ca
  } else CC=1
  return (list(AA=ca,VV=PN1,SS=CC))
}

NormAngl=function(data)
{

  if(length(data)>1)
  {
    for(i in 2:length(data))
    {
        set=c(data[i],data[i]-720,data[i]-360,data[i]+360,data[i]+720)
        rset=abs(set-data[i-1])
        ss=which(rset==min(rset))
        if(length(ss)==1)  dd=set[ss] else dd=min(set[ss])
        data[i]=dd
    }
  }
  return (data)
}
NormTo180=function(ang)
{
  ang[ang>90]=180-ang[ang>90]
  ang[ang<(-90)]=-180-ang[ang<(-90)]
  return (ang)
}

CentAngl=function(Ang, AllPos=FALSE)
{
  Ang1=Ang
  if(AllPos)
  {
    ss=Ang<0
    Ang1[ss]=360-Ang[ss]
  }
  if(max(Ang1)>360) Ang1=Ang1-360
  if(min(Ang1)<(-360)) Ang1=Ang1+360
  return (Ang1)
}

NoNegAngl=function(Ang)
{
  ss=(Ang<0)
  Ang[ss]=360-Ang[ss]
  Ang=Ang%%360
  return (Ang)
}


#rr (point): intersect of line (Axe,Ori) rotated by Ang1 
# with circle (Pt, Dist)
# Ty=TRUE larger Ori-rr distance; FALSE smaller
FIndCV=function(Axe,Ori,Pt,Dist, Angl,Ty) 
{
  Axe2=rotate_in_space(Pt, Axe, Ori, Ang)-Ori
  Axe2=as.vector(Axe2/norm_vec(Axe2))
  H1=dot(Axe2,Pt-Ori)/dot(Axe2,Axe2)*Axe2+Ori
  N1=norm_vec(H1-Pt)
  L=sqrt(Dist^2-N1^2)
  N2=norm_vec(Ori-H1)
  if(Ty) rr=Ori+(N1+N2)*Axe2 else rr=Ori+(N1-N2)*Axe2
  return (rr)
}




# Main functions-----------------------------------------------------------------------------------

GetM3=function(C1,C2,C3,mode=1) #mode: 1:center of circle, 2:mean
{
  tt=as.vector(unlist(GCenter(C1,C2,C3,TRUE)$center))
  RCFF=norm_vec(tt-C2)
  if(mode==1) cen=tt else cen=(C1+C2+C3)/3
  return (list(cen=cen,RCFF=RCFF))
}

SetNP1=function(N1,N2,N3,N4)
{
  points=rbind(N1,N2,N3,N4)
  mean = colMeans(points)
  centered_points = points - rbind(mean, mean, mean, mean)
  NP1=cross(N1-N3,N2-N4)
  try({
  cov_mat = cov(centered_points)
  eigen_vals = eigen(cov_mat)$values
  min_index = which.min(eigen_vals)
  NP1 = eigen(cov_mat)$vectors[, min_index]
  sc=cross(N1-N3,N2-N4)
  if(sc%*%NP1<0)  NP1=-NP1
  })
  Np1=NP1/norm_vec(NP1)
  return (data.frame(NP1=NP1,CC=mean))
}

# evaluate normal vector and center to porphirin nitrogen (NP1,CC1, err: fit error) and CFF (NP2,CC2)
SetNP2=function(C1,C2,C3) 
{
  NP2=cross(C1-C2,C1-C3)
  NP2=NP2/norm_vec(NP2)
  CC2=(C1+C2+C3)/3
  return (data.frame(NP2=NP2,CC2=CC2)) 
}

# evaluate normal vector to caffein plane from angle A1 and A3 
SetNormVect2=function(N1,N2,N3, N4,A0,A1,NP2) 
{
  NP1=SetNP1(N1,N2,N3,N4)$NP1
  ITC=as.vector(rotate_in_space (N4-N2, NP1, c(0,0,0),-A0))
  NP2=as.vector(rotate_in_space (NP1, ITC, c(0,0,0), -A1))
  return (data.frame(NP1=NP1,NP2=NP2,LI=ITC))
}

GetM1=function(N1,N2,N3,N4, Fe,NP1,DHO=1.66, mode=1)
{
  points=rbind(N1,N2,N3,N4)
  mean = colMeans(points)
  centered_points <- points - rbind(mean, mean, mean, mean)
  ori = points[1, ] - c(centered_points[1, ] %*% NP1) * NP1
  dist=array(0,dim=4)
  for(i in 1:4)dist[i]=abs(dot(NP1,(points[i,]-ori)))/sqrt(dot(NP1,NP1))
  err=sum(dist)/4
  O1=Fe+NP1/norm_vec(NP1)*DHO
  if (mode==1) M1=ori
  if(mode==2) M1 =Fe
  if(mode==3) M1=O1
  return(list(Parm=data.frame(M1=M1,O1=O1,NNNN=ori),Err=err))
}
# NP1 (et )normal vector oriented toward CFF) and CC1 (center) define porphirin /heme plane centered on mean(N) or Fe or O 
# NP2 (et )normal vector oriented toward He) and CC21 (center) define CFF plane 
# Fe heme iron
# D4 distance between M2 and M3
# RG distance between M3 and M4
# A3 angle M2-M3/M2-MP1
# A4 angle M3M2 M3-M4
# DHO Iron to active oxygen distance
#N1,N2,N3,N4 porphirin nitrogens
#CEN : caffein cen (center od C1,C2,C3 circle)  and C2: caffeine C2 (paper numbering)
#returned values
#M1P: intersection of normal projection of CC1 on intersect line of porphirin and caffeine planes
#M2 : normal projection on caffeine plane of M1 (either porphyrin N1-N4 cen, Fe or active O) 
#A0: angle Porphyrin 'N2-N4' with 'MP1-ori'
#A1: angle of Porphyrin and caffein planes
#A2: angle caffein 'C2-cen' with 'M2-cen'
#A3: angle of 'MP1-M2' with 'M2-cen'
#A4: angle of 'MP1-M2' with 'Mp1-cen'
#O1 active oxygen
#d2:Distance from M1 to its projection on caffeine plane (M2)
#d4: Distance from caffein cen (M4) to M1 projection on caffeine plane (M2)
#d5: Distance from caffein cen (M4) to C2
#d6: distance of MP 1to caffein cen (M4).
#local coordinate system include: A0,A1,A2, A3, A4, D2,D4 (7 parameters)
#Mode:  Reference plane position definition for Heme: 1: porphyrin N1-N4 cen, 2: Fe, 3: Active O
# Reference plane angle definition: N1,N2,N3,N4 average plane

SetMP1=function(LI,NP1,NP2,M1,D2,Way)
{
  M2=M1-D2*NP2*Way 
  c1=dot(LI,M1)
  c2=dot(NP2,M2)
  c3=dot(NP1,M1)
  C=matrix(c(c1,c2,c3),ncol=1)
  A=rbind(LI,NP2,NP1)
  MP1=as.vector(solve(A,C))
  return(data.frame(MP1=MP1,M2=M2))  
}


GetSysAngleDist <- function(M1,M2,M3,C2,NP1,NP2,N2,N4, Way,Face="L") # 1:initialize 2:use
{
  Rec=NULL
  if(!is.null(Rec)) AngRec=Rec
  LI=SetLI(NP1,NP2,N4,N2) #orienté comme n4-N2
  c1=dot(LI,M1)
  c2=dot(NP2,M3)
  c3=dot(NP1,M1)
  C=matrix(c(c1,c2,c3),ncol=1)
  A=rbind(LI,NP2,NP1)
  MP1=as.vector(solve(A,C)) 
  if(abs(NP1%*%(MP1-M1))>1e-5 || abs(NP2%*%(MP1-M3))>1e-5) STOP("MP1 error UtilA_2755")
  AX0=as.vector((MP1-M1)/norm_vec(MP1-M1))
  AX1=as.vector((M2-MP1)/norm_vec(M2-MP1))
  AX2=as.vector((M3-MP1)/norm_vec(M3-MP1))
  AX3=as.vector((M2-M3)/norm_vec(M2-M3))
  AX4=as.vector((C2-M3)/norm_vec(C2-M3))
  
  AngRec[1,]=NP1
  A0=GetAngle(N4-N2,LI,1,AngRec)
  AngRec[2,]=LI
  A1=GetAngle(NP1,NP2,2, AngRec)
  AngRec[3,]=NP2
  A2=GetAngle(AX3,AX4,3,AngRec)
  AngRec[4,]=NP2
  if(Way==1) A3=GetAngle(AX1,-AX3, 4,AngRec) else A3=GetAngle(-AX1,-AX3,4,AngRec)
  AngRec[5,]=NP2
  A4=Way*GetAngle(AX1,AX2,5, AngRec)
  d2=norm_vec(M1-M2)
  d4=norm_vec(M2-M3)
  d5=norm_vec(M3-C2)
  d6=norm_vec(MP1-M3)
  AngRec<<- AngRec
  return (list(MP1=MP1, M2=M2,LI=LI,Parm=data.frame(A0=A0,A1=A1,A2=A2,A3=A3,A4=A4,D2=d2,D4=d4,D5=d5,D6=d6)))
}


GetCoordFromCoreParm=function(Fe,N1,N2,N3,N4,D2,D4,A0,A1,A2,A3,Way, mode=1,rec=NULL)
{
  if(!is.null(Rec)) AngRec=Rec
  DHO=1.66
  RCFF=3.233081
  F2=GetM1(N1,N2,N3,N4, Fe,NP1,DHO, mode)
  M1=F2$Parm$M1
  O1=F2$Parm$O1
  F3=SetNormVect2(N1,N2,N3, N4,A0,A1)
  NP1=F3$NP1
  NP2=F3$NP2
  F4=GetParmFromAngle(M1,NP1, NP2,D2, D4,A2,A3,N2,N4, Way,RCFF,Rec)
  MP1=F4$MP1
  M2=F4$M2
  M3=F4$M3
  M4=F4$M4
  A4=F4$A4[1]
  A6=F4$A6[1]
  return(list(Coord=data.frame(M1=M1, O1=O1, MP1=MP1, M2=M2, M3=M3, M4=M4, NP1=NP1, NP2=NP2), A4=A4, A6=A6, F5=F5))
}

SetNormData=function(index,xyz, Target,AT,Way)
{
  if(!is.null(ForceTarget) )Target=ForceTarget
  Validate=TRUE
  AdjustA3A2=TRUE
  c1=xyz[index,AT$AtC1$xyz]
  c2=xyz[index,AT$AtC2$xyz]
  c3=xyz[index,AT$AtC3$xyz]
  c8=xyz[index,AT$AtC8$xyz]
  na=xyz[index,AT$AtNA$xyz]
  nc=xyz[index,AT$AtNC$xyz]
  nb=xyz[index,AT$AtNB$xyz]
  nd=xyz[index,AT$AtND$xyz]
  fe=xyz[index,AT$AtFE$xyz]
  
  DHO=1.66 #distance Fe-O
  
  if(Target=="FE") mode1=2
  if(Target=="OXY") mode1=3
  NP1=SetNP1(na,nb,nc,nd)$NP1
  F2=GetM1(na,nb,nc,nd, fe,NP1,DHO, mode1)
  M1=F2$Parm$M1
  O1=F2$Parm$O1
  
  M3=GetM3(c1,c2,c3,1)$cen
  RCFF=3.233081 
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
  if(abs(NP2%*%(M21-M3))<abs(NP2%*%(M22-M3))) M2=M21 else M2=M22
  D2=norm_vec(M1-M2)
  D4=norm_vec(M2-M3)
  
  F3=GetSysAngleDist(M1,M2,M3,c2,NP1,NP2,nb,nd,Way)
  MP1=F3$MP1
  dMP1M1=norm_vec(MP1-M1)
  dMP1M2=norm_vec(MP1-M2)
  if(abs(NP2%*%(MP1-M3))>1e-5 || abs(NP1%*%(MP1-M1))>1e-5 || abs(abs((MP1-M1)%*%(MP1-M2)/dMP1M1/dMP1M2)-abs(NP1%*%NP2))>1e-5)
  {
    stop("Error MP1 UtilA 312-")  
  }
  A0=F3$Parm$A0 # Angle of N4->N2 vector with intersect of CFF and heme plane
  A1=F3$Parm$A1 # #Angle of CFF and heme planes
  A2=F3$Parm$A2 # Angle of C2->cen with M2->cen. M2 (Proj. Fe/Oxy on CFF plane)
  A3=F3$Parm$A3 # Angle M2->Mp1 with M2->cen. MP1 (Proj. of Fe/Oxy on intersect of Heme and CFF planes)
  A4=F3$Parm$A4 # Angle M2->Mp1 with cen->Mp1. MP1 (Proj. of Fe/Oxy on intersect of Heme and CFF planes)
  D4=F3$Parm$D4
  D2=sign(Way)*F3$Parm$D2
  D5=norm_vec(MP1-M2)
  if(Target=="FE") mode=2 else mode=3
  
  if(Validate)
  {
    CCFC=GetBaseDescFromCoord(fe,na,nb,nc,nd,D2,D4,A0,A1,A2,A3,Way1,mode, Face)
    m1=CCFC$Coord$M1
    m2=CCFC$Coord$M2
    m3=CCFC$Coord$M3
    m4=CCFC$Coord$M4
    mp1=CCFC$Coord$MP1
    np2=CCFC$Coord$NP2
    err=max(norm_vec(m1-M1),norm_vec(m2-M2),norm_vec(m3-M3),norm_vec(m4-c2))
    if(err>0.2)
    {
      stop("Invalid Reverse coordinates")
    }
  }
  
  df=data.frame(A0=A0,A1=A1,A3=A3,D2=D2,D4=D4,A2=A2,D5=D5)
  return (df)
}

GetBaseDescFromCoord=function(Fe, N1,N2,N3,N4,D2,D4,A0,A1,A2,A3,Way, mode,Face)
{
  Rec=NULL
  if(!is.null(Rec)) AngRec=Rec
  DHO=1.66
  RCFF=3.233081
  NP1=SetNP1(N1,N2,N3,N4)$NP1
  F2=GetM1(N1,N2,N3,N4, Fe,NP1,DHO, mode)
  M1=F2$Parm$M1
  O1=F2$Parm$O1
  F3=SetNormVect2(N1,N2,N3, N4,A0,A1)
  NP1=F3$NP1
  NP2=F3$NP2
  F4=GetParmFromAngle(M1,NP1, NP2,D2, D4,A2,A3,N2,N4, Way,RCFF,Rec)
  
  MP1=F4$MP1
  M2=F4$M2
  M3=F4$M3
  M4=F4$M4
  A4=F4$A4[1]
  A6=F4$A6[1]
  F5=NULL
  return(list(Coord=data.frame(M1=M1, O1=O1, MP1=MP1, M2=M2, M3=M3, M4=M4, NP1=NP1, NP2=NP2), A4=A4, A6=A6, F5=F5))
}

CalibrateCFFAngle=function(C1,O1,C2,N4,C8,C3,O2)
{
  res=c(0,0,0,0,0,0,0,0,0)
  Center=as.vector(unlist(GCenter(C1,C2,C3,TRUE)$center))
  VN=cross((C1-C2),(C1-C3))
  VN=VN/norm_vec(VN)
  AngRec[6,]=VN
  res[3]=GetAngle(C2-Center,N4-Center,6, AngRec)
  res[2]=res[3]/2
  res[4]=GetAngle(C2-Center,C8-Center,6,AngRec)
  res[5]=GetAngle(C2-Center,C3-Center,6,AngRec)
  res[7]=360+GetAngle(C2-Center,O2-Center,6,AngRec)
  res[6]=(res[5]+res[7])/2
  res[8]=360+GetAngle(C2-Center,C1-Center,6,AngRec)
  res[9]=360+GetAngle(C2-Center,O1-Center,6,AngRec)
  return(res)
}

GetProjCood=function(C1,C2,C3,C8)
{
  SetAt=c("C2","C2N4","N4","C8","C3","C3O2","O2","C1","O1","C2")
  SetAn=c(0,30.613,61.227,90.297,139.83, 172.51,205.18,260.17,309.67)
 # SetAn=c(0,25.67,62.94,93.12,220.92, 170.66,210.08,258.52,310.41)
  GC=GCenter(C1,C2,C3,TRUE)
  Center=unlist(GC$center)
  Rayon=unlist(GC$Rayon)
  VN=cross((C1-C2),(C1-C3))
  VN=VN/norm_vec(VN) # normal vector to CFF
  SetCo=array(c(C2,c(0,0,0),c(0,0,0),C8,C3,c(0,0,0),c(0,0,0),C1,c(0,0,0)),dim=c(3,9))
  id=c(2,3,6,7,9)
  for(i in 1:5) SetCo[,id[i]]=rotate_in_space(C2, VN, Center, SetAn[id[i]])
  return(list(Cen=Center,Ray=Rayon,VN=VN,ASC=data.frame(C2=SetCo[1],C2N4=SetCo[2],N4=SetCo[3],C8=SetCo[4],C3=SetCo[5],
                    C3O2=SetCo[6],O2=SetCo[7],C1=SetCo[8],O1=SetCo[9])))
}

GetIminFromDist=function(C1,C2,C3,C8, At,Init=0)
{
  SetAt=c("C2","C2N4","N4","C8","C3","C3O2","O2","C1","O1","C2")
  TY=GetProjCood(C1,C2,C3,C8)
  Center=TY$Cen
  Rayon=TY$Ray
  VN=TY$VN
  C=-VN%*%Center #CFF plane : VN%*%X+C=0
  k=as.numeric(VN%*%(At-Center)) # distance P plan caffeine
  P=At-k*VN # projection of At on CFF plane
  DP=array(0,dim=9)
  for (i in 1:9) DP[i]=norm_vec(P-TY$ASC[i])
  imin=which(DP==min(DP))
  return(imin)
}

IDCffProject=function(C1,C2,C3,C8, At,Init=0,Rec=NULL)
  {
  if(!is.null(Rec)) AngRec=Rec
  GC=GCenter(C1,C2,C3,TRUE)
  Center=as.vector(unlist(GC$center))
  Rayon=unlist(GC$Rayon)
  VN=cross(C1-C2,C1-C3)
  VN=VN/norm_vec(VN) 
  C=-VN%*%Center #CFF plane : VN%*%X+C=0
  k=as.numeric(VN%*%(At-Center)) # distance P plan caffeine
  P=At-k*VN # projection of At on CFF plane
  DisAP=norm_vec(P-At)
  DisPC=norm_vec(P-Center)
  PC=Center+(P-Center)/DisPC*Rayon # closest point from At on C1,C2,C3 circumscribed circle
  PC=Center+(P-Center)/DisPC*Rayon # closest point from At on C1,C2,C3 circumscribed circle
  AngRec[7,]=VN
  Ang=GetAngle(PC-Center,C2-Center,7,AngRec)
  
  Set9=c("C2","C2N4","N4","C8","C3","C3O2","O2","C1","O1","C2") # ref for imin
  Set9A=SetAn=c(0,30.613,61.227,90.297,139.83, 172.51,205.18,260.17,309.67,360)
  Idl=array(0,dim=9)
  Angc=(-Ang)%%360
  Angc[Angc<0]=360-Angc[Angc<0]
  
  for(i in 1:9)
  {
    if(Ang>=Set9A[i] && Ang<=Set9A[i+1])
    {
      Idl[i]=(Set9A[i+1]-Angc)/(Set9A[i+1]-Set9A[i])
      if(i<=8) Idl[i+1]=1-Idl[i] else Idl[1]=1-Idl[i]
      break
    }
  }
  imin=which(abs(Set9A-Angc)==min(abs(Set9A-Angc)))
  imin=((imin-1)%%9)+1
  
  #check
  rx=norm_vec(C2-Center)
  cc=array(0,dim=c(9,3))
  dd=array(0,dim=9)
  Angd=array(0,dim=9)
  for (i in 1:9)
  {
    cc[i,]=rotate_in_space(C2, VN, Center,-Set9A[i])
    dd[i]=norm_vec(cc[i,]-PC)
    AngRec[13,]=VN
    Angd[i]=GetAngle(cc[i,]-Center,PC-Center,13,AngRec)
  }
  imind=which(dd==min(dd))[1]
  if(imind!=imin)
    {
    cat('internal error Caff_util:609')
    tk_messageBox(type="ok",'Inconsistancy in orientation/angle calculations')
  }
    
  STA=t(array(c(C1,C2,C3,C8),dim=c(3,4)))
  cc1=Idl[8]+0.5*(Idl[7]+Idl[9])
  cc2=Idl[1]+Idl[2]+0.5*Idl[9]
  cc3=Idl[5]+0.5*Idl[6]
  cc8=Idl[4]+0.5*Idl[3]
  ASTAT=array(c(cc1,cc2,cc3,cc8),dim=4)
  
  imin1=0 #C1 C2 C3 C8
  dmin=1000
  for(i in 1:4)
  {
    d=norm_vec(At-STA[i,])
    if(d<dmin)
    {
      dmin=d
      imin1=i
    }
  }

  return (list(DistPC=DisPC,DistAP=DisAP,Ang=Ang, Idl=Idl,imin=imin,id=Set9[imin],imin1=imin1,dmin=dmin, ASTAT=ASTAT))
  # DistPC: Distance proj. of At on CFF plane to CFF cen
  # DistAP: Distance proj. of At on CFF plane to ref At 
  # Idl: ponderated Orientation code  c("C2","C2N4","N4","C8","C3","C3O2","O2","C1","O1","C2")
  # Ang: normalised angle (0-360) between C2->cen et cen->proj. of At on CFF plane
  # Imin: index of larger orientation code
  # Id:  CFF atom corresponding to Imin
  # Imin1: index (C1,C2,C3,C8) of CFF oxidizable atom closer from At
  # dmin: corresponding distance
  # ASTAT: ponderated projection of  idlindex on CFF C1 C2 C3 C8atoms
}

GetIminAng=function(C1,C2,C3,C8,M1,M2,M3,imin1)
{
  STA=t(array(c(C1,C2,C3,C8),dim=c(3,4)))
  V1=STA[imin1]-M3
  V2=M2-M3
  AngRec[8,]=cross(C1-C2,C1-C3)
  A=GetAngle(V1,V2,8,AngRec)$AA
  return (A)
}

CFFOD=function(RC,D4, D2,A2) # distance CCF to be oxidized atom to active O/FE
{
  A2=A2*pi/180
  val=sqrt(RC*(RC-2*D4cos(A2))+D4^2+D2^2)
  return (val)
}

FitAngle=function(ang)
{
  rm=c(0,93.12,143.06,264.56,360,453.12,503.06,624.56)
  rd=c(3.23,2.65,3.23,3.23)
  Angr=which(abs(rm-ang)==min(abs((rm-ang))))[1]
  AA=rm[Angr]%%360
  DCC=rd[(Angr-1)%%4+1]
  return (data.frame(Ang=AA,DCC=DCC))
}

SetLI=function(NP1,NP2,N4,N2)
{
  LI=cross(NP1,NP2)
  LI=LI/norm_vec(LI)
  wayLI=LI%*%(N4-N2)
  if(wayLI<0) LI=-LI
  return (LI)
}




GetParmFromAngle=function(M1,NP1, NP2,D2, D4,A2,A3,N2,N4, Way,RCFF,Rec=NULL)
{
  if(!is.null(Rec)) AngRec=Rec
  LI=SetLI(NP1,NP2,N4,N2)
  #F1=SetMP1(LI,NP1,NP2,M1,D2,1)
  F1=SetMP1(LI,NP1,NP2,M1,D2,1)
  MP1=F1$MP1
  M2=F1$M2
  AX1=(M2-MP1)/norm_vec(M2-MP1)
  Pt=M2+sign(D2)*D4*AX1#ajout de way*ax1
  m3=as.vector(rotate_in_space(Pt, NP2, M2, -A3))  
  Ax=(M2-m3)/norm_vec(M2-m3)
  Pt=m3+Ax*RCFF
  m4=as.vector(rotate_in_space(Pt, NP2, m3, -A2))
  AngRec[9,]=LI
  A4=GetAngle(M2-MP1,m3-MP1,9,AngRec)
  A6=0
  return (data.frame(MP1=MP1, M2=M2, M3=m3,M4=m4,A4=A4,M4=m4,A6=A6,MP1=MP1))
}

data=function()
{
  N1<<-c(42.624, 27.922, 30.456)
  N2<<-c(39.920, 27.478, 29.205)
  N3<<-c(39.515, 30.364, 29.015)
  N4<<-c(42.227, 30.728, 30.175)
  Fe<<-c(41.145, 29.115, 29.587)
  O1<<-c(41.846, 29.133, 28.099)
  C1<<-c(43.904, 28.109, 23.429)
  C3<<-c(46.302, 32.379, 26.112)
  C2<<-c(45.478, 26.668, 27.958)
  DHO=1.66 # distance oxygen -iron
}



compare=function(a,b)
{
  z=abs(a-b)/abs(a+b)*2<1e-5
  return (z)
}


library(bio3d)
library(pracma)
library(stringr)
library(plyr)
library(testit)
 
 rotateFPDB=function(PDBfile,frame,Residue, angle,file, mode,mirror=0,source=NULL) # mirror axis code: 0:none 1:c1,2:c2,3c3 to cen or gcen
  {
    RCF=Residue
    target=1
    angle=angle/180*pi
    if(is.null(source)) source="PDBfile.pdb"
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
  
  # rotation of substrate
id=1 #OMM id
frame=1 # data frame
Residue=671 #CFF producer
file="MirrorC1.pdb"
mode=1 # around center C1,C2,C3
angle=0 #rotation angle
file="X:/usr/local/g16/simulation/Built/OMM_rotamers/Cxx221_mcpbpyLL0.pdb"
for (i in seq(0, 360, by=30))
{
  ff=paste(i,"_","test.pdb", sep="")
  #rrr=rotatePDB(id,frame,Residue,i,ff, mode)
  rrr=rotateFPDB(id,frame,Residue,i,ff, mode)
}

rrr=rotateFPDB(id,frame,Residue,0,"mirrorC3.pdb", mode,3)# mirror en C2-cen axis
rrr=rotateFPDB(file,frame,671,30,"test.pdb",1,0)
rrr=rotateFPDB(file,frame,671,angle,"test.pdb", 1,mirror=0)

#tests
tt=readLines(file)
s=unlist(str_split(tt[]," "))
s=s[s!=""]
temp=tt[grep("CFF",tt)]
temp1=temp[grep(as.character(671),temp)]
r=gregexpr("[0-9\\.]+" , temp1[1] )
nr=r[[1]][r[[1]]>30]
#test end

# test
OMMId=1
source="OMM01/RefPDB.pdb"
Residue=671 #producer
angle=0
file="test0.pdb"
mode=1
mirror=0
RPDB0=rotateFPDB(OMMId,frame,Residue, angle,file, mode,mirror,source) # mirror axis code: 0:none 1:c1,2:c2,3c3 to cen or gcen
#end test

#adjust TER
rrr=rotateFPDB(file,671,330,"test.pdb",1,0)
tt=readLines("test.pdb")
TerPos=c(1,7954,10564,10636,10637,10638,10662,10713,10738)
u="REMARK, BUILD BY MCPB.PY"
for(i in 1:8)
{
  for(j in (TerPos[i]:(TerPos[1+i]-1))) u=rbind(u,tt[j])
  if(i!=8) u=rbind(u,"TER")
}
writeLines(u,"test330.pdb")


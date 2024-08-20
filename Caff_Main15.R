#setwd("C:/Users/dpomp/Desktop/R-Analysis3")
source("Caff_Fun15.r") 
require(tcltk)
#options(error = traceback)
#options(trace=T)


# Main processing routine (phase I)

#"MixLL/" c(1:12,62,63)
#"MixLN/" c(15,17:20,22,25:26,30:32,49:55,57:60,70,75)
#"MixRN/" c(16,24,33:41,43:48,61,71:74)   c(16,33:41,44:48,71:74,24,33:41,43,61)
#"MixRL/" c(13,14)
#"MixRR/" c(64:69)
#"Serie669/" c(100:123)
#"Serie671/" c(76:99) 
#"Serie669_1CFF/" c(136:147)
#"Serie671_1CFF/" c(124:135)
#"MixLN+" #c(15,17:20,22,25:26,30:32,49:55,57:60,70,75,136:147,124:135)
#--------------------------------------------------------------Main parameters settings
ForceTarget='OXY'# (NULL (default) 'FE' or 'OXY')  force reference atom on heme irrespective of active oxygen presence or absence in MDs (use with caution normal value: NULL)
ActivateNorm=1 # Generate time series of standard parameters (single run at a time only)/ 0: inaction,1: raw data,2: data norm 180,3 data norm 360
AnalyseCFFInteract<<-TRUE # When TRUE analyse interaction between two CFF when present (single CFF will cause errors)
AngleCorMode=1 # 0: no orientation reference, 1: reference is the first time of each OMMs,2 reference is the fist time of first OMM of batch
UseCountMode=TRUE # TRUE for NormalMode
# Use the first OMM of a batch  or (FALSE) a new PDB  reference for each OMM for angle alignments
# !!! Do not set to TRUE when two caffeins are present in the same OMM
SinglePDBRef=FALSE # see above 
ResDataSet="MixLN/"
OMMRange= c(15,17:20,22,25:26,30:32,49:55,57:60,70,75)
SetMainParm(ResDataSet,OMMRange, UseCountMode) # set parameter for data scanning
#End main parameter settings

#-----------------------------------------------------------phase I MD scan and coordinate conversion and statistics
RunPhaseI(setrange, TRUE) #!! run data scanning routines (erasing old ResDataSet)
#-----------------------------------------------------------End of phase I

#-----------------------------------------------------    Phase II CFF selection and data conversion (Use phase I results)
NumbCFF=1
SelectCFF=1 #1:effector or not defined 2:producer
SetPhaseIIParm(NumbCFF,SelectCFF) #Set parameters for phase II
res=RunPhaseIIA() # extract phase II data
rn=rownames(res[[1]])
cat(rn)
#---------------------------------------------------------End of phase II

#------------------------------------------------------Phase III statistical analysis using phase I and II results
pdf(paste(ResDataSet,"Ph3A.pdf",sep="")) # Path for results set of phase III

mode1=c(1,1,1,1,1,1,1,1)# selected displayed parameters reduced coordinate
mode2=c(1,1,1,1,1,1,1,1)# selected displayed parameters CFFs interaction
mode3=c(1,1,1,1,1)# selected displayed parameters CFF geometries
mode=c(mode1,mode2,mode3)
Multi=3 #Number of means columns in Dataset (does nor edit without reason)
# select only one of the two next line
re=list(id=res[[7]],xlab="Runs") #Normal mode
#re=list(id=as.character(seq(0,330,30)),xlabs="Rotation angle (°)") #alternate choice for rotamers

ShowCorrelation=TRUE
SortParm="dmin"#"imin"#"NULL" # #sort parameter for following phase IIIresults
# do not edit the next lines
LabScale=1 # 1 small vertical labels,2 larger horizontal labels 
SortParmList1=c("A0","A1","A2","A3","D2","D4","dmin","imin") 
SortParmList2=c("SCP1","SCP2","ACP1","ACP2","DCP1","DCP2","DCP1P","DCP2P")
SortParmList3=c("DC12","AHC1","AHC2","ORC1O","ORC2O")
SortParmList=c(SortParmList1,SortParmList2,SortParmList3)
SortParmNb=c(8,5,11,7,6,9,12,10,1,13,2,14,3,15,4,16,17,18,19,20,21)
#-------------

if(SortParm!="NULL")
{
  CorData=RunPhaseIII(res,mode,Multi,ShowCorrelation,re, SortParmNb[which(SortParmList==SortParm)],LabScale)
}else CorData=RunPhaseIII(res,mode,Multi,ShowCorrelation,re, NULL,LabScale)
write.table(CorData,file="test.txt")
dev.off()
#----------------------------------------------------------End of Phase III



#---------------------------------Phase IV (MDS analysis of data)

DataSet="Serie671/" # Considered Dataset in analysis (pahse IV-A)
DataSet1="Serie671" #Single or multiple Considered dataset (without'/') (phase IV-B)
#pdf(paste(DataSet,"/Ph5.pdf",sep="")) # Name for the PDF results (skip for screen)

Ofset=25 # Skipped frames at the start of each  OMM (-1: none) (skip initial state and to allow relaxation)
SelCFF=2 # selected CFF (1:effector, 2: producer)

#Parm=c("CHPA","OCPD", "ACHAC","ANG","DPCCO", "COC") #(A1), (D2), (A0), (D4), (dmin), (imin)
Parm=c('A0','A1','A3','A2','D2','D4','D5','A2b','IDENA','IDENB','IDENC','dmin','imin','imin1')
pp=c("A0","D4") # selected parameters for MDS analysis plots
ShowCor=TRUE # Show correlation plot
mask=c("A0","A1","A3","A2","D2","D4",'IDENA','IDENB','dmin','imin','imin1') #item to consider for correlation plot
Parm1=which(Parm==pp[1])
Parm2=which(Parm==pp[2])
#end setting

#phase IV-A
TB=TruncKin(DataSet,SelCFF,Ofset) # data extraction
NormDataPlot(DataSet,TB,Parm1,Parm2,mask, FALSE,ShowCor)
#end phase IV-A

mask1=c("A0","A1","A3","D2","D4") #mask1=c("A1","A3","D2") mask1=c("A0","D2","A3")
Filename="MDS_XXX.txt" #Name to save MDS reduced coordinates
IminRg=c(1:9) #Filter for imin values
DminRg=c(0,15) #Filter for dmin values
sel=1 # select imin for color codes

Id=c('A0','A1','A2','A3','D2','D4','dmin','imin') # considered descriptors name
Mask=c(1,2,4,5,6) #independent variables (descriptors indexes)
Obs=3             #dependent variable   (descriptor index)
Type=c(1,1,1,1,1,1,0,0) #Descriptor normalisation 0: none(raw), 1: normalized 2: sin/cos splitted
FilterTarget=NULL#"OXY" #("OXY") or "FE" or NULL for no filter
FN=FilterNorm(DataSet1,Ofset,SelCFF,Mask,Obs,Id,Type,FilterTarget) #Data filtering and normalisation

Filename=paste(DataSet1[1],'/MDS-',length(DataSet1),'.txt',sep="")
#next: long process: run only once after changing parameters
#Filename=paste(DataSet1[1],'/MDS.txt',sep="")
BuildMSDProj(FN$Env,Filename) #Very long process (recorded) use only 1 time before call of MDSProj

#Filename=paste(DataSet1[1],'/MDS.txt',sep="")
Filter="OXY" #("OXY") or "FE" or NULL for no filter
Legend=c(-5,3) #legend placement on figures
BarLegend="topright"
MDSProj(Filename,FN, IminRg,DminRg,Filter,Legend,BarLegend) #Projection of MDS reduced coordinates with imin color codes

# Plot filtered MDS map
IminRg=c(1:9) #selected orientations (NULL for All)
SelOMM=NULL #selected OMMs (NULL for All)
ThD=c(0,15) #selected dmin values (NULL for All)
ClustNum=1
SU=FilterMDS(DataSet,TB,Type,OfSet,ThD,SelOMM,IminRg,ClustNum)
dev.off()
# ---------------------------------------------------------end of phase V

#----------------------------------------- rotamer analysis (Phase V)
pdf(paste(ResDataSet,"Ph3B.pdf",sep="-"))
MyItem=11 # select parameter to plot
Sync=TRUE # match angle with initial values (use only for angles)
res1=SetPhaseIIBParm(MyItem,SetPrefix,SelectCFF) #set rainbow on plots parameters
res2=RunPhaseIIB(res,res1,MyItem,Sync) #prepare data
res=res2[[5]]
Mode=2 # mode of plot 1:rainbow, 2:violon, 3: special count base plot for CFF orientations

RunPhaseIIBPlot(Mode,MyItem,res,res1,res2)# Plot angular analysis
SavePhaseIIBPlotParm(MyItem,SelectCFF,res,res1,res2)# Save current plot parameters in database
dev.off()
#------------------------------------------End rotamer analysis

#-------------------------------------rotamer generation
# rotation of substrate
id=1 #OMM id
frame=1 # data frame
Residue=671 #CFF producer
file="test.pdb"
mode=1 # around center C1,C2,C3
angle=120 #rotation angle
file="X:/usr/local/g16/simulation/Built/OMM_rotamers/Cxx221_mcpbpyLL0.pdb"
for (i in seq(0, 360, by=30))
{
  ff=paste(i,"_","test.pdb", sep="")
  #rrr=rotatePDB(id,frame,Residue,i,ff, mode)
  rrr=rotateFPDB(id,frame,Residue,i,ff, mode)
}

rrr=rotatePDB(id,frame,Residue,0,"mirror_test.pdb", mode,4)# mirror en C2-cen axis
rrr=rotateFPDB(file,671,30,"test.pdb",1,0)
rrr=rotateFPDB(file,671,angle,"test.pdb", 1,mirror=0)

#tests
tt=readLines(file)
s=unlist(str_split(tt[]," "))
s=s[s!=""]
temp=tt[grep("CFF",tt)]
temp1=temp[grep(as.character(671),temp)]
r=gregexpr("[0-9\\.]+" , temp1[1] )
nr=r[[1]][r[[1]]>30]
#test end


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
#----------------------------------------end rotamer generation


#--------------------------------------------------------------------------Accessory code (to revalidate)
# plot specific correlation between cordata
a=4
b=3
plot(t(CorData)[,a], TCorData[,b], xlab=colnames(TCorData)[a], ylab=colnames(TCorData)[b], main=paste("p-value=",test$p[a,b]))
lines(lowess(TCorData[,a],TCorData[,b]), col="red")
#end

# rainbow parm edit frame
Limits="c(-500,500)" 
ScaleY="c(0,360)"
Ofset=0
AsAngle=FALSE
AsCOC=FALSE
AsStaking=FALSE
ylab=YlabV[2]
#end


par(mar=c(4,4,3,5)) 
plot(mat[10,],ty='b',axes=F,xlab="",ylab="", pch=16)
axis(side=4, at=c(1:8))
axis(side=1, at=c(1:12),labels=c("0","30","60","90","120","150","180","210","240","270","300","330"))
mtext("Orientation index",side=4,line=2.5)
par(new=T)
plot(mat[12,],axes=F,col="red", ty='b',ylim=c(0,10),xlab="",ylab="",lwd=3)
lines(mat[9,],col="blue",ty='b',lwd=2)
lines(mat[6,],col="green",ty='b',lwd=2, pch=15)
mtext("Distance ( Å )",side=2,line=2.5,pch=18)
axis(2, ylim=c(0,10))
mtext("Rotation angle (°)",side=1,col="black",line=2.5)


print(rownames(CorData))
a=mat[2,]
a[a>90]=180-a[a>90]
b=mat[14,]
b[b>90]=180-b[b>90]
plot(a,b)

plot(mat[1,],mat[12,])
rt=identify(mat[1,],mat[13,],re)
plot(mat[10,],mat[24,])


# start of correlation job
set1=1:dim(CorData)[2]
set2=Sets[[1]]
set3=Sets[[2]]
set4=CorData[12,]<5
ast=set1
M=cor(CorData[,ast])
colnames(M)=re[ast]
rownames(M)=re[ast]

#dendogramme
library(fastcluster)
library(ggplot2)
library(dendextend)
library("factoextra")
dist_mat <- dist(scale(M), method = 'euclidean')
arbre = hclust(dist_mat, method ='average')
plot(arbre)
inertie =sort(arbre$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
fviz_dend(arbre, k = 4,color_labels_by_k = TRUE,rect = TRUE)
NC=dim(M)[1]
classe = cutree(arbre, k = 4)
id=names(classe)
val=as.vector(classe)
ss=sort(val, index.return=TRUE)
S1=M[,ss$ix]
S2=S1[ss$ix,]
corrplot(S2)

#par(mfrow=c(2,2))
Mt=cor(t(CorData[,ast]))
PMA=colMeans(Mt)
rownames(Mt)=c("Stacking overlap1","CFF-PHE1 plane angle","CFF-PHE1 plane dist.","CFF-PHE1 Proj Dist","CFF-HEM plane angle","Fe/O CFFP dist.",
              "CFFP-Na_Nc angle","CFFP-Nb_Nd angle","CFF-Fe/O Proj. Dist.", "Orientation code","Orientation angle","Shorter Fe/O to CFF atom",
              "Stacking overlap2","CFF-PHE2 plane angle","CFF-PHE2 plane dist.","CFF-PHE2 Proj Dist")
colnames(Mt)=c("Stacking overlap","CFF-PHE plane angle","CFF-PHE plane dist","CFF-PHE Proj. Dist.","CFF-HEM plane angle","Fe/O CFFP dist.",
              "CFFP-Na_Nc angle","CFFP-Nb_Nd angle","CFF-Fe/O Proj. Dist.", "Orientation code","Orientation angle","Shorter Fe/O to CFF atom",
              "Stacking overlap2","CFF-PHE2 plane angle","CFF-PHE2 plane dist.","CFF-PHE2 Proj Dist")

Mt=Mt[!is.na(Mt[,1]),]
Mt=Mt[,!is.na(Mt[1,])]
corrplot(Mt)
# end of correlation job


plot(t(CorData[10,ast]),t(CorData[9,ast]))

par(mfrow=c(2,2))# bar mean +/- std deviation
PrintRec(5, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"Phe-CFF stacking overlap")
PrintRec(6, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"Phe-CFF plane diedral Angle")
PrintRec(7, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"Dist. CFFcen to Phe plane")
PrintRec(8, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"Shift dist. Phe and CFF cen")

par(mfrow=c(2,2))
PrintRec(9, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar)
PrintRec(10, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar)
PrintRec(11, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar)
PrintRec(12, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar)
par(mfrow=c(2,2))

PrintRec(13, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"CFF_cen to Fe/O proj. distance")
PrintRec(14, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"CFF orientation code")
PrintRec(17, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"Shorter distance Fe/O to CFF atoms")
PrintRec(16, PdfSet, WorkSet, asel$sel,  asel$sel1,Sort, Bar,"CFF orientation Angle")

par(mfrow=c(1,1))


# correlation dynamique (single frame)  
OMM="61"
CFF=2
WorkSetFile="Files/WorkSet.txt"
IndexFile="Files/Index.txt"

FullMode=TRUE # When true override the frame windows defined in "worset.txt
AnalyseCFFInteract=TRUE
DoPrint=TRUE #detailled statistic printing
PrintMode1=1 # NULL  No detaiiled Print mode any value details On
SetReadRefPDB=FALSE  # FALSE en standard
MaskType=NULL  # 0 if ReadRefPDB=TRUE or  =NULL en standard
MaskOverlay=NULL# Use NULL for setting on Worset.txt values
FullMode=TRUE # When true override  "worset.txt" and MaskOverlay setting
ws=read.table("Files/WorkSet.txt",header=TRUE,sep='\t')
idx=read.table(IndexFile,header=TRUE,sep='\t')
Set=which(ws[,3]==paste("OMM",OMM,sep="") & ws[,4]==CFF)

tt=file.exists(paste(ws$OMM[Set],"RefPDB.pdb",sep="/"))
if(SetReadRefPDB && tt) ReadRefPDB=TRUE else ReadRefPDB=FALSE
if(ReadRefPDB) MaskType=0
code=idx$Struct[ws$Id[Set]]
if(length(grep("c",code))==1 && length(grep("C",code))==1 && AnalyseCFFInteract) twoCFF=TRUE else twoCFF=FALSE

if(!is.null(MaskOverlay)) mmode=MaskOverlay
if(FullMode) mmode="1:9999" else mmode=ws$mode[u]

FN=paste(Set,"_",ws$OMM[Set],"CFF_",ws$CFF[Set],".pdf",sep="")
FN1=paste(Set,ws$OMM[Set],"_Id_",ws$Id[Set],"CFF_",ws$CFF[Set],".txt",sep="")
pdf(file=FN, width=7,height=5)

#appel de  fonction a mettre à jour !!
res=AnalyseCFF(ws$Id[Set],ws$Target[Set],ws$CFF[Set],mmode, FN,twoCFF,MaskType,PrintMode1, ReadRefPDB,DoPrint, FN1)
ss=PrintRes(res, paste(ws$OMM[Set],"_CFF",ws$CFF[Set],".txt", sep=""))

dev.off()
data1=res[[2]]
data2=res[[3]]
data3=res[[4]]
data4=res[[5]]
rdim=length(data1$AHC)
val=cbind(data1$AHC,data1$DFP,data1$AAC,data1$ADB,data1$DPHPC,data1$imin1,data1$Dist,
      data2$rec,data2$ang,data2$dist,data2$dist2,
      data3$rec,data3$ang,data3$dist,data3$dist2,
      data4$DCAB, data4$AHCA, data4$AHCB, data4$AHAB, data4$VCENA, data4$VCENB, data4$IDENA,data4$IDENB,data4$IDENC)
dmat=matrix(val, nrow=rdim,ncol=24)
M=cor(dmat)
rownames(M)=c("CFF-HEM plane angle","Fe/O CFFP dist.","CFFP-Na_Nc angle","CFFP-Nb_Nd angle","CFF-Fe/O Proj. Dist.","Orientation code",
              "Shorter Fe/O to CFF atom","Stacking overlap 1","CFF-PHE1 plane angle","Stacking distance1","Stacking ofset1",
              "Stacking overlap 2","CFF-PHE2 plane angle","Stacking distance2","Stacking ofset2",
              "Dist. CFF cenA-cenB","Ang. hem/CFFA planes","Ang. hem/CFFB planes","CFFA/CFFB plane ang.","Ang. C2-cenA/cenA-Proj.O","Ang. C2-cenB/cenB-Proj.O",
              "Closest CFFA atom to CFFB cen","Closest CFFB atom to CFFA cen","Closest CFF* atom to FE/O")
colnames(M)=rownames(M)
ss=c(16:24)
M=M[ss,ss]
D=dmat[,ss]
par(mfrow=c(1,1))
corrplot(M)

V1=D[,5]
V2=D[,4]
plot(V1,V2)
cor.test(V1,V2,use="complete.obs")
rg=lm(V2~V1)
abline(rg, col = 'red')

r=vector()
for(i in 1:10)
{
  for (j in (i+1):11)
  {
    r=c(r,M[i,j])
  }
}
barplot(r)

# rotation of substrate
id=1 #OMM id
frame=1 # data frame
Residue=671 #CFF producer
file="test.pdb"
mode=1 # around center C1,C2,C3
angle=120 #rotation angle
file="X:/usr/local/g16/simulation/Built/OMM_rotamers/Cxx221_mcpbpyLL0.pdb"
for (i in seq(0, 360, by=30))
{
  ff=paste(i,"_","test.pdb", sep="")
  #rrr=rotatePDB(id,frame,Residue,i,ff, mode)
  rrr=rotateFPDB(id,frame,Residue,i,ff, mode)
}

rrr=rotatePDB(id,frame,Residue,0,"mirror_test.pdb", mode,4)# mirror en C2-cen axis
rrr=rotateFPDB(file,671,30,"test.pdb",1,0)
rrr=rotateFPDB(file,671,angle,"test.pdb", 1,mirror=0)

#tests
tt=readLines(file)
s=unlist(str_split(tt[]," "))
s=s[s!=""]
temp=tt[grep("CFF",tt)]
temp1=temp[grep(as.character(671),temp)]
r=gregexpr("[0-9\\.]+" , temp1[1] )
nr=r[[1]][r[[1]]>30]
#test end


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



# Caff-1A2-Rscripts
R-Script linked to publication "Deep learning assisted molecular dynamics addresses control of P450 1A2 metabolic region-specificity"

Folder includes all R scripts used to build data of paper "Deep learning assisted molecular dynamics addresses control of P450 1A2 metabolic region-specificity" and requiered txt files used by scripts to describe data of molecular dynamics.

Molecular dynamic files by themselves are not included due to size, can be accessed separatly and are structured as follow for each individual dynamic:
- Folder name : OMMxxx where xxx is the index number of the dynamic. These index numbers are used to build dataset used in R scripts as defined for example in lines 10 to 19 of CAFF_Main15.r script. Dataset anme and correspopnding  list of molecular dynamic runs have to be set in lines 29 and 30 of the Main R-script.
- Folder content: Each folder must contain at least tree files: 1/ a "data.pdb" file giving the PDB structure associated with a dynamic (can be any frame of it).2/ a "data.dcd" file : a dcd file (OpenMM) describing the dynamic. 3/ a "readme.txt file" describing the dynamic. Complementary informations can be retrieved  on the optionally included "OMM_solv.prmtop" files, on extracted frames (as (.pdb) and on the Github includes indextxt, TimeSeries.txt, Workset.txt, plotParm.txt files using as key the folder index number.

Runing scripts generate many intermediate and results files that must be manually transfered into a new folder named using the previously defined dataset name. These intermediate files can be used by scripts to rerun parts of scripts within requiring runing again all scripts. 

Script organisation:
Many parameters present in the different parts of scripts mustbe adjusted before run as follow:
Main Script
-  line 21: Used either "NULL" (no filter), "FE" (only P450 on iron III or II states) or "OXY" (only P450 on OXY ferrous state) are considered in datasets.
- line 22 : defaut value 1 must be used excepted for special features.
- line 23 : TRUE for 2 CFF at active site FALSE in the contrary
- line 24 : defaut value 1 must be used excepted for special features.
- line 25 : defaut value TRUE must be used excepted for special features.
- line 28 : default value FALSE (see line 26 and 27 comments)
- line 29 and 30 : to be set to define batch of data (see previously)
- line 39 and 40 : number and id (1 or 2) of CFF at active site to consider.
- line 50-53 : keep default values (or read scripts to change display details)
- line 59 : default TRUE
- line 60: Select the default parameter for sorting results of phase III calculation
- line 82, 83 : adjust for dataset according to comment.
- line 86 : number of skipped frames at the begining of each time series (allow for equilibation). default 25.
- line 87 : see comment
- line 90 : menbers of list can be used to set line 60
- line 91-93 : set according to comments
- line 103 : define mask for phase IVA (see examples on comments)
- line 105-107 : set ranges of considered Imin and Dmin values
- lines 109-113,122-124;128-131,138-139,143 : set display filters according to comments.
WithKeras15 script: this part supposes that datasets have been previously defined and processed using phase I-IV of Main scripts and that suitable computing environment (including python)  was previously set to run keras and tensorflow libraries. 
- line 42 : list of datasets used in calculation

# Caff-1A2-Rscripts
R-Script linked to publication "Deep learning assisted molecular dynamics addresses control of P450 1A2 metabolic region-specificity"

Folder includes all R scripts used to build data of paper "Deep learning assisted molecular dynamics addresses control of P450 1A2 metabolic region-specificity" and requiered txt files used by scripts to describe data of molecular dynamics.

Molecular dynamic files by themselves are not included due to size, can be accessed separatly and are structured as follow for each individual dynamic:
- Folder name : OMMxxx where xxx is the index number of the dynamic. These index numbers are used to build dataset used in R scripts as defined for example in lines 10 to 19 of CAFF_Main15.r script. Dataset anme and correspopnding  list of molecular dynamic runs have to be set in lines 29 and 30 of the Main R-script.
- Folder content: Each folder must contain at least tree files: 1/ a "data.pdb" file giving the PDB structure associated with a dynamic (can be any frame of it).2/ a "data.dcd" file : a dcd file (OpenMM) describing the dynamic. 3/ a "readme.txt file" describing the dynamic. Complementary informations can be retrieved  on the optionally included "OMM_solv.prmtop" files, on extracted frames (as (.pdb) and on the Github includes indextxt, TimeSeries.txt, Workset.txt, plotParm.txt files using as key the folder index number.

 

Script organisation:
Caff_Main15.R and withkeres15.R are the two entry points, the second being dependant on results from the first. Other R files are utilities files. Provided .txt files provides essential informations to process primary data that must be separatly loaded as a set of folders and placed in the same folder that the R-Scripts. The provided .txt files must be placed in a "Files" folder also placed in the same folder than scripts.
Main Script: this part is divised on several "phases" commented on the script that can be run independantly provided to be run at least one time in order to build intermediate result files that are reuused by the next phases. Runing scripts with a set of primary data generate many intermediate and result files that must be manually transfered into a new folder named using the script defined dataset name. These intermediate files can be used by scripts to rerun parts of scripts without requiring runing again all scripts. However, all scripts for a given phase must be run at a whole. Phase I that process primary data to intermediate files is the longer to run (could take hours) but need only to be run once for each defined dataset (list of molecular dynamic runs). Take care to conserve parameters consistancy between phases.

Many parameters present in the different parts of scripts mustbe adjusted before run as follow:
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
WithKeras15 script: this part which involmves deep learning processes supposes that datasets have been previously defined and processed using phase I-IV of Main scripts and that suitable computing environment (including python)  was previously set to run keras and tensorflow libraries.
This script can be run / rerun by parts: it is recommended to read it and understand roles of parts before using. Executing at a whole migth mask intermediate results and graphs.
- line 42 : list of datasets used in calculation
- line 43 : corresponding names used in graphs.
- line 44 : selected CFF for each dataset
- line 45 : see comment and main script line 86
- line 46-48 : select parameter index and parameter type for learning and prediction based on ids of line 49
- line 50 : learning and testing fractions of data
- line 51, 52 : filter on a specific P450 iron redox state
- Line 54-56 : perform or not randomisation of the order of input or output data
- Lines 57-58 : data ponderation and color sets
- lines 60-62: filters  (learn, test, ignore) of batches defined according to comments
- lines 64-65, 128-129 : set according to comments


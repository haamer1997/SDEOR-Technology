# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 09:02:39 2020

@author: ahass16
"""
#%% Import packages and libraries
import numpy as np
import os
#%% Generate MATLAB Files and HPC scripts
for filename in os.listdir("../"): #loop through EOR MATLAB scripts in parent folder
    if filename.endswith(".m"): #grab only MATLAB file scripts
        file = open("../"+filename, 'r') #open EDFM MATLAB Files
        Lines = file.readlines() #read EDFM MATLAB File
        NewMATLABName =  filename[:filename.find('.')] + '_pEDFM' +".m" #
        Lines[121] ='     fracplanes(i).perm=1*nano*darcy; \n'
        Lines[137] = 'G = pMatFracNNCs3D(G,tol); \n'
        Lines[139] = '%TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol); \n'
        Lines[141] = 'TPFAoperators = setupPEDFMOpsTPFA(G, G.rock, tol); \n'
        s = Lines[285] 
        Lines [285] = s[:s.find('.')] + '_pEDFM' +".mat']; \n"
        f = open(NewMATLABName, "x")
        for line in Lines: 
            f.write(line)
        f.close()

        NewPBSName = NewMATLABName[:NewMATLABName.find('.')] +".pbs"
        Lines = ['#!/bin/bash \n','#PBS -q workq \n','#PBS -A loni_unconvrs20 \n', '#PBS -l nodes=1:ppn=2 \n', '#PBS -l walltime=10:30:00 \n']
        Lines.append('#PBS -o ' + NewMATLABName[NewMATLABName.find('_')+1:NewMATLABName.find('.')] + ' \n')
        Lines.append('#PBS -j oe \n')
        Lines.append('#PBS -N ' + NewMATLABName[NewMATLABName.find('_')+1:NewMATLABName.find('.')] +' \n')
        Lines.append('#PBS -m e \n')
        Lines.append('#PBS -M ahass16@lsu.edu \n')
        Lines.append('\n')
        Lines.append('cd /scratch/ahass16/mrst-2020a \n')
        Lines.append('module load matlab/r2019b \n')
        Lines.append('matlab -nodisplay -nosplash -r startup \n')
        Lines.append('matlab -nodisplay -nosplash -r "addpath(genpath(' + "'/scratch/ahass16/mrst-2020a/modules/shale_chapter')); " + NewPBSName[:NewPBSName.find('.')] + ';" \n')
        f = open(NewPBSName, "x")
        for line in Lines: 
            f.write(line)
        f.close()  
#%% Generate PBS Files for HPC
#for filename in os.listdir("../"):
#    if filename.endswith(".pbs"): 
#        file = open("../"+filename, 'r') 
#        Lines = file.readlines()
#        NewPBSName = filename[:filename.find('.')] + '_pEDFM' +".pbs"
#        Lines[5] = '#PBS -o 8NF_73_pEDFM \n'
#        Lines[7] = '#PBS -N 8NF_73_pEDFM \n'
#        Lines[14] = 'matlab -nodisplay -nosplash -r "addpath(genpath(' + "'/scratch/ahass16/mrst-2020a/modules/shale_chapter')); " + NewMATLABName[:NewMATLABName.find('.')] + ';" \n'
#        f = open(NewPBSName, "x")
#        for line in Lines: 
#            f.write(line)
#        f.close()    
#%% Generate Bash Script for Easy run
directory = r'C:\Users\ahass16\Desktop\Research Work\mrst-2020a\modules\shale_chapter\eor_examples\Percolation Caclutions\Percolation Cases\EDFM\Refined mesh and CPR-AMGCL\HPC Cases\Uncertainty Analysis\pEDFM'
f = open("Uncertainty_pEDFM.sh", "x")
f.write("#!/bin/bash \n")
for filename in os.listdir(directory):
    if filename.endswith(".pbs"):
        f.write("dos2unix " + filename + "; \n")
        f.write("qsub " +  filename +"; \n")
f.close()        
#%%Test one MATLAB file first
#filename = 'Percolation_8NF_EOR_73.m'
#file = open(filename, 'r') 
#Lines = file.readlines()
#NewMATLABName =  filename[:filename.find('.')] + '_pEDFM' +".m"
#Lines[121] ='     fracplanes(i).perm=1*nano*darcy; \n'
#Lines[137] = 'G = pMatFracNNCs3D(G,tol); \n'
#Lines[139] = '%TPFAoperators = setupShaleEDFMOpsTPFA(G, G.rock, tol); \n'
#Lines[141] = 'TPFAoperators = setupPEDFMOpsTPFA(G, G.rock, tol); \n'
#s = Lines[285] 
#Lines [285] = s[:s.find('.')] + '_pEDFM' +".mat']; \n"
#f = open(NewName, "x")
#for line in Lines: 
#    f.write(line)
#f.close()
##%%TEST ONE HPC Script
#filename = "8NF_73.pbs"
#file = open(filename, 'r') 
#Lines = file.readlines()
#NewPBSName = filename[:filename.find('.')] + '_pEDFM' +".pbs"
#Lines[5] = '#PBS -o 8NF_73_pEDFM \n'
#Lines[7] = '#PBS -N 8NF_73_pEDFM \n'
#Lines[14] = 'matlab -nodisplay -nosplash -r "addpath(genpath(' + "'/scratch/ahass16/mrst-2020a/modules/shale_chapter')); " + NewMATLABName[:NewMATLABName.find('.')] + ';" \n'
#f = open(NewPBSName, "x")
#for line in Lines: 
#    f.write(line)
#f.close()
##%% Generate one bash script for HPC
#f = open("UncertaintypEDFM.sh", "x")
#f.write("#!/bin/bash \n")
#f.write("dos2unix " + NewPBSName + ";\n")
#f.write("qsub " +  NewPBSName +";\n")
#f.close()
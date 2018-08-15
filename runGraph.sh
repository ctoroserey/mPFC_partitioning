#!/bin/bash -l

# load modules within SCC
module load R/3.4.3
module list

SubjID=$1

# Setups
#$ -l mem_total=48G
#$ -pe omp 16
#$ -j y

# Prep subject
sh prepSubj.sh $SubjID

# Run the R script that calculates community detection for a set of participants
Rscript graphAnalysis_solo.R $SubjID

# Clean up
rm -f ${SubjID}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii
rm -f ${SubjID}_rfMRI_REST1_LR.dtseries.nii
rm -f ${SubjID}_rfMRI_REST1_RL.dtseries.nii
rm -f ${SubjID}_rfMRI_REST2_LR.dtseries.nii
rm -f ${SubjID}_rfMRI_REST2_RL.dtseries.nii
rm -f ${SubjID}.ptseries.nii
rm -f ${SubjID}.corrThickness.32k_fs_LR.dscalar.nii
rm -f ${SubjID}.curvature.32k_fs_LR.dscalar.nii
rm -f ${SubjID}.MyelinMap_BC.32k_fs_LR.dscalar.nii

mv ${SubjID}_myelin.txt ./MyelinMaps
mv ${SubjID}_curvature.txt ./Curvature
mv ${SubjID}_thickness.txt ./Thickness

#!/bin/bash -l

# this script will copy a subject's dense time series, turn it into a csv, and run community detection on it
# its output is just the eigenvectors and eigenvalues of SP, plus the determined modularity membership

#$ -l mem_total=94G
#$ -N commDetection
module load connectomewb/1.2.3

subj=$1

# copy data
cp /projectnb/connectomedb/Q6_postprocessing/${subj}/MNINonLinear/Results/concatenated_rfMRI/rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii \
    /restricted/projectnb/cd-lab/Claudio/Community/${subj}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii

# transform to csv
wb_command -cifti-convert -to-text ${subj}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii \
                                    temp_${subj}_timeSeries.csv

# get only the cortical vertices
awk 'NR <= 59412 { print }' temp_${subj}_timeSeries.csv >> ${subj}_timeSeries.csv

# remove original and temp time series, and move time series to the appropriate folder
rm -f ${subj}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii
rm -f temp_${subj}_timeSeries.csv
mv -f ${subj}_timeSeries.csv tseries

# run the R script that will produce the final counts
Rscript community_detection.R $subj



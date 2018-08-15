#!bin/bash/

subj=$1

echo "Preparing data for subject ${subj}"

# Copy data
cp /projectnb/connectomedb/Q6/${subj}/MNINonLinear/fsaverage_LR32k/${subj}.corrThickness.32k_fs_LR.dscalar.nii /restricted/projectnb/cd-lab/Claudio/Community/
cp /projectnb/connectomedb/Q6/${subj}/MNINonLinear/fsaverage_LR32k/${subj}.curvature.32k_fs_LR.dscalar.nii /restricted/projectnb/cd-lab/Claudio/Community/
cp /projectnb/connectomedb/Q6/${subj}/MNINonLinear/fsaverage_LR32k/${subj}.MyelinMap_BC.32k_fs_LR.dscalar.nii /restricted/projectnb/cd-lab/Claudio/Community/

cp /projectnb/connectomedb/Q6_postprocessing/${subj}/MNINonLinear/Results/concatenated_rfMRI/rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii /restricted/projectnb/cd-lab/Claudio/Community/${subj}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii
cp /projectnb/connectomedb/Q6_postprocessing/${subj}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.int.bpf.residMGSR.censor.demean.dtseries.nii /restricted/projectnb/cd-lab/Claudio/Community/${subj}_rfMRI_REST1_LR.dtseries.nii
cp /projectnb/connectomedb/Q6_postprocessing/${subj}/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_hp2000_clean.int.bpf.residMGSR.censor.demean.dtseries.nii /restricted/projectnb/cd-lab/Claudio/Community/${subj}_rfMRI_REST1_RL.dtseries.nii
cp /projectnb/connectomedb/Q6_postprocessing/${subj}/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_hp2000_clean.int.bpf.residMGSR.censor.demean.dtseries.nii /restricted/projectnb/cd-lab/Claudio/Community/${subj}_rfMRI_REST2_LR.dtseries.nii
cp /projectnb/connectomedb/Q6_postprocessing/${subj}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_Atlas_hp2000_clean.int.bpf.residMGSR.censor.demean.dtseries.nii /restricted/projectnb/cd-lab/Claudio/Community/${subj}_rfMRI_REST2_RL.dtseries.nii

# Convert structural data
wb_command -cifti-convert -to-text ${subj}.MyelinMap_BC.32k_fs_LR.dscalar.nii ${subj}_myelin.txt
wb_command -cifti-convert -to-text ${subj}.curvature.32k_fs_LR.dscalar.nii ${subj}_curvature.txt
wb_command -cifti-convert -to-text ${subj}.corrThickness.32k_fs_LR.dscalar.nii ${subj}_thickness.txt

echo "Parcellating time series..."
wb_command -cifti-parcellate ${subj}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii  \
                              Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii \
                              COLUMN ${subj}.ptseries.nii -method MEAN

# Get each parcel's mean timeseries
wb_command -cifti-convert -to-text ${subj}.ptseries.nii ${subj}_ptSeries.txt

# Get the overall time series
echo "Storing the whole time series..."
wb_command -cifti-convert -to-text ${subj}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii \
                                    temp_${subj}_timeSeries.csv
awk 'NR <= 59412 { print }' temp_${subj}_timeSeries.csv >> ${subj}_timeSeries.csv

rm -f temp_${subj}_timeSeries.csv

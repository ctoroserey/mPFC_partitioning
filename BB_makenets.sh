#!/bin/bash -l

# grab a seed from a dconn.nii file, produce its dense scalar form, then transform to txt
# using cifti math produces unreadable dscalars (it assumes they are dconns), which is why I add txt transforms

subj=$1
rowA=$2 # row index for net A
rowB=$3 # row index for net B

# get seed for both nets
wb_command -cifti-math 'x' ${subj}_A.dscalar.nii -var x ${subj}_corr.dconn.nii -select 1 $rowA
wb_command -cifti-math 'x' ${subj}_B.dscalar.nii -var x ${subj}_corr.dconn.nii -select 1 $rowB

# combine
wb_command -cifti-math 'B - A' ${subj}_temp.dscalar.nii -var A ${subj}_A.dscalar.nii -var B ${subj}_B.dscalar.nii 

# to text
wb_command -cifti-convert -to-text ${subj}_temp.dscalar.nii ${subj}_combined.txt

# from text
wb_command -cifti-convert -from-text ${subj}_combined.txt sample.dscalar.nii ${subj}_combined.dscalar.nii

rm -f ${subj}_temp.dscalar.nii
rm -f ${subj}_combined.txt
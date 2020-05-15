#!bin/bash/

# this will take ALL the outputs from HCPOut (from communityDetection.R) and turn them into CIFTIs
# these can then be viewed on wb_command with a 60k surface model from the HCP

for dataFile in *dataforCifti.txt ; do

  fname=$(echo $dataFile | cut -f 1,2,3,4 -d '_')

  wb_command -cifti-convert -from-text ${dataFile} ./Basic_files/template.dscalar.nii ${fname}.dscalar.nii

done

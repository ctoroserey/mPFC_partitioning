#!bin/bash/

for dataFile in *dataforCifti.txt ; do

  fname=$(echo $dataFile | cut -f 1,2 -d '_')

  wb_command -cifti-convert -from-text ${dataFile} 100307.MyelinMap_BC.32k_fs_LR.dscalar.nii ${fname}.dscalar.nii

done

#!bin/bash/

# This will produce text files for each subject's CIFTI within a directory
# First input to the function should be the general filename to follow (e.g. .MyelinMap_BC.32k_fs_LR.dscalar.nii)
# Second should be the suffix of the generated file (e.g. myelin)
dataName=$1
dataType=$2

# echo *${dataName}
# echo $dataType

# Turn each file into a text file
for dataFile in *${dataName} ; do

  fname=$(echo $dataFile | cut -f 1 -d '.')

  echo $fname >> ${dataType}_subjList.txt

  wb_command -cifti-convert -to-text $dataFile ${fname}_${dataType}.txt

done

# And combine
# Note: this creates a tab-separated file.

# Increase the upper limit of open files if necessary (for paste)
nFiles=$(ls -l | wc -l)
Limit=2
if [ "$nFiles" -gt 500 ] ; then
  ulimit -Hn 500 # The hard limit
  ulimit -Sn 500 # The soft limit
fi

paste *_${dataType}.txt >> ${dataType}All.txt

#!/bin/bash -l

# this function will take a set of labels (from Glasser, .csv) and create an atlas that contains only those ROIs
# no statistics are computed, it's just a visualization mask

dirName=BBC

#function createAtlas () {
echo "#!ascii label, from condition ${dirName}" >> ${dirName}_lh_atlas.label
echo "#!ascii label, from ondition ${dirName}" >> ${dirName}_rh_atlas.label
totalLeft=0
totalRight=0
while IFS=$'\t' read -r -a label ; do
  	if [[ $label == L* ]] ; then
			totalLines=$(expr 2 + $(awk 'NR == 2 {print $1}' $HOME/fMRI/Glasser/label/lh/lh.${label}.label))
    			totalLeft=$(expr $totalLines + $totalLeft)
    			sed -n 3,${totalLines}p $HOME/fMRI/Glasser/label/lh/lh.${label}.label >> ${dirName}_lh_atlas.label
    	else
    			totalLines=$(expr 2 + $(awk 'NR == 2 {print $1}' $HOME/fMRI/Glasser/label/rh/rh.${label}.label))
    			totalRight=$(expr $totalLines + $totalRight)
    			sed -n 3,${totalLines}p $HOME/fMRI/Glasser/label/rh/rh.${label}.label >> ${dirName}_rh_atlas.label
    	fi
done < ${dirName}_threshLabels.csv
sed -i "2i "$totalLeft"" ${dirName}_lh_atlas.label
sed -i "2i "$totalRight"" ${dirName}_rh_atlas.label
#}

view $dirName

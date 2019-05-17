#!/bin/bash -l

# This script will take volumetric coordinates from either DN or SV literatures and project them to fsaverage
# USE: `sh metacoordsVol2Surf.sh metacoordsDN` (or metacoordsPOS for subjective value)
# IMPORTANT: THIS PROCESS TAKES TIME, ESPECIALLY WHEN PRODUCING SPHERES.

echo "TRANSFORMING METANALYSIS COORDINATES INTO CORTICAL PARCELS (Glasser et al., 2016)"

# file name will become a directory
# read filename
filename=$1 
echo $filename


# create dirs, including those to make FreeSurfer commands compatible
# important for step 7 and annot2label transformation
mkdir -p $PWD/$filename/results
mkdir -p $PWD/Glasser/surf
mkdir -p $PWD/Glasser/label/lh
mkdir -p $PWD/Glasser/label/rh

if [ ! -f ./Glasser/label/lh.HCP_MMP1p0.annot ] ; then
    cp ./Basic_files/*.HCP_MMP1p0.annot ./Glasser/label # make sure to do it for both hemis, I'm using Sean's HCP annotation file
    cp $SUBJECTS_DIR/fsaverage/surf/*h.white ./Glasser/surf/
fi
if [ "$(ls -A ./Glasser/label/lh)" ] ; then
    echo
    echo "HCP labels exist..."
    echo
else
    echo
    echo "Creating HCP label files from HCP annot..."
    echo
     # --border <fname>.gii creates a single overlay of binary borders
    mri_annotation2label --subject Glasser --hemi lh --sd $PWD --annotation HCP_MMP1p0 --outdir ./Glasser/label/lh
    mri_annotation2label --subject Glasser --hemi rh --sd $PWD --annotation HCP_MMP1p0 --outdir ./Glasser/label/rh
fi

# copy MNI template to new dir and cd to it to create spheres
if [ ! -f ./$filename/MNI152_T1_1mm_brain.nii.gz  ] ; then
    cp ./Basic_files/MNI152_T1_1mm_brain.nii.gz ./$filename
    cd ./$filename
else
    cd ./$filename
fi

# parameters
sphereSze=10 ; # size of the sphere to be projected to the cortical surface
labelThresh=1; # Minimum n of vertices needed to map single coordinate spheres to a Glasser parcel (previously 10, but decided to avoid this threshold)
printf "%s," Study Vertices_count | awk '{print $0}' >> noMatch.csv

# index so that the MNI atlas is only multiplied by 0 once
startIndex=0

### loop reading each study and converting it ###
for coordFile in coords*.txt ; do
    fname=$(echo $coordFile | cut -f 1 -d '.')
    echo
    echo "Working on ${fname}..."
    echo

    ## 1. generate an MNI template with the coordinate points
    if [ ! -f ${fname}_Points.nii.gz ] ; then
		# get the first set of coordinates, each num = (x,y,z)
		while IFS=$'\t' read -r -a coords ; do
			num1=$(printf "%.0f" "${coords[0]}")
			num2=$(printf "%.0f" "${coords[1]}")
			num3=$(printf "%.0f" "${coords[2]}")
			# sanity check
			echo $num1
      	    echo $num2
      	    echo $num3
            # create points using AFNI
			# MNI template is first multiplied by 0 to ensure that only the relevant coordinates are present 
      	    if [ $startIndex = 0 ] ; then
            	echo
            	echo "Setting coordinate foci..."
       			echo
            	3dcalc -LPI -a MNI152_T1_1mm_brain.nii.gz -prefix points.nii.gz -expr "a*0 + (equals(x,$num1)*equals(y,$num2)*equals(z,$num3))"
            else
            	3dcalc -LPI -a ${fname}_Points.nii.gz -prefix points.nii.gz -expr "a + (equals(x,$num1)*equals(y,$num2)*equals(z,$num3))"
      	    fi
      	    mv points.nii.gz ${fname}_Points.nii.gz
      	    startIndex=1
      	done <${coordFile}
    else
      	echo
      	echo "Points file for ${fname} already exists..."
      	echo
    fi

    ## 2. convert points to spheres using FSLMATHS
    if [ ! -f ${fname}_spheres.nii.gz ] ; then
      	echo
      	echo "Converting foci to spheres..."
      	echo
      	fslmaths ${fname}_Points.nii.gz -kernel sphere $sphereSze -dilF -bin ${fname}_spheres.nii.gz
    else
      	echo
      	echo "Sphere file already exists..."
      	echo
    fi

    ## 3. volumetric to cortical surface conversion using white matter folding
    ## surface templates in $FREESURFER_HOME/subjects/fsaverage/surf/
    ## note: appending "--nvox ${fname}_nvox.csv --srchit ${fname}_vox.nii.gz" gives num of surviving converted voxels (not specific to mask though)
    if [ ! -f ${fname}_lh.gii ] ; then
      	echo
      	echo "Converting volumetric mask into cortical vertices..."
      	echo
      	mri_vol2surf --src ${fname}_spheres.nii.gz --o ${fname}_lh.gii --srcreg $FREESURFER_HOME/average/mni152.register.dat --hemi lh --surf white --projfrac-max 0 1 0.25
      	mri_vol2surf --src ${fname}_spheres.nii.gz --o ${fname}_rh.gii --srcreg $FREESURFER_HOME/average/mni152.register.dat --hemi rh --surf white --projfrac-max 0 1 0.25
    else
      	echo
      	echo "Surface files already exist..."
      	echo
    fi

    ## 4. extracting boolean vector (is vertex in sphere, yes/no?) for all vertices to match with Glasser parcels
    ## wb_command properly transforms data into ASCII, so don't use FS, AFNI, or FSL for this
	## personal note: potentially no need to use xmllint if wb_command -gifti-convert is used instead
    if [ ! -f ${fname}_lh_binaries.csv ] ; then
      	echo
      	echo "Converting GIFTI into ASCII for binary extraction..."
      	echo
      	wb_command -gifti-convert ASCII ${fname}_lh.gii ${fname}_lh_binaries
      	wb_command -gifti-convert ASCII ${fname}_rh.gii ${fname}_rh_binaries
      	echo "cat //GIFTI/DataArray/Data" | xmllint --shell ${fname}_lh_binaries | sed '/^\/ >/d' | sed 's/<[^>]*.//g' >> ${fname}_lh_binaries.csv
      	echo "cat //GIFTI/DataArray/Data" | xmllint --shell ${fname}_rh_binaries | sed '/^\/ >/d' | sed 's/<[^>]*.//g' >> ${fname}_rh_binaries.csv
		sed -i 1d ${fname}_lh_binaries.csv
      	sed -i 1d ${fname}_rh_binaries.csv
    else
      	echo
      	echo "Binary files already exist..."
      	echo
    fi

    ## 5. match boolean vector to labels and create parcel list
    if [ ! -f ${fname}_matchLabels.csv ] ; then
      	echo
      	echo "Matching mask vertices to labels..."
      	echo
      	paste ${fname}_lh_binaries.csv ./Glasser_labels_lh.csv > ${fname}_lh_allLabels.csv
      	paste ${fname}_rh_binaries.csv ./Glasser_labels_rh.csv > ${fname}_rh_allLabels.csv
      	(awk '$1 == 1' ${fname}_lh_allLabels.csv ; awk '$1 == 1' ${fname}_rh_allLabels.csv) > ${fname}_matchLabels.csv
      	sort -u -o ${fname}_u_matchLabels.csv ${fname}_matchLabels.csv && sed -i 1d ${fname}_u_matchLabels.csv
		empty=$(awk -F "\"*,\"*" '{print $1}' ${fname}_matchLabels.csv | cut -c10 | grep -c -e '^$')
		printf '%s,' ${fname} $((empty / 2)) >> noMatch.csv # store the per-study match-less vertex #
      	#rm ${fname}_lh_binaries.csv ${fname}_rh_binaries.csv ${fname}_rh_allLabels.csv ${fname}_lh_allLabels.csv
    else
      	echo
      	echo "Matched list already exists..."
      	echo
    fi

    ## 6. threshold labels to include (visually checked with Excel query of ${fname}_matchLabels.csv)
	## this step counts the number of vertices in each parcel, and compares this to the prespecified threshold (see labelThresh in 'parameters')
    if [ ! -f ${fname}_threshParcels.csv ] ; then
      	echo
      	echo "Picking out thresholded parcels..."
      	echo
      	printf "%s," Label Vertices_count | awk '{print $0}' >> ${fname}_summary.csv
      	while IFS=$'\t' read -r -a label ; do
      	    count="$(grep -c ${label[1]} ${fname}_matchLabels.csv)"
      	    if [ "$count" -ge "$labelThresh" ]  ; then
				echo ${label[1]} >> ${fname}_threshParcels.csv
				printf "%s," ${label[1]} $count | awk '{print $0}' >> ${fname}_summary.csv
      	    fi
      	done < ${fname}_u_matchLabels.csv
    else
      	echo
      	echo "Thresholded parcel file already exists..."
      	echo
    fi

    ## 7. merge labels to create a study-specific label atlas to be visualized in FreeView
    if [ ! -f ./results/${fname}_lh_atlas.label ] ; then
      	echo
      	echo "Creating study atlas..."
      	echo
      	echo "#!ascii label, from subject ${fname}" >> ${fname}_lh_atlas.label # header
      	echo "#!ascii label, from subject ${fname}" >> ${fname}_rh_atlas.label
      	totalLeft=0
      	totalRight=0
      	while IFS=$'\t' read -r -a label ; do
      	    if [[ $label == L* ]] ; then
				totalLines=$(expr 2 + $(awk 'NR == 2 {print $1}' $HOME/fMRI/Glasser/label/lh/lh.${label}.label))
				totalLeft=$(expr $totalLines + $totalLeft)
				sed -n 3,${totalLines}p $HOME/fMRI/Glasser/label/lh/lh.${label}.label >> ${fname}_lh_atlas.label
      	    else
				totalLines=$(expr 2 + $(awk 'NR == 2 {print $1}' $HOME/fMRI/Glasser/label/rh/rh.${label}.label))
				totalRight=$(expr $totalLines + $totalRight)
				sed -n 3,${totalLines}p $HOME/fMRI/Glasser/label/rh/rh.${label}.label >> ${fname}_rh_atlas.label
      	    fi
      	done < ${fname}_threshParcels.csv
      	sed -i "2i "$totalLeft"" ${fname}_lh_atlas.label
      	sed -i "2i "$totalRight"" ${fname}_rh_atlas.label
        # move final summary files to the results folder
        mv -f -n ${fname}_summary.csv ./results/
        mv -f -n ${fname}_*h_atlas.label ./results/
    else
      	echo
      	echo "Atlas files already exist..."
      	echo
    fi

    echo
    echo "Done with ${fname}..."
    echo

    startIndex=0
done

###### Summarize condition #########

cd ./results

dirName=$(basename $(dirname $(pwd)))
echo
echo "Summarizing ${dirName}..."
echo

numStudies=0

# concatenate the surviving parcels across studies
for file in *summary.csv ; do
	awk -F "\"*,\"*" '{print $1}' $file | sed 1d >> ${dirName}_studyLabels.csv
	numStudies=$(expr 1 + $numStudies)
done

# get the proportion of times each parcel was reported across studies
while IFS=$',' read -r -a label ; do
    	count="$(grep -c $label ${dirName}_studyLabels.csv)"
   	if [ "$count" -ge 1 ] ; then
		fraction=$(bc <<< "scale = 2; $count / $numStudies")
		echo $label >> ${dirName}_threshLabels.csv
   	 	printf "%s," $label $fraction | awk '{print $0}' >> ${dirName}_summary.csv
	else
		echo 0 >> ${dirName}_summary.csv # if parcel was not reported in the literature
  	fi
done < $HOME/fMRI/Glasser_labels.csv
############

# create a label file to visualize which areas were reported in the literature
# first comes the header, then for each label you add the vertices (and their coords) associated with that label from the Glasser label file
echo "#!ascii label, from condition ${dirName}" >> ${dirName}_lh_atlas.label
echo "#!ascii label, from condition ${dirName}" >> ${dirName}_rh_atlas.label
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



############ Create overlays ##########

# Unlike labels, these overlays contain vertex-specific information (specifically, proportion of reports across studies)
## Note: the contents of *_*h_overlay.csv have to be pasted into templateOverlay.gii
## Note on note: there has to be a better way to do that.. it's just a .txt -> xml conversion that FS must have a function for

echo
echo "Creating overlays for ${dirName}..."
echo

touch ${dirName}_lh_overlay.csv
touch ${dirName}_rh_overlay.csv

while IFS=$'\t' read -r -a label ; do
    if [ -z "$label" ] ; then
        #echo "Cell is empty"
        echo 0 >> ${dirName}_lh_overlay.csv
        #echo 0 >> check.csv
    elif [[ $(grep "$label" ${dirName}_summary.csv) != "" ]] ; then
        #echo "$label in summary"
        temp=$(grep $label ${dirName}_summary.csv)
        count=(${temp//,/ })
        echo ${count[1]} >> ${dirName}_lh_overlay.csv
        #echo ${count[0]} >> check.csv
    else
        #echo "Label not in summary"
        echo 0 >> ${dirName}_lh_overlay.csv
        #echo 0 >> check.csv
    fi
done < $HOME/fMRI/Glasser_labels_lh.csv

while IFS=$'\t' read -r -a label ; do
    if [ -z "$label" ] ; then
        #echo "Cell is empty"
        echo 0 >> ${dirName}_rh_overlay.csv
        #echo 0 >> check.csv
    elif [[ $(grep "$label" ${dirName}_summary.csv) != "" ]] ; then
        #echo "$label in summary"
        temp=$(grep $label ${dirName}_summary.csv)
        count=(${temp//,/ })
        echo ${count[1]} >> ${dirName}_rh_overlay.csv
        #echo ${count[0]} >> check.csv
    else
        #echo "Label not in summary"
        echo 0 >> ${dirName}_rh_overlay.csv
        #echo 0 >> check.csv
    fi
done < $HOME/Basic_files/Glasser_labels_rh.csv




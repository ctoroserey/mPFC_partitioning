#!/bin/bash

# This script takes the summary of the permutation analysis, and assigns the corresponding (signed) chi squared value that came from the comparison between conditions.

# loop through all summary files
for comparison in *summary.csv ; do

    fname=$(echo $comparison | cut -d '_' -f 1)
    touch ${fname}_lh_overlay.csv
    touch ${fname}_rh_overlay.csv

    echo "Working on ${fname} lh..."

    # Go through each label from Glasser and insert its chi2 value on the respective vertex
    # Change the label system if needed (i.e. DMN-space analysis).
    while IFS=$'\t' read -r -a label ; do

        if [ -z "$label" ] ; then
            #echo "Cell is empty"
            echo 0 >> ${fname}_lh_overlay.csv
            #echo 0 >> check.csv
        elif [[ $(grep "$label" ${fname}_summary.csv) != "" ]] ; then
            #echo "$label in summary"
            temp=$(grep $label ${fname}_summary.csv)
            count=(${temp//,/ })
            echo ${count[1]} >> ${fname}_lh_overlay.csv
        else
            #echo "Label not in summary"
            echo 0 >> ${fname}_lh_overlay.csv
            #echo 0 >> check.csv
        fi

    done < Glasser_labels_lh.csv

    echo "Working on ${fname} rh..."

    while IFS=$'\t' read -r -a label ; do

        if [ -z "$label" ] ; then
            #echo "Cell is empty"
            echo 0 >> ${dirName}_rh_overlay.csv
            #echo 0 >> check.csv
        elif [[ $(grep "$label" ${fname}_summary.csv) != "" ]] ; then
            #echo "$label in summary"
            temp=$(grep $label ${fname}_summary.csv)
            count=(${temp//,/ })
            echo ${count[1]} >> ${fname}_rh_overlay.csv
            #echo ${count[0]} >> check.csv
        else
            #echo "Label not in summary"
            echo 0 >> ${fname}_rh_overlay.csv
            #echo 0 >> check.csv
        fi

    done < Glasser_labels_rh.csv

  done

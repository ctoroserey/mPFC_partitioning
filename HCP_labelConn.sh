#!bin/bash/

# This is based on Sean's suggestions (see his NeuroImage paper, 2017)
# First, grab an MGSR dtseries.nii data file from an HCP subject (from SCC). This file contains the residuals from regressing the mean greyordinate value from the series.
# Second, get the mean of the timeseries for a given parcel using 'Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii' from Glasser's paper
# Note: using -cifti-convert -to-text on the resulting file will generate a csv with the mean time series (column) for each parcel (row, N=360).
#       This can be used to sed-vertex correlations, once the vertex index for columns have been sorted out
# Third, correlate the parcellated ptseries file. This will create a 360x360 pconn file.
# Fourth, convert the pconn file to text (as in the note) to get the 360x360 correlation matrix as a csv
# The advantage here is that the RS are well preprocessed, allowing us to bypass the merging of the 4 RS scans.

dataName=$1

# Parcellate
for dataFile in *${dataName} ; do

  filename=$(echo $dataFile | cut -f 1 -d '_')

  echo "Working on ${filename}"

  # wb_command -cifti-parcellate ${filename}_rfMRI_REST1+2_LR+RL.residMGSR.dtseries.nii \
  #                               Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii \
  #                               COLUMN ${filename}.ptseries.nii -method MEAN

  echo "Parcellating time series"
  wb_command -cifti-parcellate ${dataFile} \
                                Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii \
                                COLUMN ${filename}.ptseries.nii -method MEAN

  # Get each parcel's mean timeseries
  wb_command -cifti-convert -to-text ${filename}.ptseries.nii ${filename}_ptSeries.txt

  # Get the overall time series
  echo "Storing the whole time series"
  wb_command -cifti-convert -to-text ${dataFile} ${filename}_timeSeries.csv

done

# series=$(expr $series + 1)

# done


# z-score

# average across participants


##------ Notes (before 2/19/18)

# 1) The current issue is that the label for 7m is not working. I'm trying to use -cifti-label-to-roi, but it's not taking rh.R_7m_ROI.label in
#    Try with regular -label commands instead, or extracting from the official Glasser scene for the 1.2k release.
#
# 2) try using -cifti-convert -to-text *.dconn.nii adjacency.txt when done
#
# 3) This worked: wb_command -cifti-correlation merged.dtseries.nii test.dconn.nii \
#                            -roi-override -cifti-roi Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii
#
#    Problem is that it did so for all labels in the cifti parcellation. Still have to work out a way to do just 7m. Resulting file is ~30GB, and transforming to text
#    using point number 2 makes it a multi gig text file. Might be able to extract relevant info from that.

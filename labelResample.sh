#!/bin/bash

# 1. To transform the label data from Glasser et al. (2016) into fsaverage space:
# a. download the following dataset: http://brainvis.wustl.edu/workbench/standard_mesh_atlases.zip
#    steps from: https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP.pdf?version=1&modificationDate=1471988055831&api=v2
#    you also need to grab the Glasser label files

# b. separate the .dlabel.nii file into .gii for both hemispheres:
wb_command -cifti-separate Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN  \
                           -label CORTEX_LEFT Glasser_pre_labels_lh.label.gii  \
                           -label CORTEX_RIGHT Glasser_pre_labels_rh.label.gii  \

# c. resample the labels into fsaverage space (twice, one for each hemisphere):
wb_command -label-resample Glasser_pre_labels_lh.label.gii  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fsaverage6_std_sphere.L.41k_fsavg_L.surf.gii \
                           ADAP_BARY_AREA \
                           Glasser_post_labels_lh.label.gii  \
                           -area-metrics  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fsaverage6.L.midthickness_va_avg.41k_fsavg_L.shape.gii
mris_convert --annot $HOME/git_clones/fMRI/Glasser_post_labels_lh.label.gii $SUBJECTS_DIR/fsaverage/surf/lh.white $HOME/git_clones/fMRI/lh.glasser_aparc.annot


wb_command -label-resample Glasser_pre_labels_rh.label.gii  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fsaverage6_std_sphere.R.41k_fsavg_R.surf.gii \
                           ADAP_BARY_AREA \
                           Glasser_post_labels_rh.label.gii  \
                           -area-metrics  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii  \
                           $HOME/Documents/HBM/Glasser_et_al_2016_HCP_MMP1.0_RVVG/standard_mesh_atlases/resample_fsaverage/fsaverage6.R.midthickness_va_avg.41k_fsavg_R.shape.gii
mris_convert --annot $HOME/git_clones/fMRI/Glasser_post_labels_rh.label.gii $SUBJECTS_DIR/fsaverage/surf/rh.white $HOME/git_clones/fMRI/rh.glasser_aparc.annot

freeview -f $SUBJECTS_DIR/fsaverage6/surf/lh.inflated:annot=$PWD/Glasser_post_labels_lh.label.gii \
         -f $SUBJECTS_DIR/fsaverage6/surf/rh.inflated:annot=$PWD/Glasser_post_labels_rh.label.gii

#!/bin/bash

freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:annot=lh.HCP_MMP1p0.annot:overlay=DECvsDMN_lh_overlay.gii:overlay=DECvsNEG_lh_overlay.gii:overlay=DECvsPOS_lh_overlay.gii:overlay=DMNvsNEG_lh_overlay.gii:overlay=DMNvsPOS_lh_overlay.gii:overlay=NEGvsPOS_lh_overlay.gii \
            $SUBJECTS_DIR/fsaverage/surf/rh.inflated:annot=rh.HCP_MMP1p0.annot:overlay=DECvsDMN_rh_overlay.gii:overlay=DECvsNEG_rh_overlay.gii:overlay=DECvsPOS_rh_overlay.gii:overlay=DMNvsNEG_rh_overlay.gii:overlay=DMNvsPOS_rh_overlay.gii:overlay=NEGvsPOS_rh_overlay.gii

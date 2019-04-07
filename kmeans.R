#!/usr/bin/env Rscript
SubjID <- 101006
nClusts <- 3
thresh <- T

library(data.table)

# To store files ready to be converted to CIFTI
HCPOut <- function(Data = dmnval7mCommunities[[1]], MOI = "Membership", SubjID = "100307", padding = 0, n = 12){
  
  nVertices <- 59412
  tempVec <- rep(padding, nVertices)
  temp <- grep(MOI, colnames(Data))
  tempVec[Data$Vertex] <- Data[[temp]]
  write.table(file = paste(SubjID,'_Kmeans_', n, 'clusters_dataforCifti.txt', sep=""), tempVec, row.names = F, col.names = F, dec = ".")
  
}

# load data
setwd('./tseries')
temp <- list.files(pattern = paste(SubjID, "_Session*", sep=""))
timeSeries_sess <- lapply(temp, fread, header=F)
timeSeries_sess <- lapply(timeSeries_sess, data.matrix)

tseries <- do.call(cbind, timeSeries_sess)

setwd('../BB_comparison/')

if (thresh) {
  # ROIs
  yeoLbls <- c("L_25_ROI",
               "L_OFC_ROI",
               "L_10v_ROI",
               "R_25_ROI",
               "R_OFC_ROI",
               "R_10v_ROI",
               "L_s32_ROI",
               "L_RSC_ROI",
               "R_RSC_ROI",
               "R_23d_ROI",
               "R_d23ab_ROI",
               "R_31a_ROI",
               "R_31pv_ROI",
               "R_31pd_ROI",
               "R_7m_ROI",
               "R_v23ab_ROI",
               "R_p24_ROI",
               "R_d32_ROI",
               "R_9m_ROI",
               "R_p32_ROI",
               "R_a24_ROI",
               "R_10r_ROI",
               "R_10d_ROI",
               "L_23d_ROI",
               "L_d23ab_RO",
               "L_31a_ROI",
               "L_31pv_ROI",
               "L_31pd_ROI",
               "L_7m_ROI",
               "L_v23ab_ROI",
               "L_p24_ROI",
               "L_d32_ROI",
               "L_9m_ROI",
               "L_p32_ROI",
               "L_a24_ROI",
               "L_10r_ROI",
               "L_10d_ROI",
               "R_s32_ROI",
               "R_9a_ROI",
               "L_PCV_ROI",
               "L_POS1_ROI", # from this one on, it's Lauren D's extended space
               "L_PreS_ROI",
               "L_H_ROI",
               "L_POS2_ROI",
               "L_7Pm_ROI",
               "L_PGs_ROI",
               "L_PGi_ROI",
               "L_TPOJ2_ROI",
               "L_TPOJ3_ROI",
               "L_PGp_ROI",
               "L_TE2a_ROI",
               "L_MST_ROI",
               "L_MT_ROI",
               "L_9-46d_ROI",
               "L_8Ad_ROI",
               "L_46_ROI",
               "R_PGp_ROI",
               "R_PGs_ROI",
               "R_PGi_ROI",
               "R_TPOJ3_ROI",
               "R_TPOJ2_ROI",
               "R_MST_ROI",
               "R_MT_ROI",
               "R_TE2a_ROI",
               "R_TE1a_ROI",
               "R_9-46d_ROI",
               "R_46_ROI",
               "R_8Ad_ROI",
               "R_POS1_ROI",
               "R_POS2_ROI",
               "R_7Pm_ROI",
               "L_23c_ROI",
               "R_23c_ROI",
               "R_PreS_ROI",
               "R_H_ROI",
               "L_STSdp_ROI",
               "L_STSvp_ROI",
               "L_STSva_ROI",
               "L_TE1m_ROI",
               "L_TE1a_ROI",
               "L_TPOJ1_ROI",
               "L_8Av_ROI",
               "L_8C_ROI",
               "L_PHA2_ROI",
               "L_44_ROI",
               "L_PSL_ROI",
               "R_STV_ROI",
               "R_TPOJ1_ROI",
               "R_TE1m_ROI",
               "R_STSdp_ROI",
               "R_STSda_ROI",
               "R_8Av_ROI",
               "R_8C_ROI",
               "R_PCV_ROI",
               "R_PHA1_ROI",
               "R_p47r_ROI")   
  
  # load labels and ensure they are numeric
  labelCoords_vertex <- read.table('labelCoords_vertex.csv', sep = ",", header = T)
  
  # get ROI indices
  indx <- sapply(yeoLbls, function(lbl) {grep(lbl, labelCoords_vertex$Label)}) 
  indx <- do.call(c, indx)
  
  # reduce to this search space
  tseries <- tseries[indx, ]
  
  # run K-means as a clustering control (no seeding, this is just proving a point)
  km <- kmeans(tseries, nClusts, iter.max = 100)$cluster
  
  # df to write out in HCP format
  df <- data.frame(Vertex = indx,
                   kmeans = km)
  
  # write to text files ready to be transformed to CIFTIs
  HCPOut(df, MOI = "kmeans", SubjID = SubjID, padding = 0, n = nClusts)
  
} else {
  # run K-means as a clustering control (no seeding, this is just proving a point)
  km <- kmeans(tseries, nClusts, iter.max = 100)$cluster
  
  # store the labels
  write.table(file = paste(SubjID, "Kmeans_fullBrain.txt", sep = "_"), km, row.names = F, col.names = F, dec = ".")
}












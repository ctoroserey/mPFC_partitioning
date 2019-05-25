#!/usr/bin/env Rscript

# IMPORTANT: THIS FUNCTION ASSUMES THAT THE TIME SERIES FOR THE SUBJECT ARE HOUSE
# IN THE SAME DIRECTORY AS A CSV FILE, CALLED '<SUBJID>_timeSeries.csv``
SubjID <- commandArgs(trailingOnly = TRUE)
write(paste("Analyzing data for subject", SubjID), stdout())

library(igraph)
library(data.table)
library(parallel)

## tanh-z transformation (variance stabilizing Fisher) and p-values (adjusted and not)
# This takes either a matrix of correlation values (vectors too, but manually compute pvals)
# Normalization approach suggested in network textbook (equation 7.20)
# This transformation is approximately normal with mean 0 and sd = sqrt(1 / (n - 3))
fisherTanh <- function(Data = padjMatrix, preThresh = NA){
  
  transformed <- list()
  
  if (is.numeric(preThresh)) {
    Data[Data == 1] <- 0.999
    Data[Data < preThresh] <- 0
    transformed$tanhZ <- 0.5 * log((1 + Data) / (1 - Data))
  } else {
    # tanh
    Data[Data == 1] <- 0.999
    transformed$tanhZ <- 0.5 * log((1 + Data) / (1 - Data))
    
    # p-vals
    if (is.matrix(Data)) {
      z.vec <- transformed$tanhZ[upper.tri(transformed$tanhZ)]
      n <- dim(Data)[1]
    } else if (is.vector(Data)) {
      z.vec <- transformed$tanhZ
      n <- length(Data)
    }
    transformed$pvals <- 2 * pnorm(abs(z.vec), 0 , sqrt(1 / (n-3)), lower.tail=F)
    
    # adjust pvals
    transformed$adjustPvals <- p.adjust(transformed$pvals, "BH")
    
    if (is.matrix(Data)) {
      # get pvals and their adjustment into a symetric matrix form
      # regular
      tempMat <- matrix(0, dim(Data)[1], dim(Data)[2])
      tempMat[upper.tri(tempMat)] <- transformed$pvals
      tempMat[lower.tri(tempMat)] <- transformed$pvals
      dimnames(tempMat) <- list(rownames(Data), rownames(Data))
      transformed$pvals <- tempMat
      
      # adjusted
      tempMat <- matrix(0, dim(Data)[1], dim(Data)[2])
      tempMat[upper.tri(tempMat)] <- transformed$adjustPvals
      tempMat[lower.tri(tempMat)] <- transformed$adjustPvals
      dimnames(tempMat) <- list(rownames(Data), rownames(Data))
      transformed$adjustPvals <- tempMat
      
      # retain only the significant adjusted pval Z scores
      transformed$tanhZ[transformed$adjustPvals > 0.05] <- 0
    }
  }
  
  return(transformed)
  
}

# To store files ready to be converted to CIFTI
HCPOut <- function(Data = dmnval7mCommunities[[1]], MOI = "Membership", SubjID = "100307", padding = 0){
  
  nVertices <- 59412
  tempVec <- rep(padding, nVertices)
  temp <- grep(MOI, colnames(Data))
  tempVec[Data$Vertex] <- Data[[temp]]
  write.table(file = paste(SubjID,"_",MOI,'_dataforCifti.txt', sep=""), tempVec, row.names = F, col.names = F, dec = ".")
  
}

# if a minimum corr value is desired prior to clustering (otherwise set NA for p-value based thresholding)
corrThresh <- NA


# ROIs used on the paper
ROI <- c("L_25_ROI",
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
             "L_PCV_ROI")   

# load labels and ensure they are numeric
labelCoords_vertex <- read.table('labelCoords_vertex.csv', sep = ",", header = T)

# get ROI indices
indx <- sapply(ROI, function(lbl) {grep(lbl, labelCoords_vertex$Label)}) 
indx <- do.call(c, indx)

## load and prep data
write("Loading data", stdout())
Data <- fread(paste(SubjID, "_timeSeries.csv", sep = ""), header = F)
Data <- data.matrix(Data)

# append the ROI label to each vertex
dimnames(Data) <- list(labelCoords_vertex$Label, seq(ncol(Data)))

# select ROIs
indx <- sapply(ROI, function(lbl) {grep(lbl, rownames(Data))}) 
indx <- do.call(c, indx)

# reduce time series to only search space
nVerts <- length(indx)
tseries <- Data[indx, ]

## analyze
# correlate activity among all vertices, fisher transform, and filter out non-sig corrs.
write("Correlating data", stdout())
corrMat <- cor(t(tseries))
transfMat <- fisherTanh(Data = corrMat, preThresh = corrThresh)
corrMat <- transfMat$tanhZ

diag(corrMat) <- 0 # avoid self-loops when creating the graph
corrMat <- exp(corrMat) # ensure all connections are positive while maintaining the ordinal ranks
corrMat[corrMat == 1] <- 0 # return filtered-out vertices to 0

# transform network weight matrix to a graph object
tempGraph <- graph_from_adjacency_matrix(corrMat, weighted = T, mode = "undirected")


# Run a fastgreedy modularity community detection on ROIs
# this is usually performed as a validity check
write("Running modularity", stdout())
tempCommunity <- fastgreedy.community(tempGraph)


# now spectral partitioning
write("Running SP", stdout())
tempLap <- laplacian_matrix(tempGraph, normalized = T)
tempEigen <- eigen(tempLap)
eigenVals <- tempEigen$values
eigenVecs <- tempEigen$vectors

# store all outputs in this data frame, so they can be written afterwards in HCP format  
df <- data.frame(Vertex = indx,
                 eigenVal = eigenVals,
                 FV = eigenVecs[ , ncol(eigenVecs) - 1],
                 modularity = modularity$V1)

# write to text files ready to be transformed to CIFTIs
# these can be turned into CIFTI files using wb_command:
# wb_command -cifti-convert -from-text <HCPOut output> <templace CIFTI.dscalar.nii> <new CIFTI.dscalar.nii>
HCPOut(df, MOI = "FV", SubjID = SubjID, padding = 0)
HCPOut(df, MOI = "modularity", SubjID = SubjID, padding = -1)

# store the dtaframe just in case (ARI comparisons, for ex)
write.csv(df, file = paste(SubjID, "commDetection.csv", sep = "_"), sep = ",")


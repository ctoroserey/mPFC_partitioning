#!/usr/bin/env Rscript
SubjID <- commandArgs(trailingOnly=TRUE)
write(paste("Analyzing data for subject", SubjID), stdout())

### Run a graph theoretic fMRI analysis on the pre-processed time series of an individual

###------- Libraries and functions ---------------
write('Loading functions and libraries...', stdout())

library(igraph)
#library(tidyverse)
library(ggnetwork)
library(ggplot2)
library(gridExtra)
library(data.table)
library(plyr)
library(mcclust)
library(lme4)
library(parallel)


# tanh-z transformation (variance stabilizing Fisher) and p-values (adjusted and not)
fisherTanh <- function(Data = padjMatrix, preThresh = NA){
  # This takes either a matrix of correlation values (vectors too, but manually compute pvals)
  # Normalization approach suggested in network textbook (equation 7.20)
  # This transformation is approximately normal with mean 0 and sd = sqrt(1 / (n - 3))
  
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
    transformed$pvals <- 2 * pnorm(abs(z.vec), 0 , sqrt(1 / (n - 3)), lower.tail = F)
    
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
      
    }
    
    # retain only the significant adjusted pval Z scores
    transformed$tanhZ[transformed$adjustPvals > 0.05] <- 0
  }
  
  return(transformed)
  
}


# get ROI coords & index
getCoords <- function(Labels = DMN_labels, Coords = labelCoords_parcel, TimeSeries = FALSE){
  # The point here is to reduce the summary dframes from community detection to show only ROIs
  # Should work for extracting any label-indexed dframe though
  # I wanted to also get the index in case I want to extract specific rows from parcel/vertex coord dframes
  
  indx <- numeric()
  
  # If you want to select time series from raw data
  if (TimeSeries) {
    for (ROI in Labels) {
      indx <- c(indx, grep(ROI, rownames(Coords)))
    }
  } else { # for the summary output of the community detection output
    for (ROI in Labels) {
      indx <- c(indx, grep(ROI, Coords$Label))
    }
  }
  
  results <- list()
  results$Index <- indx
  results$Coords <- Coords[indx, ]
  
  return(results)
  
}


# Run community detection on ROIs
communityDetection <- function(Data = NA, ROIS = "None", modularity = T, extras = T) {
  # This function relies on having the timeSeries data uploaded, and labelCoords_vertex as the row names for ROI selection
  # Extras dictates whether the community object + correlation matrix should also be extracted
  #
  # Parameter definitions:
  #
  #   Data: raw time series data
  #
  #   ROIS: a list. If the data needs to be reduced to specific ROIs  
  #
  #   modularity: TRUE if you would like to also run the fastgreedy modularity community detection
  #
  #   extras: if you want to return the original and transformed correlation matrices  

  if (ROIS != "None") {
    # To store the vertex indices corresponding to the ROIs
    indx <- numeric()
    
    # had to place dashes on each side because grep grabbed strings containing the names (i.e. 47m, a24pr)
    for (ROI in ROIS) {
      indx <- c(indx, grep(ROI, rownames(Data)))  
    }
    
    # reduce time series to only include ROIs
    nVerts <- length(indx)
    Data <- Data[indx, ]
  }
  
  # This used to be done with the for loop, but it was too slow. cor() speeds up the process by a lot
  corrMat <- cor(t(Data))
  
  # name the dimensions of the matrix according to the surface vertex index
  rownames(corrMat) <- indx
  
  # store untransformed correlation matrix for later
  corrMatrix <- corrMat
  
  # transform to Fisher's (think of thresholding)
  transfMat <- fisherTanh(Data = corrMat)
  
  # Store Fisher transformed vals for graphing
  corrMat <- transfMat$tanhZ
  
  # diagonals of 1 could be interpreted as self-loops
  diag(corrMat) <- 0
  
  # Exponentiate to preserve distribution while ensuring positive weights
  corrMat <- exp(corrMat)
  corrMat[corrMat == 1] <- 0
  
  # community detection
  # I initially used the absolute value of the correlation, but the exp preserves the distribution 
  # create a graph from adjacency matrix
  tempGraph <- graph_from_adjacency_matrix(corrMat, weighted = T, mode = "undirected")
  
  # spectral partitioning
  tempLap <- laplacian_matrix(tempGraph, normalized = T)
  tempEigen <- eigen(tempLap)
  fvec <- tempEigen$vectors[, length(tempEigen$values) - 1]
  binarized <- as.factor(ifelse(fvec > 0, 1, 0)) # binarized Fiedler Vector 
  
  if (modularity) {
    # modularity
    tempCommunity <- fastgreedy.community(tempGraph)
    
    # put summary together
    summary <- data.frame(Vertex = indx,
                          Label = tempCommunity$names,
                          x = labelCoords_vertex[indx, "x"],
                          y = labelCoords_vertex[indx, "y"],
                          z = labelCoords_vertex[indx, "z"],
                          Hemisphere = substring(tempCommunity$names, 1, 1),
                          Membership = tempCommunity$membership,
                          Modularity = tempCommunity$modularity,
                          EigenVal = tempEigen$values,
                          FiedlerVec = fvec,
                          FiedlerBinary = binarized)
  } else {
    # put summary together
    summary <- data.frame(Vertex = indx,
                          Label = tempCommunity$names,
                          x = labelCoords_vertex[indx, "x"],
                          y = labelCoords_vertex[indx, "y"],
                          z = labelCoords_vertex[indx, "z"],
                          Hemisphere = substring(tempCommunity$names, 1, 1),
                          EigenVal = tempEigen$values,
                          FiedlerVec = fvec,
                          FiedlerBinary = binarized)
  }

  # Get the final components together
  # extras should be TRUE if you want to keep the original and transformed correlation matrices
  if (extras) {
    communityResults <- list(CorrMatrix = corrMatrix,
                             TransfMatrix = transfMat,
                             Summary = summary)
  } else {
    communityResults <- summary
  }
  
  return(communityResults)
  
}


# Permutation for 2 groups
permute <- function(group1 = 1, group2 = 2, statType = mean, nPerms = 5000, paired = FALSE){
  
  # prep data
  summaryPerm <- list()
  lOne <- length(group1)
  lTwo <- length(group2)
  bigSample <- c(group1,group2)  
  
  if (paired == FALSE) {
    
    
    for (i in 1:nPerms){
      
      # relabel samples
      tempBig <- sample(bigSample)
      tempOne <- tempBig[1:lOne]
      tempTwo <- tempBig[(lOne+1):length(bigSample)]
      
      # stats
      tempDiffs <- statType(tempOne,na.rm=T) - statType(tempTwo,na.rm=T)
      summaryPerm$jointDist[i] <- tempDiffs # statType(tempDiffs, na.rm = T) 
      
    }  
    
  } else {
    
    for (i in 1:nPerms){
      
      # shift labels in a pairwise fashion
      tempDiffs <- statType((-1)^rbinom(lOne,1,0.5) * (group1 - group2))
      summaryPerm$jointDist[i] <- tempDiffs
      
    }
    
  }
  
  # get the observed difference
  diffs <- statType(group1,na.rm=T) - statType(group2,na.rm=T)
  observedAbs <- abs(diffs) # maybe leave it as means here
  observed <- diffs
  summaryPerm$Pval <- 2 * (1 - ecdf(summaryPerm$jointDist)(observedAbs))
  if (length(unique(abs(summaryPerm$jointDist))) == 1) {summaryPerm$Pval <- 1} # if the difference is always the same, then p = 1
  summaryPerm$Observed <- observed
  
  return(summaryPerm)
  
}


# Non-parametric Bootstrap for a single group
bootstrap <- function(group = 1, statType = mean, B = 5000){
  
  # prep param
  bootStats <- rep(0,B)
  
  # iterate
  for(b in 1:B){
    
    # wait group
    x <- sample(group,length(group),replace=T)  
    bootStats[b] <- statType(x,na.rm = T)
    
  }
  
  return(bootStats)
  
}


# Sliding window
slidingWindow <- function(subjTS = NA, mins = 15, jump = 1) {
  # This is a fairly specific function. It takes the time series from a participant and preps/runs community detection at each specified time window. 
  # Returns the summaries for each window (based on communityDetection function)
  # The selection of a window size is based on the fact that a TR = 0.720s and an hour is 5000 TRs (HCP-based)
  # Since each subject has a different amount of time points, and all are slightly under 1 hr, I chose to round down the number of window moves to avoid unevenness
  # 
  # Parameter definitions:
  #   
  #   subjTS: Subject time series. Note: if you want to look at more than vmpfc and 7m, add a parameter for the labels that can be passed to the community detection function
  # 
  #   mins: Size of the window that will slide through the data.
  #   
  #   jump: Steps (in mins) advanced per slide
  #   
  #   Spectral: Whether to compute spectral partitioning as well. Takes significantly longer, but might be more useful for bisections
  #
  # Right now this takes ~27 mins per subject to run. Think of ways to improve that.
  
  
  # Adapted so it works with the output from getCoords
  # One wouldn't really apply a full-brain analysis of this sort anyways, too computationally intensive
  indx <- subjTS$Index
  subjTS <- subjTS$Coords
  
  # Window sizing (length)
  # Think about incorporating custom TRs
  TS <- ncol(subjTS) # time series for the subject
  WS <- seq(834 * (mins / 10)) # window from the first TR up to 834 (~10 mins) times the desired multiplier
  jump <- 84 * jump # 84 ~ 1 min, times the number of mins that the window moves
  nJumps <- floor((TS - length(WS)) / jump) # number of jumps to be performed, based on the selected parameters
  
  ##------- using lapply
  winData <- mclapply(seq(nJumps), function(x) {subjTS[, WS + (jump * (x - 1))]})
  commTS <- mclapply(winData, communityDetection, ROIS = "None", modularity = F, extras = F)
  
  return(commTS)
  
}


# compare sliding window data to overall communities
slideCompare <- function(subjData = slideCommunities[[1]], template = NA, func = "RI", comm = "Spectral") {
  # This function compares the community partition from each window slide to the one derived from the whole data set
  # 
  # Parameter definitions:
  #   
  #   subjData: A participant's output from slideCommunities
  #   
  #   template: The partition from the whole time series
  # 
  #   func: Which function to use for comparing ("RI" for adjusted RI, "VI" for variation of information)
  # 
  #   comm: The subject data might contain modularity and spectral partitions. Choose which to use.
  
  # alternative way:
  tempComparison <- sapply(subjData, function(data) {arandi(data$FiedlerBinary, template$FiedlerBinary, adjust = T)})
  
  # For storing RIs or VIs
  tempComparison <- numeric()
  
  # How many jumps does the original data contain?
  nJumps <- length(subjData)
  
  # This loop is technically backwards. I should technically divide by function, then partition method, and then run the window comparisons
  # It's still fast, so I won't worry.

    if (comm == "Modularity") {
      if (func == "RI") {
        tempComparison[Win] <- arandi(subjData[[Win]]$Membership, template$Membership, adjust = T)
      } else if (func == "VI") {
        tempComparison[Win] <- vi.dist(subjData[[Win]]$Membership, template$Membership)
      }
    } else if (comm == "Spectral") {
      if (func == "RI") {
        tempComparison <- sapply(subjData, function(data) {arandi(data$FiedlerBinary, template$FiedlerBinary, adjust = T)})
      } else if (func == "VI") {
        tempComparison[Win] <- vi.dist(subjData[[Win]]$FiedlerBinary, template$FiedlerBinary)
      }
    }
  
  return(tempComparison)
  
}


# Ensure that all spectral communities associated with 7m (i.e. DMN) have the same label value of 1 
evenSpectral <- function(Data = slideCommunities[[1]][[7]]) {
  
  # Get only the values for 7m
  shortData <- Data[grep("_7m_", Data$Label), ]
  
  # Get 7m's most probable affiliation
  affil <- mean(as.numeric(shortData$FiedlerBinary)-1)
  
  # If it is close to 1, then invert the labeling. 
  # Note: This has no effect on RI or VI, since they are insensitive to actual labeling
  # This is for visualization purposes only
  if (affil < 0.5) {
    
    UD <- as.numeric(Data$FiedlerBinary) - 1
    Data$FiedlerBinary <- as.factor((UD - 1)^2)
    Data$FiedlerVec <- Data$FiedlerVec * -1
    
  }
  
  return(Data)
  
}


# Create a vector ready to be used for HCP data (32k CIFTI surface)
HCPOut <- function(Data = NA, MOI = "Membership", SubjID = "100307"){
  # The input should be the summary from community partitioning
  # Once this is created, go to the terminal and input something like this
  # wb_command -cifti-convert -from-text dataforCifti.txt 100307.MyelinMap_BC.32k_fs_LR.dscalar.nii testCifti.dscalar.nii
  # Where the myelin file here is just a templace. It can be any dscalar.nii with the right surface size
  
  nVertices <- 59412
  tempVec <- rep(-1, nVertices)
  temp <- grep(MOI, colnames(Data))
  tempVec[Data$Vertex] <- Data[[temp]]
  write.table(file = paste(SubjID,"_",MOI,'_dataforCifti.txt', sep=""), tempVec, row.names = F, col.names = F, dec = ".")
  
}


# Perform pairwise comparisons of clustering outcomes on all subjects
comparePartitions <- function(Data = NA, MOI = "FiedlerBinary", Index = "VI", nSubjects = nSubj, subjNames = subjList) {
  # This function will compare the community partitions from all subjects and create a 'comparison matrix' for every pairwise combination of subjects
  # Alternatively, if a second MOI is added 
  # Inputs
  # 
  # Data: the list of summaries produced by the script
  # 
  # MOI: measure of interest (usually the binarized Fiedler vector). If a vector, compares across algorithms per subject
  # 
  # Index: VI for variation of information, RI for the adjusted rand index
  
  # Get the column position of the MOI(s)
  Columns <- colnames(Data[[1]])
  indx <- as.numeric(Columns %in% MOI)
  MOI_indx <- which(indx==1)
  
  if (length(MOI_indx) < 2) {
    # Combine the measures of interest
    allVecs <- do.call(cbind, lapply(Data, "[[", MOI))
    
    # Create empty matrix
    indexMatrix <- matrix(data = NA, nrow = nSubjects, ncol = nSubjects)
    dimnames(indexMatrix) <- list(subjNames, subjNames)
    
    # Run every pairwise comparison with the index of interest on the measure of interest
    if (Index == "VI") {
      for (subj in seq(nSubjects)) {
        for (subj2 in seq(nSubjects)) {
          indexMatrix[subj, subj2] <- vi.dist(allVecs[, subj], allVecs[, subj2])
        }
      }
    } else if (Index == "RI") {
      for (subj in seq(nSubjects)) {
        for (subj2 in seq(nSubjects)) {
          indexMatrix[subj, subj2] <- arandi(allVecs[, subj], allVecs[, subj2], adjust = T)
        }
      }
    }
  } else {
    indexMatrix <- data.frame(SubjID = as.character(subjList),
                              Index = rep(0, nSubjects))
    # Run every pairwise comparison with the index of interest on the measure of interest
    if (Index == "VI") {
      for (subj in seq(nSubjects)) {
        indexMatrix$Index[subj] <- vi.dist(Data[[subj]][, MOI_indx[1]], Data[[subj]][, MOI_indx[2]])
      }
    } else if (Index == "RI") {
      for (subj in seq(nSubjects)) {
        indexMatrix$Index[subj] <- arandi(Data[[subj]][, MOI_indx[1]], Data[[subj]][, MOI_indx[2]], adjust = T)
      }
    }
  }
  
  return(indexMatrix)
  
}


###------- Setups ---------------
write('Setting up variables for future computations...', stdout())

## Choose what to analyze
Vertex <- T # ROI-based vertex analysis
Sliding <- T

# colors to differentiate other things
#Cols <- c("aquamarine4","#D9541A",rgb(190,190,190,100, maxColorValue = 255)) # left, right, interhemisphere
Cols <- c("DMN" = "aquamarine4","Valuation" = "#D9541A",rgb(190,190,190,100, maxColorValue = 255)) # left, right, interhemisphere

# Label groups for specific analyses (DMN and DEC)
# These come from an arbitrary division based on where Mackey & Petrides (2014) draw their line
# This is computed below, but if we want to save time it can be loaded instead.
#bothLbls <- as.character(read.table('bothLbls.csv', header = F)$V2)

# Based on Yeo et al.
bothLbls <- c("L_25_ROI",
              "L_OFC_ROI",
              "L_10v_ROI",
              "R_25_ROI",
              "R_OFC_ROI",
              "R_10v_ROI",
              "L_s32_ROI",
              "L_RSC_ROI",
              "R_RSC_ROI",
              "R_8Ad_ROI",
              "R_9p_ROI",
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
              "L_8Ad_ROI",
              "L_9p_ROI",
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
              "L_p10p_ROI",
              "L_PCV_ROI")            

PCC_labels <- c("_7m_",
                    "_23d_",
                    "_31a_",
                    "_31pd_",
                    "_31pv_",
                    "_PCV_",
                    "_POS1_",
                    "_RSC_",
                    "_d23ab_",
                    "_v23ab_")

mPFC_labels <- c("_a24_",
               "_d32_",
               "_p32_",
               "_10r_",
               "_9m_",
               "_10v_",
               "_OFC_",
               "_s32_",
               "_p24_",
               "_10d_",
               "_25_")

# Little test to ensure that the ROIs match before running the whole thing
if (FALSE %in% (mPFC_labels %in% c(mPFC_labels, "Blah"))) {
  stop("mPFC-only and Overlapping ROIs don't match")
}

# Cortical surface for plotting
labelPerVertex <- read.table('labelPerVertex.csv', header = F)
labelCoords_vertex <- read.csv2('labelCoords_vertex.csv', sep = ",")[, 2:6]
labelCoords_parcel <- read.csv2('labelCoords_parcel.csv', sep = ",")

# Ensure that coordinates are numeric once loaded
labelCoords_vertex <- transform(labelCoords_vertex, x = as.numeric(as.character(x)))
labelCoords_vertex <- transform(labelCoords_vertex, y = as.numeric(as.character(y)))
labelCoords_vertex <- transform(labelCoords_vertex, z = as.numeric(as.character(z)))

###------- Load HCP Data ---------------
# Prep subject data from HCP
write("Loading vertex time series...", stdout())  
  
# Non-demeaned tseries for SNR calculations
# It is feasible to load these and demeaned them here, but I'm not sure if additional preprocessing took place post-demean (i.e. MGSR)
temp <- list.files(pattern = paste(SubjID, "_SNRSession*", sep=""))
timeSeries_SNR <- lapply(temp, fread, header=F)
timeSeries_SNR <- lapply(timeSeries_SNR, data.matrix)  

# Time series per session
temp <- list.files(pattern = paste(SubjID, "_Session*", sep=""))
timeSeries_sess <- lapply(temp, fread, header=F)
timeSeries_sess <- lapply(timeSeries_sess, data.matrix)

# Time series per day
timeSeries_halves <- list(cbind(timeSeries_sess[[1]], timeSeries_sess[[2]]),
                          cbind(timeSeries_sess[[3]], timeSeries_sess[[4]]))

# Concatenate
timeSeries <- do.call(cbind, timeSeries_sess)

# Populate each subject's column and row names
for (Win in seq_along(timeSeries_sess)) {
  dimnames(timeSeries_sess[[Win]]) <- list(labelCoords_vertex$Label, seq(ncol(timeSeries_sess[[Win]])))
  dimnames(timeSeries_SNR[[Win]]) <- list(labelCoords_vertex$Label, seq(ncol(timeSeries_SNR[[Win]])))
}
dimnames(timeSeries_halves[[1]]) <- list(labelCoords_vertex$Label, seq(ncol(timeSeries_halves[[1]])))
dimnames(timeSeries_halves[[2]]) <- list(labelCoords_vertex$Label, seq(ncol(timeSeries_halves[[2]])))
dimnames(timeSeries) <- list(labelCoords_vertex$Label, seq(ncol(timeSeries)))

rm(temp)



###------- Community Detection: all ---------------
write("Computing community detection for the ROI vertices...", stdout())

# compute spectral partitioning and modularity
communities <- communityDetection(Data = timeSeries, ROIS = bothLbls, modularity = T, extras = F)

# Append the tSNR
# communities$Summary$tSNR <- tSNR

# Ensure that DMN is assciated with the positive values of the eigenvector
communities <- evenSpectral(communities)


###------- Measures of variance -----------

write("Computing each vertex's variance...", stdout())

## Raw mPFC and PCC SNR
# Reduce the 4 SNR-dedicated tseries to the ROIs
timeSeries_SNR_ROI <- lapply(timeSeries_SNR, getCoords, Labels = bothLbls, TimeSeries = T)
rm(timeSeries_SNR)

# Compute tSNR
SNR_mean <- lapply(timeSeries_SNR_ROI, function(data) rowMeans(data$Coords))
SNR_sd <- lapply(timeSeries_SNR_ROI, function(data) apply(data$Coords, 1, sd))
tSNR <- do.call(cbind, lapply(seq_along(length(SNR_mean)), function(x) SNR_mean[[x]] / SNR_sd[[x]]))

# Get the mean coefficient across sessions, and append it to the overall summary
communities$tSNR <- rowMeans(tSNR)

rm(SNR_mean, SNR_sd, timeSeries_SNR_ROI)

## Now compute the variance explained by the mean time series from the DMN community on each DMN community vertex
## This should be high, since the communities are correlation-based

# Get the vertices from the DMN community
index <- communities$Vertex[communities$FiedlerBinary == 1]
tempY <- communities$y[communities$FiedlerBinary == 1]
tempZ <- communities$z[communities$FiedlerBinary == 1]

# Select their corresponding time series, and compute the average tseries of the whole community
DMN_commTSeries <- timeSeries[index, ] 
DMN_commTSeries_mean <- colMeans(DMN_commTSeries)

# Run a simple linear model per vertex time series
DMN_commTSeries_models <- apply(DMN_commTSeries, 1, function(vertexTimeseries) {lm(vertexTimeseries ~ DMN_commTSeries_mean)})

# Get the R-squares and T-values, and store them along with the vertex index
Rsquared <- sapply(DMN_commTSeries_models, function(data) {summary(data)$r.squared})
Tvalue <- sapply(DMN_commTSeries_models, function(data) {summary(data)$coefficients[2,3]})
DMN_commVrtxvsOverallCorr <- data.frame(index, tempY, tempZ, Rsquared, Tvalue)
colnames(DMN_commVrtxvsOverallCorr)[1] <- "Vertex"

rm(index, Rsquared, Tvalue, DMN_commTSeries, DMN_commTSeries_models, DMN_commTSeries_mean)
  
###------- Community Detection: Day 1 vs Day 2 -----------
  
write("Computing PCC-seed correlations and community detection for the ROI vertices, Day 1 vs Day 2...", stdout())

# Correlate Yeo 2011 PCC ROI with mPFC
# Get vertex indices
Yeo_PCC_indx <- getCoords(Labels = PCC_labels, Coords = timeSeries_halves[[1]], TimeSeries = T)$Index
mPFC_indx <- getCoords(Labels = mPFC_labels, Coords = timeSeries_halves[[1]], TimeSeries = T)$Index

# Correlate the mean PCC time series with mPFC
corHalves <- mclapply(timeSeries_halves, function(data) as.numeric(cor(t(data[mPFC_indx,]), colMeans(data[Yeo_PCC_indx,]))))

# Compute the day 1/2 communities 
communities_halves <- mclapply(timeSeries_halves, communityDetection, ROIS = bothLbls, modularity = T, extras = F)

# Ensure that the binarized partition is consistent across windows (i.e. ~7m = 1)\
communities_halves <- lapply(communities_halves, evenSpectral)

rm(timeSeries_halves)
  

###------- Community Detection: session ---------------
  
write("Computing community detection for the ROI vertices, session...", stdout())

# Modularity
communities_sess <- mclapply(timeSeries_sess, communityDetection, ROIS = bothLbls, modularity = T, extras = F)

# Ensure that the binarized partition is consistent across windows (i.e. ~7m = 1)\
communities_halves <- lapply(communities_halves, evenSpectral)

rm(timeSeries_sess)
  


###------- Vertex communities: Sliding window ---------------
if (Sliding) {
  
  write("Computing sliding window analysis for vertex communities...", stdout())
  
  # Prep data
  ROI_timeSeries <- getCoords(Coords = timeSeries, Labels = bothLbls, TimeSeries = T)
  
  # Sliding window analysis for each subject (using defaults)
  slideCommunities <- slidingWindow(subjTS = ROI_timeSeries, mins = 20, jump = 1)
  
  # Ensure that the binarized partition is consistent across windows (i.e. ~7m = 1)
  slideCommunities <- lapply(slideCommunities, evenSpectral)
  
  # Store the community structure per sliding window 
  slidingVals <- do.call(cbind, lapply(slideCommunities, "[[", "FiedlerBinary")) # do call will perform a function on a list. lapply creates a sublist from the original data
  slidingVals <- slidingVals - 1 # Lazy way of fixing the fact that the factorized FiedlerBinary becomes numeric as 1 and 2, thus making slideProp below messy (but fixable)

  # now compare to the template 
  # Get the index for each comparison
  tempCompare <- slideCompare(subjData = slideCommunities,
                              template = communities,
                              comm = "Spectral",
                              func = "RI")
  
  # Store as data frame 
  allSlideCompared <- data.frame(SubjID = rep(SubjID, length(tempCompare)),
                          Window = seq(length(tempCompare)),
                          Index = tempCompare)

  # Add proportion of times a node was associated with DNM to the final summary
  communities$slidePropDMN <- rowMeans(slidingVals)
    
}
###------- Save Data ---------------
  
write("Saving summary files...", stdout())

# Produce a final summary for each participant
write.csv(communities, paste(SubjID, "_finalSummary.csv", sep=""), row.names = F)

# Final summary for each half
lapply(seq_along(communities_halves), function(x) {write.csv(communities_halves[[x]], paste(SubjID, "_H", x,"_finalSummary.csv", sep = ""), row.names = F)})
lapply(seq_along(corHalves), function(x) {write.table(corHalves[[x]], paste(SubjID, "_H", x,"_corHalves.csv",sep=""), row.names = F, col.names = F)})

# Final summary for each session
lapply(seq_along(communities_sess), function(x) {write.csv(communities_sess[[x]], paste(SubjID, "_S", x,"_finalSummary.csv",sep=""), row.names = F)})

# The variance explained by the mean DMN time series on each vertex
write.csv(DMN_commVrtxvsOverallCorr, paste(SubjID, "_DMN_VrtxvMean.csv", sep=""), row.names = F)

# Save the agreement between each slide and the overall partition
if (Sliding) {
  write.csv2(allSlideCompared, file = paste(SubjID, "_slideWindowComparisons.csv", sep=""), row.names = F)

  # Write the actual sliding affiliations 
  write.csv2(slidingVals, file = paste(SubjID, "_slideWindowValues.csv", sep=""), row.names = F)
}

# Write files for HCP
HCPOut(Data = communities, MOI = "FiedlerBinary", SubjID = SubjID)  
HCPOut(Data = communities, MOI = "FiedlerVec", SubjID = SubjID) 
HCPOut(Data = communities, MOI = "tSNR", SubjID = SubjID) 
  













#!/usr/bin/env Rscript
SubjID <- commandArgs(trailingOnly=TRUE)
write(paste("Analyzing data for subject", SubjID), stdout())

### Run a graph theoretic fMRI analysis

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
#library(lmerTest)

## log-log degree distribution
logDegreeDist <- function(Graph = padjMat) {
  
  tempD <- degree(Graph)
  tempDD <- degree.distribution(Graph) # frequency of occurrencies of certain degrees
  d <- (0:(max(tempD)-1))
  ind <- (tempDD!=0)
  plot(d[ind], tempDD[ind], log = "xy", col = "blue", 
       xlab = "Log-Degree", ylab = "Log-Intensity",  
       main = "Log-Log Degree Distribution")
  
}


## Create correlation matrix between vertices of 2 ROIs
# Requires labeled time series matrix
vertexCorrMat <- function(ROI_1 = 'R_7m_ROI', ROI_2 = 'L_7m_ROI'){
  
  # extract the time series from the ROIs
  indexing_1 <- which(rownames(timeSeries)==as.character(ROI_1))
  indexing_2 <- which(rownames(timeSeries)==as.character(ROI_2))
  nVerts_1 <- length(indexing_1)
  nVerts_2 <- length(indexing_2)
  ROI_1_tseries <- timeSeries[indexing_1, ]
  ROI_2_tseries <- timeSeries[indexing_2, ]
  
  # Create empty matrix
  corrMat <- matrix(data = 0, nrow = nVerts_1, ncol = nVerts_2)
  
  # 
  for (vertex_1 in seq(nVerts_1)){
    
    # temp series 1
    tempOne <- ROI_1_tseries[vertex_1, ]
    
    for(vertex_2 in seq(nVerts_2)){
      
      # temp series 2
      tempTwo <- ROI_2_tseries[vertex_2, ]
      
      # store the Pearson correlation in the corr matrix
      corrMat[vertex_1, vertex_2] <- cor(tempOne, tempTwo, method = "pearson") 
      
    }
  }
  
  # name the dimensions of the matrix according to the surface vertex index
  rownames(corrMat) <- indexing_1
  colnames(corrMat) <- indexing_2
  
  # plot (optional?)  
  # unclustered heatmap
  heatmap(corrMat,
          Rowv = NA,
          Colv = NA,
          scale = "row",
          col = pallette(1000))
  
  # clustered heatmap
  heatmap(corrMat,
          scale = "row",
          col = pallette(1000))
  
  # print min-max corr vals
  print(range(corrMat))
  
  return(corrMat)
  
}


## Correlation matrix between an ROI and the Glasser parcels
parcelCorrMat <- function(ROI = 'R_7m_ROI', lbels = lookup$V1){
  
  # extract the time series from the ROIs
  indexing <- which(rownames(timeSeries)==as.character(ROI))
  nVerts <- length(indexing)
  ROI_tseries <- timeSeries[indexing, ]
  
  # create empty matrix to store values
  corrMat <- matrix(data = 0, nrow = nVerts, ncol = 360)
  
  # Loop through every combination of label and 7m vertex
  for (Parcel in lbels){
    
    # Get an ROI index to retrieve timeseries and store correlation vals
    indxROI <- which(lbels==as.character(Parcel))
    
    for(vertex in seq(nVerts)){
      
      # temp series
      tempROI <- as.numeric(ptSeries[indxROI, ])
      tempVertex <- ROI_tseries[vertex, ]
      
      # store the Pearson correlation in the corr matrix
      corrMat[vertex, indxROI] <- cor(tempROI, tempVertex, method = "pearson") 
      
    }
  }
  
  # name the columns according to labels
  colnames(corrMat) <- lbels
  rownames(corrMat) <- indexing
  
  # plot (optional?)  
  # unclustered heatmap
  heatmap(corrMat,
          Rowv = NA,
          Colv = NA,
          scale = "none",
          col = pallette(1000))
  
  # clustered heatmap
  heatmap(corrMat,
          scale = "none",
          col = pallette(1000))
  
  # print min-max corr vals
  print(range(corrMat))
  
  return(corrMat)
  
}


## Correlation matrix between an Parcel and every other grayordinate
# THIS CAN BE OPTIMIZED WITH A BETTER USE OF COR()
partoverCorrMat <- function(ROI = 'R_7m_ROI', lbels = verts[[1]]){
  
  # extract the time series from the ROI
  indexing <- which(rownames(ptSeries)==as.character(ROI))
  ROI_tseries <- as.numeric(ptSeries[indexing, ])
  
  # get the remaining number of vertices
  nVerts <- length(lbels) 
  
  # Create data frame to store values
  parcelIndx <- which(labelCoords_parcel$Label==ROI)
  corrVec <- data.frame(Seed = rep(ROI, nVerts),
                        Labels = lbels,
                        xstart = rep(labelCoords_parcel$x[parcelIndx], nVerts), # 7.33324 for R_7m
                        ystart = rep(labelCoords_parcel$y[parcelIndx], nVerts), # -63.142
                        zstart = rep(labelCoords_parcel$z[parcelIndx], nVerts)) # 42.5578
  
  for(vertex in seq(nVerts)){
    
    # temp series
    tempVertex <- timeSeries[vertex, ] 
    
    # store the Pearson correlation in the corr matrix
    #corrVec$Correlation[vertex] <- cor(ROI_tseries, tempVertex, method = "pearson") 
    
    # Another option, so that FDR can be applied
    tempCor <- cor.test(ROI_tseries, tempVertex, method = "pearson")
    corrVec$Correlation[vertex] <- tempCor$estimate
    corrVec$pval[vertex] <- tempCor$p.value
    
  }
  
  # Turn ROI vertices into 0s
  indexing <- which(corrVec$Labels==as.character(ROI))
  corrVec$Correlation[indexing] <- 0
  
  # Create a column with p-vals corrected for multiple comparisons
  corrVec$adjPval <- p.adjust(corrVec$pval, "BY")
  
  # normalize correlation vals for analysis (see below)
  tempTanh <- fisherTanh(Data = corrVec$Correlation)
  corrVec$tanhZ <- tempTanh$tanhZ
  corrVec$tanhPvals <- tempTanh$pvals
  corrVec$tanhPAdjusted <- tempTanh$adjustPvals
  
  return(corrVec)
  
}


## tanh-z transformation (variance stabilizing Fisher) and p-values (adjusted and not)
# This takes either a matrix of correlation values (vectors too, but manually compute pvals)
# Normalization approach suggested in network textbook (equation 7.10)
fisherTanh <- function(Data = padjMatrix){
  
  transformed <- list()
  
  # tanh
  Data[Data == 1] <- 1 * 0.999
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
  }
  
  return(transformed)
  
}


## Extract the centroid-most vertex from each parcel
parcelCentroid <- function(ROI = 'R_7m_ROI', all_coordinates = labelCoords_vertex){
  
  # The coordinate file must contain a Label column
  
  # Get the ROI-specific vertices
  indx <- grep(ROI, all_coordinates$Label)
  
  # Grab the coordinates
  tempCoords <- all_coordinates[indx, c("x","y","z")]
  rownames(tempCoords) <- indx
  
  # Compute the distance among all vertices, and turn into a matrix
  tempDist <- dist(tempCoords, diag = T)
  tempDist <- as.matrix(tempDist)
  
  # Compute the sum of distances for each vertex, and get the minimum
  sumDists <- colSums(tempDist)
  minDist <- which(sumDists == min(sumDists))
  minDist <- indx[minDist]
  
  return(minDist)
  
}


## This function takes the output from partoverCorrMat() and adds/removes the characteristics we decided on.
prepROItoVer <- function(Data = R_7m_allCorr, Coordinates = labelCoords_vertex){
  
  # Add vertex coordinates
  tempDframe <- cbind(Data, Coordinates[,1:3])
  
  # Remove rows with non-significant adjusted pvalues
  indx <- tempDframe$adjPval < 0.05
  tempDframe <- tempDframe[indx, ]
  
  # Characterize correlation direction and round up vals for plotting
  tempDframe$Relation <- ifelse(tempDframe$Correlation < 0, "Negative", "Positive")
  tempDframe$Correlation <- round(tempDframe$Correlation, digits = 2)
  
  # Add a column to differentiate between left/right hemispheres
  tempDframe$Hemisphere <- substring(tempDframe$Labels,1,1)
  
  return(tempDframe)
  
}


## Plot correlation from ROI to rest of vertices
plotROItoVertex <- function(Data = R_7m_allCorr, ROI = 'R_7m_ROI', ColRange = Cols, View = "Axial", Legends = TRUE){
  
  ROIindx <- grep(ROI, Data$Labels)
  ROIvertices <- Data[ROIindx,c("x","y","z")] # grab ROI-specific vertices to black out
  labelCoord <- ROIvertices[1,] # just to place the label
  
  if (View == "Axial") {
    ggplot() +
      geom_point(data = Data, aes(x=x, y=y, alpha=.5, color = Correlation), show.legend = Legends) +
      geom_point(data = labelCoords_vertex, aes(x=x, y=y), alpha=0.01) +
      geom_nodes(data = ROIvertices, aes(x=x, y=y)) +
      geom_nodelabel_repel(aes(x=labelCoord$x, y=labelCoord$y, label = sub("_ROI","",ROI))) +
      scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
      theme_blank()
  } else if (View == "Sagittal") {
    ggplot() +
      geom_point(data = Data, aes(x=y, y=z, alpha=.5, color = Correlation), show.legend = Legends) +
      geom_point(data = labelCoords_vertex, aes(x=y, y=z), alpha=0.01) +
      geom_nodes(data = ROIvertices, aes(x=y, y=z)) +
      geom_nodelabel_repel(aes(x=labelCoord$y, y=labelCoord$z, label = sub("_ROI","",ROI))) +
      scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
      theme_blank()
  } else if (View == "Coronal") {
    ggplot() +
      geom_point(data = Data, aes(x=x, y=z, alpha=.5, color = Correlation), show.legend = Legends) +
      geom_point(data = labelCoords_vertex, aes(x=x, y=z),  alpha=0.01) +
      geom_nodes(data = ROIvertices, aes(x=x, y=z), alpha=.5) +
      geom_nodelabel_repel(aes(x=labelCoord$x, y=labelCoord$z, label = sub("_ROI","",ROI))) +
      scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
      theme_blank()
  } else if (View == "Medial Right") {
    if (substring(ROI,1,1) == "R") {
      oneHemi <- Data[grep("R", Data$Hemisphere), ]
      oneHemi$Medial <- ifelse(oneHemi$x < 20, "Medial", "Other")
      medialHemi <- oneHemi[which(oneHemi$Medial=="Medial"), ]
      medialHemi$y <- medialHemi$y * -1
      ROIvertices$y <- ROIvertices$y * -1
      labelCoord$y <- labelCoord$y * -1
      oneHemi_all <- labelCoords_vertex[grep("R",labelCoords_vertex$Hemisphere), ] # whole surface
      oneHemi_all$Medial <- ifelse(oneHemi_all$x < 20, "Medial", "Other")
      medialHemi_all <- oneHemi_all[grep("Medial", oneHemi_all$Medial), ]
      medialHemi_all$y <- medialHemi_all$y * -1
      ggplot() +
        geom_point(data = medialHemi, aes(x=y, y=z, alpha=.1, color = Correlation), show.legend = Legends) +
        #geom_point(data = medialHemi_all, aes(x=y, y=z), alpha=0.07) +
        scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
        geom_nodes(data = ROIvertices, aes(x=y, y=z), alpha=.5) +
        geom_nodelabel_repel(aes(x=labelCoord$y, y=labelCoord$z, label = sub("_ROI","",ROI))) +
        theme_blank()
    } else {
      oneHemi <- Data[grep("R", Data$Hemisphere), ]
      oneHemi$Medial <- ifelse(oneHemi$x < 20, "Medial", "Other")
      medialHemi <- oneHemi[which(oneHemi$Medial=="Medial"), ]
      medialHemi$y <- medialHemi$y * -1
      ROIvertices$y <- ROIvertices$y * -1
      labelCoord$y <- labelCoord$y * -1
      oneHemi_all <- labelCoords_vertex[grep("R",labelCoords_vertex$Hemisphere), ] # whole surface
      oneHemi_all$Medial <- ifelse(oneHemi_all$x < 20, "Medial", "Other")
      medialHemi_all <- oneHemi_all[grep("Medial", oneHemi_all$Medial), ]
      medialHemi_all$y <- medialHemi_all$y * -1
      ggplot() +
        geom_point(data = medialHemi, aes(x=y, y=z, alpha=.1, color = Correlation), show.legend = Legends) +
        #geom_point(data = medialHemi_all, aes(x=y, y=z), alpha=0.07) +
        scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
        theme_blank()
    }
  } else if (View == "Medial Left") {
    if (substring(ROI,1,1) == "L") {
      oneHemi <- Data[grep("L", Data$Hemisphere), ]
      oneHemi$Medial <- ifelse(oneHemi$x > -20, "Medial", "Other")
      medialHemi <- oneHemi[which(oneHemi$Medial=="Medial"), ]
      oneHemi_all <- labelCoords_vertex[grep("R",labelCoords_vertex$Hemisphere), ] # whole surface
      oneHemi_all$Medial <- ifelse(oneHemi_all$x < 20, "Medial", "Other")
      medialHemi_all <- oneHemi_all[grep("Medial", oneHemi_all$Medial), ]
      ggplot() +
        geom_point(data = medialHemi, aes(x=y, y=z, alpha=.1, color = Correlation), show.legend = Legends) +
        #geom_point(data = medialHemi_all, aes(x=y, y=z), alpha=0.07) +
        scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
        geom_nodes(data = ROIvertices, aes(x=y, y=z), alpha=.5) +
        geom_nodelabel_repel(aes(x=labelCoord$y, y=labelCoord$z, label = sub("_ROI","",ROI))) +
        theme_blank()
    } else {
      oneHemi <- Data[grep("L", Data$Hemisphere), ]
      oneHemi$Medial <- ifelse(oneHemi$x > -20, "Medial", "Other")
      medialHemi <- oneHemi[which(oneHemi$Medial=="Medial"), ]
      oneHemi_all <- labelCoords_vertex[grep("R",labelCoords_vertex$Hemisphere), ] # whole surface
      oneHemi_all$Medial <- ifelse(oneHemi_all$x < 20, "Medial", "Other")
      medialHemi_all <- oneHemi_all[grep("Medial", oneHemi_all$Medial), ]
      ggplot() +
        geom_point(data = medialHemi, aes(x=y, y=z, alpha=.1, color = Correlation), show.legend = Legends) +
        #geom_point(data = medialHemi_all, aes(x=y, y=z), alpha=0.07) +
        scale_color_gradient2(low = ColRange[1], mid = "white", high = ColRange[2], limits = c(-1,1)) +
        theme_blank()
    }
  }
  
}


## Run a fastgreedy modularity community detection on ROIs
# This function relies on having the timeSeries data uploaded, and labelCoords_vertex 
# Extras dictates whether the community object + correlation matrix should also be extracted
communityDetection <- function(Data = parcelBins$First, ROIS = "None", Type = "vertex", thresh = F, extras = F) {
  
  print(paste('Computing modularity based on', Type))
  
  if (ROIS == "None"){
    
    # This will just do SP for now, for sliding window
    print("Previously concatenated data")
    
    corrMat <- cor(t(Data))
    corrMatrix <- corrMat
    transfMat <- fisherTanh(Data = corrMat)
    if (thresh == T) {
      transfMat$tanhZ[transfMat$adjustPvals > 0.05] <- 0
    }
    corrMat <- transfMat$tanhZ
    
    diag(transfMat$tanhZ) <- 0    
    diag(corrMat) <- 0
    
    corrMat <- exp(corrMat)
    corrMat[corrMat==1] <- 0
    transfMat$tanhZ <- corrMat
    
    tempGraph <- graph_from_adjacency_matrix(corrMat, weighted = T, mode = "undirected")
    tempLap <- laplacian_matrix(tempGraph, normalized=T)
    tempEigen <- eigen(tempLap)
    f.vec <- length(tempEigen$values) - 1
    tempEigen$binarized <-  as.factor(ifelse(tempEigen$vectors[,f.vec] > 0, 1, 0)) # binarized Fiedler Vector 
    
    summary <- data.frame(Label = colnames(corrMat),
                          Hemisphere = substring(colnames(corrMat),1,1),
                          EigenVal = tempEigen$values,
                          FiedlerVec = tempEigen$vectors[, (length(tempEigen$values) - 1)],
                          FiedlerBinary = tempEigen$binarized)
    
    ## Get the final components
    if (extras == T) {
      modularityResults <- list(CorrMatrix = corrMatrix,
                                TransfMatrix = transfMat,
                                Summary = summary)
    } else {
      modularityResults <- list(Summary = summary)
    }
    
  } else {
    
    if (Type == "vertex") {
      # To store the vertex indices corresponding to the ROIs
      indx <- numeric()
      
      # had to place dashes on each side because grep grabbed strings containing the names (i.e. 47m, a24pr)
      #ROIS <- c("_7m_", vmPFC_labels) # c("_a24_", "_7m_")
      for (ROI in ROIS) {
        indx <- c(indx, grep(ROI, rownames(Data)))  
      }
      
      nVerts <- length(indx)
      ROI_tseries <- Data[indx, ]
      
      # This used to be done with the for loop, but it was too slow. cor() speeds up the process by a lot
      corrMat <- cor(t(ROI_tseries))
      
      # Check for outliers. 
      # If found, remove vertices from indx, coords, and data, and recompute
      # tempStr <- colSums(corrMat)
      # outliers <- boxplot.stats(tempStr)$out
      # if (length(outliers) > 0) {
      #   outIndx <- which(tempStr %in% outliers)
      #   indx <- indx[! outIndx]
      #   ROI_tseries <- Data[indx, ]
      #   corrMat <- cor(t(ROI_tseries))
      # }
      
      # name the dimensions of the matrix according to the surface vertex index
      rownames(corrMat) <- indx
      
      # for storing later
      corrMatrix <- corrMat
      
      
      # transform to Fisher's (think of thresholding)
      transfMat <- fisherTanh(Data = corrMat)
      
      # Determine if edges should be thresholded or not
      if (thresh == T) {
        transfMat$tanhZ[transfMat$adjustPvals > 0.05] <- 0
      }
      
      # Store Fisher transformed vals for graphing
      corrMat <- transfMat$tanhZ
      
      # diagonals of 1 could be interpreted as self-loops
      diag(transfMat$tanhZ) <- 0    
      diag(corrMat) <- 0
      
      # Exponentiate to preserve distribution while ensuring positive weights
      # I'm keeping corrMat and transfMat$tanhZ separate in case I want to uncorrect transfMat in the future
      corrMat <- exp(corrMat)
      corrMat[corrMat==1] <- 0
      transfMat$tanhZ <- corrMat
      
      # community detection
      # I initially used the absolute value of the correlation, but the exp preserves the distribution 
      # Next, try using the fisher transform
      
      tempGraph <- graph_from_adjacency_matrix(corrMat, weighted = T, mode = "undirected")
      tempCommunity <- fastgreedy.community(tempGraph)
      
      # community object
      #modularityResults$community <- tempCommunity
      
      # correlation matrix with transformed values
      #modularityResults$corrMat <- corrMatrix
      
      # get coordinate info from selected regions (useful for ggplot)
      summary <- data.frame(Vertex = indx,
                            Label = tempCommunity$names,
                            Membership = tempCommunity$membership,
                            Modularity = tempCommunity$modularity,
                            x = labelCoords_vertex[indx, "x"],
                            y = labelCoords_vertex[indx, "y"],
                            z = labelCoords_vertex[indx, "z"],
                            Hemisphere = substring(tempCommunity$names,1,1))
      
      
    } else if (Type == "parcels") {
      
      # This used to be done with the for loop, but it was too slow. cor() speeds up the process by a lot
      corrMat <- cor(t(Data))
      
      # for storing later
      corrMatrix <- corrMat
      
      # transform to Fisher's (think of thresholding)
      transfMat <- fisherTanh(Data = corrMat)
      
      # Determine if edges should be thresholded or not
      if (thresh == T) {
        transfMat$tanhZ[transfMat$adjustPvals > 0.05] <- 0
      }
      
      # Store Fisher transformed vals for graphing
      corrMat <- transfMat$tanhZ
      
      # diagonals of 1 could be interpreted as self-loops
      diag(transfMat$tanhZ) <- 0
      diag(corrMat) <- 0
      
      # Exponentiate to preserve distribution while ensuring positive weights
      # I'm keeping corrMat and transfMat$tanhZ separate in case I want to uncorrect transfMat in the future
      corrMat <- exp(corrMat)
      corrMat[corrMat==1] <- 0
      transfMat$tanhZ <- corrMat
      
      # community detection
      # I initially used the absolute value of the correlation, but the exp preserves the distribution 
      # Next, try thresholding it by the adjusted p-vals
      tempGraph <- graph_from_adjacency_matrix(corrMat, weighted = T, mode = "undirected")
      tempCommunity <- fastgreedy.community(tempGraph)
      
      # community object
      #modularityResults$community <- tempCommunity
      
      # correlation matrix with transformed values
      #modularityResults$corrMat <- corrMatrix
      
      # get coordinate info from selected regions (useful for ggplot)
      summary <- data.frame(Label = tempCommunity$names,
                            Membership = tempCommunity$membership,
                            Modularity = tempCommunity$modularity,
                            x = labelCoords_parcel[ ,"x"],
                            y = labelCoords_parcel[ ,"y"],
                            z = labelCoords_parcel[ ,"z"],
                            Hemisphere = substring(tempCommunity$names,1,1))
      
    }
    
    
    ## Get the final components
    modularityResults <- list(Community = tempCommunity,
                              CorrMatrix = corrMatrix,
                              TransfMatrix = transfMat,
                              Summary = summary)
    
  }
  
  return(modularityResults)
  
}


## plot the communities from communityDetection
plotCommunities <- function(Data = modularityResults, Hemi = "R", type = "Membership", Legends = T, bground=0.1, Cols = c("aquamarine4", "#D9541A")) {
  # This takes the output from the communityDetection function and plots them on the medial wall (since we're interested in medial now)
  # type is the name of the column to plot
  if (Hemi == "R") {
    # "Medial Right"
    oneHemi <- Data[grep("R", Data$Hemisphere), ] # from the ROIs
    oneHemi$Medial <- ifelse(oneHemi$x < 20, "Medial", "Other")
    medialHemi <- oneHemi[grep("Medial", oneHemi$Medial), ]
    medialHemi$y <- medialHemi$y * -1
    oneHemi_all <- labelCoords_vertex[grep("R",labelCoords_vertex$Hemisphere), ] # whole surface
    oneHemi_all$Medial <- ifelse(oneHemi_all$x < 20, "Medial", "Other")
    medialHemi_all <- oneHemi_all[grep("Medial", oneHemi_all$Medial), ]
    medialHemi_all$y <- medialHemi_all$y * -1
    ggplot() + 
      geom_point(data = medialHemi, aes_string(x="y", y="z", color = type), show.legend=Legends) +
      geom_point(data = medialHemi_all, aes(x=y, y=z), alpha=bground) +
      scale_color_gradient(low = Cols[1], high = Cols[2]) +
      theme_blank() 
  } else if (Hemi == "L") {
    # "Medial Left"
    oneHemi <- Data[grep("L", Data$Hemisphere), ]
    oneHemi$Medial <- ifelse(oneHemi$x > -20, "Medial", "Other")
    medialHemi <- oneHemi[which(oneHemi$Medial=="Medial"), ]
    oneHemi_all <- labelCoords_vertex[grep("L",labelCoords_vertex$Hemisphere), ] # whole surface
    oneHemi_all$Medial <- ifelse(oneHemi_all$x > -20, "Medial", "Other")
    medialHemi_all <- oneHemi_all[grep("Medial", oneHemi_all$Medial), ]
    ggplot() + 
      geom_point(data = medialHemi, aes_string(x="y", y="z", color = type), show.legend=Legends) +
      geom_point(data = medialHemi_all, aes(x=y, y=z), alpha=bground) +
      scale_color_gradient(low = Cols[1], high = Cols[2]) +
      theme_blank() 
  }
}


## get ROI coords & index
# The point here is to reduce the summary dframes from community detection to show only ROIs
# Should work for extracting any label-indexed dframe though
# I wanted to also get the index in case I want to extract specific rows from parcel/vertex coord dframes
getCoords <- function(Labels = DMN_labels, Coords = labelCoords_parcel, TimeSeries = FALSE){
  
  indx <- numeric()
  
  # If you want to select time series from raw data
  if (TimeSeries == TRUE) {
    
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


# eigen value community detection
eigenVals <- function(Data) {
  
  # This will produce a vector of values, ready to be plotted
  tempGraph <- graph_from_adjacency_matrix(Data$TransfMatrix$tanhZ, mode = "undirected", weighted = T)
  tempGraph <- laplacian_matrix(tempGraph, normalized=T)
  tempGraph <- eigen(tempGraph)
  f.vec <- length(tempGraph$values) - 1
  tempGraph$binarized <-  as.factor(ifelse(tempGraph$vectors[,f.vec] > 0, 1, 0)) # binarized Fiedler Vector
  
  return(tempGraph)
  
}


# Attempt at setting up data for confusion matrices and Jaccard index calculations
confusionMatrix <- function(partition1 = parcelCommunities[[1]]$Membership, partition2 = parcelCommunities[[1]]$FiedlerBinary) {
  
  # Divide partitions to evaluate
  # Usually 2 will be fiedler
  
  # Check if any partition has 0s (since I binarize the Fiedler vector)
  if (0 %in% partition1) {
    
    partition1[grep(0, partition1)] <- 1
    partition1[grep(1, partition1)] <- 2
    
  } 
  
  if (0 %in% partition2) {
    
    partition2[grep(1, partition2)] <- 2   
    partition2[grep(0, partition2)] <- 1
    
  }
  
  # Community sizes
  commSizes1 <- table(partition1)
  commSizes2 <- table(partition2)
  
  # Number of communities per partition
  nComms1 <- length(commSizes1)
  nComms2 <- length(commSizes2)
  
  # Number of vertices
  n <- length(partition1)
  
  # Putting together elements of the confusion matrix
  confMatrix <- matrix(nrow = nComms1,
                       ncol = nComms2)
  
  for (i in as.numeric(case.names(commSizes1))) {
    
    # Vertices belonging to community qX of partition X
    tempComm1 <- partition1 == i
    
    for (j in as.numeric(case.names(commSizes2))) {
      
      # Vertices belonging to community qY of partition Y
      tempComm2 <- partition2 == j
      
      # Populate matrix
      confMatrix[i,j] <- sum(tempComm1 & tempComm2) 
      
    }
  }
  
  if (sum(confMatrix) != n) {warning('Sum of the confusion matrix is not equal to number of vertices')}
  
  
  return(confMatrix)
  
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


# number of times a node changes affiliation in a time series
flexibility <- function(Data = cbind(1:10, 6:15)) {
  # Calculation of flexibility based on Garcia et al., 2018  
  # This function takes in a data frame or matrix in which columns are the community affiliation 
  # It needs at least 2 colums, although it's pointless for that  
  
  # How many jumps can there be?
  nJumps <- dim(Data)[2] - 1
  
  # Matrix to store jumps
  jumpCount <- matrix(nrow = dim(Data)[1], ncol = nJumps)
  
  # For each transition, get the which nodes changed affiliation
  for (jump in seq(nJumps)) {
    
    jumpCount[, jump] <- Data[, jump] != Data[, jump + 1]
    
  }
  
  # How many times did every node jump?
  totalJumps <- rowSums(jumpCount)
  
  # Calculate the flexibility per node
  flexibility <- totalJumps / nJumps
  
  return(flexibility)
  
}


# Sliding window
slidingWindow <- function(subjTS = ROI_timeSeries[[1]], mins = 15, jump = 1, Spectral = T, Modularity = F, ROIs = c(vmPFC_labels)) {
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
  TS <- dim(subjTS)[2] # time series for the subject
  WS <- seq(834 * (mins/10)) # window from the first TR up to 834 (~10 mins) times the desired multiplier
  jump <- 84 * jump # 84 ~ 1 min, times the number of mins that the window moves
  nJumps <- floor((TS - length(WS)) / jump) # number of jumps to be performed, based on the selected parameters
  
  ##------- using lapply
  winData <- mclapply(seq(nJumps), function(x) subjTS[,WS+(jump*(x-1))])
  commTS <- mclapply(winData, communityDetection, ROIS = "None", Type = "vertex", thresh = T, extras = F)
  return(commTS)
  
}


# compare sliding window data to overall communities
slideCompare <- function(subjData = slideCommunities[[1]], template = vmpfc7mCommunities[[1]], func = "RI", comm = "Spectral") {
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
  
  # For storing RIs or VIs
  tempComparison <- numeric()
  
  # How many jumps does the original data contain?
  nJumps <- length(subjData)
  
  # This loop is technically backwards. I should technically divide by function, then partition method, and then run the window comparisons
  # It's still really fast, so I won't worry.
  for (Win in seq(nJumps)) {
    if (comm == "Modularity") {
      if (func == "RI") {
        tempComparison[Win] <- arandi(subjData[[Win]]$Membership, template$Membership, adjust = T)
      } else if (func == "VI") {
        tempComparison[Win] <- vi.dist(subjData[[Win]]$Membership, template$Membership)
      }
    } else if (comm == "Spectral") {
      if (func == "RI") {
        tempComparison[Win] <- arandi(subjData[[Win]]$FiedlerBinary, template$FiedlerBinary, adjust = T)
      } else if (func == "VI") {
        tempComparison[Win] <- vi.dist(subjData[[Win]]$FiedlerBinary, template$FiedlerBinary)
      }
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
# The input should be the summary from community partitioning
# Once this is created, go to the terminal and input something like this
# wb_command -cifti-convert -from-text dataforCifti.txt 100307.MyelinMap_BC.32k_fs_LR.dscalar.nii testCifti.dscalar.nii
# Where the myelin file here is just a templace. It can be any dscalar.nii with the right surface size
HCPOut <- function(Data = dmnval7mCommunities[[1]], MOI = "Membership", SubjID = "100307"){
  
  nVertices <- 59412
  tempVec <- rep(-1, nVertices)
  temp <- grep(MOI, colnames(Data))
  tempVec[Data$Vertex] <- Data[[temp]]
  write.table(file = paste(SubjID,"_",MOI,'_dataforCifti.txt', sep=""), tempVec, row.names = F, col.names = F, dec = ".")
  
}

# Perform pairwise comparisons of clustering outcomes on all subjects
comparePartitions <- function(Data = dmnval7mCommunities, MOI = "FiedlerBinary", Index = "VI", nSubjects = nSubj, subjNames = subjList) {
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
          indexMatrix[subj,subj2] <- vi.dist(allVecs[, subj], allVecs[, subj2])
        }
      }
    } else if (Index == "RI") {
      for (subj in seq(nSubjects)) {
        for (subj2 in seq(nSubjects)) {
          indexMatrix[subj,subj2] <- arandi(allVecs[, subj], allVecs[, subj2], adjust = T)
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
Parcellated <- F # Parcellated whole brain community detection
Vertex <- T # ROI-based vertex analysis
Plots <- F # If running on the SCC, plots might be unhelpful. Have a separate script for that instead
Sliding <- T

# colors to differentiate other things
#Cols <- c("aquamarine4","#D9541A",rgb(190,190,190,100, maxColorValue = 255)) # left, right, interhemisphere
Cols <- c("DMN" = "aquamarine4","Valuation" = "#D9541A",rgb(190,190,190,100, maxColorValue = 255)) # left, right, interhemisphere

# Label groups for specific analyses (DMN and DEC)
# These come from an arbitrary division based on where Mackey & Petrides (2014) draw their line
# This is computed below, but if we want to save time it can be loaded instead.
#bothLbls <- as.character(read.table('bothLbls.csv', header = F)$V2)

vmPFC_labels <- c("_10pp_",
                  "_10r_",
                  "_10v_",
                  "_25_",
                  "_OFC_",
                  "_a24_",
                  "_p32_",
                  "_pOFC_",
                  "_s32_")

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

Yeo_PCC_labels <- c("_7m_",
                    "_23d_",
                    "_31a_",
                    "_31pd_",
                    "_31pv_",
                    "_PCV_",
                    "_POS1_",
                    "_RSC_",
                    "_d23ab_",
                    "_v23ab_")

# mPFC_only <- c("_a24_",
#                "_d32_",
#                "_p32_",
#                "_10r_",
#                "_9m_",
#                "_10v_",
#                "_OFC_",
#                "_s32_",
#                "_p24_")

mPFC_only <- c("_a24_",
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
if (FALSE %in% (mPFC_only %in% c(mPFC_only, "Blah"))) {
  stop("mPFC-only and Overlapping ROIs don't match")
}

# Cortical surface for plotting
labelPerVertex <- read.table('labelPerVertex.csv', header = F)
labelCoords_vertex <- read.csv2('labelCoords_vertex.csv', sep = ",")[,2:6]
labelCoords_parcel <- read.csv2('labelCoords_parcel.csv', sep = ",")

# Ensure that coordinates are numeric once loaded
labelCoords_vertex <- transform(labelCoords_vertex, x = as.numeric(as.character(x)))
labelCoords_vertex <- transform(labelCoords_vertex, y = as.numeric(as.character(y)))
labelCoords_vertex <- transform(labelCoords_vertex, z = as.numeric(as.character(z)))

# Data frame to store the graph densities of each network
densities <- data.frame(SubjID = rep(SubjID, 7),
                        Type = c("Overall", "Day 1", "Day 2", "Session 1", "Session 2", "Session 3", "Session 4"),
                        Density = rep(NA, 7))



###------- Structural Data ---------------
write('Loading structural brain data from HCP subjects...', stdout())

#setwd('../../../restricted/projectnb/cd-lab/Claudio/Community/')
# Myelin density maps (all HCP subjects)
myelinData <- read.csv(paste(SubjID,'_myelin.txt', sep=""), header = F)
myelinData$Label <- labelCoords_vertex$Label

# Curvature (update to all HCP subjects)
curvatureData <- read.csv(paste(SubjID,'_curvature.txt', sep=""), header = F)
curvatureData$Label <- labelCoords_vertex$Label

# Cortical thickness, with curvature regressed out (update to all HCP subjects)
thicknessData <- read.csv(paste(SubjID,'_thickness.txt', sep=""), header = F)
thicknessData$Label <- labelCoords_vertex$Label

# Get myelin density for ROIs
myelinROI <- getCoords(Labels = bothLbls, Coords = myelinData, TimeSeries = F)$Coords

# Get curvature for ROIs
curvatureROI <- getCoords(Labels = bothLbls, Coords = curvatureData, TimeSeries = F)$Coords

# Get thickness for ROIs
thicknessROI <- getCoords(Labels = bothLbls, Coords = thicknessData, TimeSeries = F)$Coords

###------- Load HCP Data ---------------
# Prep subject data from HCP
# Load everysubjects' time series
if (Parcellated == T) {
  
  write("Loading parcellated time series...", stdout())
  
  temp <- list.files(pattern="*_ptSeries.txt")
  ptSeries <- mclapply(temp, fread, sep="\t", header=F)
  ptSeries <- mclapply(ptSeries, data.matrix)
  
  # Populate each subject's column and row names
  for (i in seq(length(ptSeries))) {
    nTPoints <- seq(dim(ptSeries[[i]])[2])
    dimnames(ptSeries[[i]]) <- list(labelCoords_parcel$Label, nTPoints)
  }
  
}

if (Vertex == T) {
  
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
  
}

rm(temp)



###------- Community Detection: all ---------------
## Parcellated
if (Parcellated == T) {
  
  write("Computing community detection for the parcellated brain...", stdout())
  
  # Modularity
  parcelCommunities <- mclapply(ptSeries, communityDetection, ROIS = bothLbls, Type = "parcels", thresh = T)
  
  # Spectral
  tempEigen <- mclapply(parcelCommunities, eigenVals)
  
  # Store SP attributes and compute descriptives
  for (i in seq(length(tempEigen))) {
    
    # SP
    parcelCommunities[[i]]$Summary$EigenVal <- tempEigen[[i]]$values
    parcelCommunities[[i]]$Summary$FiedlerVec <- tempEigen[[i]]$vectors[, (length(tempEigen[[i]]$values) - 1)] # grab the Fiedler vector
    parcelCommunities[[i]]$Summary$FiedlerBinary <- tempEigen[[i]]$binarized
    
    # Descriptives
    tempGraph <- graph_from_adjacency_matrix(parcelCommunities[[i]]$TransfMatrix$tanhZ,
                                             mode = "undirected",
                                             weighted = T)
    parcelCommunities[[i]]$Summary$Strength <- strength(tempGraph)
    parcelCommunities[[i]]$Summary$Betweenness <- betweenness(tempGraph)
    
  }
  
  rm(tempStrength, tempGraph)
  
}



## Vertex-wise
if (Vertex) {
  
  write("Computing community detection for the ROI vertices...", stdout())
  
  # Modularity
  dmnval7mCommunities <- communityDetection(Data = timeSeries, ROIS = bothLbls, Type = "vertex", thresh = T)
  
  # Spectral
  tempEigen <- eigenVals(dmnval7mCommunities)
  
  # Store SP attributes and compute descriptives
  # SP
  dmnval7mCommunities$Summary$EigenVal <- tempEigen$values
  dmnval7mCommunities$Summary$FiedlerVec <- tempEigen$vectors[, (length(tempEigen$values) - 1)] # grab the Fiedler vector
  dmnval7mCommunities$Summary$FiedlerBinary <- tempEigen$binarized
  
  # Descriptives
  tempGraph <- graph_from_adjacency_matrix(dmnval7mCommunities$TransfMatrix$tanhZ,
                                                mode = "undirected",
                                                weighted = T)
  
  # Store the graph density
  densities$Density[grep("Overall", densities$Type)] <- graph.density(tempGraph)
  
  # Descriptives (Betweenness takes a long time, disregard)
  dmnval7mCommunities$Summary$Strength <- strength(tempGraph)
  
  # Append the coefficient of variation
  # dmnval7mCommunities$Summary$varCoeff <- varCoeff
  
  # Retain summary only from community detection (this is for now, since the matrices might be useful in the future)
  dmnval7mCommunities <- dmnval7mCommunities$Summary
  
  # Make sure the coordinates are numeric
  dmnval7mCommunities <-transform(dmnval7mCommunities, x = as.numeric(as.character(x)), y = as.numeric(as.character(y)), z = as.numeric(as.character(z)))
  
  # Make sure the binarized Fiedler Vector is a factor
  dmnval7mCommunities <- transform(dmnval7mCommunities, FiedlerBinary = as.factor(FiedlerBinary))
  
  # And ensure that DMN is assciated with the positive values of the eigenvector
  dmnval7mCommunities <- evenSpectral(dmnval7mCommunities)
  
  # Myelin
  dmnval7mCommunities$Myelin <- myelinROI[ , 1]
  
  # Curvature
  dmnval7mCommunities$Curvature <- curvatureROI[ , 1]
  
  # Thickness
  dmnval7mCommunities$Thickness <- thicknessROI[ , 1]
  
  rm(tempGraph)
  
}



###------- Measures of variance -----------
if (Vertex) {
  
  write("Computing each vertex's variance...", stdout())
  
  ## Raw mPFC and PCC SNR
  # Reduce the 4 SNR-dedicated tseries to the ROIs
  timeSeries_SNR_ROI <- lapply(timeSeries_SNR, getCoords, Labels = bothLbls, TimeSeries = T)
  rm(timeSeries_SNR)
  
  # Compute coefficient of variation
  SNR_mean <- lapply(timeSeries_SNR_ROI, function(data) rowMeans(data$Coords))
  SNR_sd <- lapply(timeSeries_SNR_ROI, function(data) apply(data$Coords, 1, sd))
  varCoeff <- do.call(cbind, lapply(seq_along(length(SNR_mean)), function(x) SNR_mean[[x]] / SNR_sd[[x]]))
  
  # Get the mean coefficient across sessions, and append it to the overall summary
  dmnval7mCommunities$varCoeff <- rowMeans(varCoeff)
  
  rm(SNR_mean, SNR_sd, timeSeries_SNR_ROI)
  
  ## Now compute the variance explained by the mean time series from the DMN community on each DMN community vertex
  ## This should be high, since the communities are correlation-based
  
  # Get the vertices from the DMN community
  # tempComm <- dmnval7mCommunities[dmnval7mCommunities$y > 0, ]
  index <- dmnval7mCommunities$Vertex[dmnval7mCommunities$FiedlerBinary == 1]
  tempY <- dmnval7mCommunities$y[dmnval7mCommunities$FiedlerBinary == 1]
  tempZ <- dmnval7mCommunities$z[dmnval7mCommunities$FiedlerBinary == 1]
  
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
  
  rm(tempComm, index, Rsquared, Tvalue, DMN_commTSeries, DMN_commTSeries_models, DMN_commTSeries_mean)
  
}
###------- Community Detection: Day 1 vs Day 2 -----------
if (Vertex) {
  
  write("Computing PCC-seed correlations and community detection for the ROI vertices, Day 1 vs Day 2...", stdout())
  
  # Correlate Yeo 2011 PCC ROI with mPFC
  # Get vertex indices
  Yeo_PCC_indx <- getCoords(Labels = Yeo_PCC_labels, Coords = timeSeries_halves[[1]], TimeSeries = T)$Index
  mPFC_indx <- getCoords(Labels = mPFC_only, Coords = timeSeries_halves[[1]], TimeSeries = T)$Index
  
  # Correlate the mean PCC time series with mPFC
  corHalves <- mclapply(timeSeries_halves, function(data) as.numeric(cor(t(data[mPFC_indx,]), colMeans(data[Yeo_PCC_indx,]))))
  
  # Fisher transform and FDR threshold (same as with the community detection)
  # corHalves_corrected <- lapply(corHalves, fisherTanh)
  # corHalves_corrected <- lapply(corHalves_corrected, function(data) ifelse(data$tanhZ == 0, 0, 1))
  # 
  # rm(corHalves)
  
  # Compute the day 1/2 communities too
  dmnval7mCommunities_halves <- lapply(timeSeries_halves, communityDetection, ROIS = bothLbls, Type = "vertex", thresh = T)
  
  # Spectral
  tempEigen_halves <- mclapply(dmnval7mCommunities_halves, eigenVals)
  
  # Store SP attributes and compute descriptives
  tempGraph <- list()
  for (i in seq(length(tempEigen_halves))) {
    
    # SP
    dmnval7mCommunities_halves[[i]]$Summary$EigenVal <- tempEigen_halves[[i]]$values
    dmnval7mCommunities_halves[[i]]$Summary$FiedlerVec <- tempEigen_halves[[i]]$vectors[, (length(tempEigen_halves[[i]]$values) - 1)] # grab the Fiedler vector
    dmnval7mCommunities_halves[[i]]$Summary$FiedlerBinary <- tempEigen_halves[[i]]$binarized
    
    # Descriptives
    # It takes too long to compute for the extended analyses.
    tempGraph[[i]] <- graph_from_adjacency_matrix(dmnval7mCommunities_halves[[i]]$TransfMatrix$tanhZ,
                                                  mode = "undirected",
                                                  weighted = T)
    
  }
  
  # Density
  densities$Density[grep("Day", densities$Type)] <- sapply(tempGraph, graph.density)
  
  # Descriptives (Betweenness takes a long time, disregard)
  # This code is horrible. Optimize with apply/do.call, etc.
  tempStrength <- mclapply(tempGraph, strength)
  
  for (i in seq(length(tempEigen_halves))) {
    
    dmnval7mCommunities_halves[[i]]$Summary$Strength <- tempStrength[[i]]
    
  }
  
  # Ensure that the binarized partition is consistent across windows (i.e. ~7m = 1)
  for (Win in seq(length(dmnval7mCommunities_halves))) {
    dmnval7mCommunities_halves[[Win]] <- dmnval7mCommunities_halves[[Win]]$Summary
    dmnval7mCommunities_halves[[Win]] <- transform(dmnval7mCommunities_halves[[Win]], x = as.numeric(as.character(x)), y = as.numeric(as.character(y)), z = as.numeric(as.character(z)))
    dmnval7mCommunities_halves[[Win]] <- transform(dmnval7mCommunities_halves[[Win]], FiedlerBinary = as.factor(FiedlerBinary))
    dmnval7mCommunities_halves[[Win]] <- evenSpectral(dmnval7mCommunities_halves[[Win]])
  }
  
  rm(tempStrength, tempGraph, timeSeries_halves)
  
}
###------- Community Detection: session ---------------
## Vertex-wise
if (Vertex == T) {
  
  write("Computing community detection for the ROI vertices, session...", stdout())
  
  # Modularity
  dmnval7mCommunities_sess <- lapply(timeSeries_sess, communityDetection, ROIS = bothLbls, Type = "vertex", thresh = T)
  
  # Spectral
  tempEigen_sess <- mclapply(dmnval7mCommunities_sess, eigenVals)
  
  # Store SP attributes and compute descriptives
  tempGraph <- list()
  for (i in seq(length(tempEigen_sess))) {
    
    # SP
    dmnval7mCommunities_sess[[i]]$Summary$EigenVal <- tempEigen_sess[[i]]$values
    dmnval7mCommunities_sess[[i]]$Summary$FiedlerVec <- tempEigen_sess[[i]]$vectors[, (length(tempEigen_sess[[i]]$values) - 1)] # grab the Fiedler vector
    dmnval7mCommunities_sess[[i]]$Summary$FiedlerBinary <- tempEigen_sess[[i]]$binarized
    
    # Descriptives
    # For MA703 I won't worry about this, since I'm not using it. 
    # It takes too long to compute for the extended analyses.
    tempGraph[[i]] <- graph_from_adjacency_matrix(dmnval7mCommunities_sess[[i]]$TransfMatrix$tanhZ,
                                                  mode = "undirected",
                                                  weighted = T)
    
  }
  
  # Density
  densities$Density[grep("Session", densities$Type)] <- sapply(tempGraph, graph.density)
  
  # Descriptives (Betweenness takes a long time, disregard)
  tempStrength <- mclapply(tempGraph, strength)
  
  for (i in seq(length(tempEigen_sess))) {
    
    dmnval7mCommunities_sess[[i]]$Summary$Strength <- tempStrength[[i]]
    
  }
  
  # Ensure that the binarized partition is consistent across windows (i.e. ~7m = 1)
  for (Win in seq(length(dmnval7mCommunities_sess))) {
    dmnval7mCommunities_sess[[Win]] <- dmnval7mCommunities_sess[[Win]]$Summary
    dmnval7mCommunities_sess[[Win]] <- transform(dmnval7mCommunities_sess[[Win]], x = as.numeric(as.character(x)), y = as.numeric(as.character(y)), z = as.numeric(as.character(z)))
    dmnval7mCommunities_sess[[Win]] <- transform(dmnval7mCommunities_sess[[Win]], FiedlerBinary = as.factor(FiedlerBinary))
    dmnval7mCommunities_sess[[Win]] <- evenSpectral(dmnval7mCommunities_sess[[Win]])
  }
  
  rm(tempStrength, tempGraph, timeSeries_sess)
  
}
###------- Parcel Communities: Descriptive summaries ---------------
# THINK ABOUT WRITING A CSV SUMMARY FOR THIS
if (Parcellated == T) {
  
  # Some cleaning up and averaging to plot the parcellated strength onto the brain
  tempStr <- cbind(parcelCommunities[[1]]$Strength,
                   parcelCommunities[[2]]$Strength,
                   parcelCommunities[[3]]$Strength)
  tempBtw <- cbind(parcelCommunities[[1]]$Betweenness,
                   parcelCommunities[[2]]$Betweenness,
                   parcelCommunities[[3]]$Betweenness)
  summaryParcels <- data.frame(Label = metaStuff$Label,
                               x = metaStuff$x,
                               y = metaStuff$y,
                               z = metaStuff$z,
                               Strength = rowMeans(tempStr),
                               Betweenness = rowMeans(tempBtw))
  
  # Append whether a label is DMN or not
  indx <- grep(paste(bothLbls, collapse="|"), parcelCommunities[[1]]$Label)
  indx2 <- seq(360)[-indx]
  parcelCommunities[[1]]$lblComb[indx] <-  "DMN/Val"
  parcelCommunities[[2]]$lblComb[indx] <-  "DMN/Val"
  parcelCommunities[[3]]$lblComb[indx] <-  "DMN/Val"
  parcelCommunities[[1]]$lblComb[indx2] <-  "Other"
  parcelCommunities[[2]]$lblComb[indx2] <-  "Other"
  parcelCommunities[[3]]$lblComb[indx2] <-  "Other"
  
  # Shortened versions
  summaryParcels$lblComb <- parcelCommunities[[1]]$lblComb
  summaryParcels_s <- getCoords(Labels = DMN_labels, Coords = summaryParcels)
  summaryParcels_comb <- getCoords(Labels = combinedLbls, Coords = summaryParcels)
  summaryParcels_both <- getCoords(Labels = bothLbls, Coords = summaryParcels)
  summaryParcels_dmn <- getCoords(Labels = dmnLbls, Coords = summaryParcels)
  
  # Short communities
  parcelCommunities_s <- list()
  parcelCommunities_s[[1]] <- getCoords(Labels = bothLbls, Coords = parcelCommunities[[1]])
  parcelCommunities_s[[2]] <- getCoords(Labels = bothLbls, Coords = parcelCommunities[[2]])
  parcelCommunities_s[[3]] <- getCoords(Labels = bothLbls, Coords = parcelCommunities[[3]])
  
}
###------- Parcel communities: Inter-subject agreement ---------------
# WRITE AS CSV, AND NOTE THAT THE VIs ARE UNADJUSTED (DIVIDE BY THE LOG OF THE NUMBER OF VERTICES)
if(Parcellated == T) {
  
  write("Computing within and between subject partition comparisons for parcellated brain...", stdout())
  
  # Inter subject SP
  # Rand index
  interRI_p <- comparePartitions(Data = parcelCommunities, Index = "RI")
  
  # Variation of information
  interVI_p <- comparePartitions(Data = parcelCommunities, Index = "VI")
  
  # Within subject cross-method
  # Rand index
  withinRI_p <- comparePartitions(Data = parcelCommunities, index = "RI", MOI = c("FiedlerBinary","Membership"))
  
  # Variation of information  
  withinVI_p <- comparePartitions(Data = parcelCommunities, index = "VI", MOI = c("FiedlerBinary","Membership"))
  
  save(interRI_p, interVI_p, withinRI_p, withinVI_p, file = "parcelComparisons.RData")
  
}
###------- Vertex communities: Sliding window ---------------
if (Vertex & Sliding) {
  
  write("Computing sliding window analysis for vertex communities...", stdout())
  
  # Prep data
  ROI_timeSeries <- getCoords(Coords = timeSeries, Labels = bothLbls, TimeSeries = T)
  
  # Sliding window analysis for each subject (using defaults)
  slideCommunities <- slidingWindow(subjTS = ROI_timeSeries, ROIs = bothLbls, Spectral = T, Modularity = F, mins = 20, jump = 1)
  
  # Ensure that the binarized partition is consistent across windows (i.e. ~7m = 1)
  for (Win in seq(length(slideCommunities))) {
    slideCommunities[[Win]] <- evenSpectral(Data = slideCommunities[[Win]]$Summary)
  }
  
  # Store the community structure per sliding window 
  slidingVals <- do.call(cbind, lapply(slideCommunities, "[[", "FiedlerBinary")) # do call will perform a function on a list. lapply creates a sublist from the original data
  slidingVals <- slidingVals - 1 # Lazy way of fixing the fact that the factorized FiedlerBinary becomes numeric as 1 and 2, thus making slideProp below messy (but fixable)

  # now compare to the template 
  # Get the index for each comparison
  tempCompare <- slideCompare(subjData = slideCommunities,
                              dmnval7mCommunities,
                              comm = "Spectral",
                              func = "RI")
  
  # Store as data frame 
  allSlideCompared <- data.frame(SubjID = rep(SubjID, length(tempCompare)),
                          Window = seq(length(tempCompare)),
                          Index = tempCompare) #scale(tempCompare, center = T, scale = F))
  

  # Add other relevant values to the final summary
  # With this info, we can potentially run a mixed effects GLM to estimate the probability of a vertex being mostly DMN (propDMN) based on myelin density
  # Flexibility
  testy <- slidingVals 
  dmnval7mCommunities$slideFlexibility <- flexibility(testy)
  
  # Proportion of times a node was associated with DNM
  dmnval7mCommunities$slidePropDMN <- rowMeans(testy)
  
  # Another one is to look at the fiedler vector and plot that. Perhaps nodes closer to 0 are also the most variable ones.
  testy <- do.call(cbind, lapply(slideCommunities, "[[", "FiedlerVec")) 
  dmnval7mCommunities$slideFiedVec_mean <- rowMeans(testy)
  
  # or the SD
  dmnval7mCommunities$slideFiedVec_sd <- transform(testy, SD=apply(testy,1, sd, na.rm = TRUE))$SD
    
}
###------- Save Data ---------------
if (Vertex) {
  
  write("Saving summary files...", stdout())
  
  # Produce a final summary for each participant
  write.csv(dmnval7mCommunities, paste(SubjID, "_finalSummary.csv", sep=""), row.names = F)
  
  # Final summary for each half
  lapply(seq_along(dmnval7mCommunities_halves), function(x) {write.csv(dmnval7mCommunities_halves[[x]], paste(SubjID, "_H", x,"_finalSummary.csv",sep=""), row.names = F)})
  lapply(seq_along(corHalves), function(x) {write.table(corHalves[[x]], paste(SubjID, "_H", x,"_corHalves.csv",sep=""), row.names = F, col.names = F)})
  
  # Final summary for each session
  lapply(seq_along(dmnval7mCommunities_sess), function(x) {write.csv(dmnval7mCommunities_sess[[x]], paste(SubjID, "_S", x,"_finalSummary.csv",sep=""), row.names = F)})
  
  # The variance explained by the mean DMN time series on each vertex
  write.csv(DMN_commVrtxvsOverallCorr, paste(SubjID, "_DMN_VrtxvMean.csv", sep=""), row.names = F)
  
  # Save the agreement between each slide and the overall partition
  if (Sliding) {
    write.csv2(allSlideCompared, file = paste(SubjID, "_slideWindowComparisons.csv", sep=""), row.names = F)
  
    # Write the actual sliding affiliations 
    write.csv2(slidingVals, file = paste(SubjID, "_slideWindowValues.csv", sep=""), row.names = F)
  }
  
  # Densities
  if (NA %in% densities$Density) {
    warning('Some densities not stored properly')
  } else {
    write.csv(densities, paste(SubjID, "_Densities.csv", sep=""), row.names = F)
  }
  
  # Write files for HCP
  HCPOut(Data = dmnval7mCommunities, MOI = "FiedlerBinary", SubjID = SubjID)  
  HCPOut(Data = dmnval7mCommunities, MOI = "FiedlerVec", SubjID = SubjID) 
  HCPOut(Data = dmnval7mCommunities, MOI = "Strength", SubjID = SubjID) 
  HCPOut(Data = dmnval7mCommunities, MOI = "varCoeff", SubjID = SubjID) 
  
}












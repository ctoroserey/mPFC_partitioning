###-------- Libraries and functions ---------------
Cols <- c("aquamarine4","#D9541A",rgb(190,190,190,100, maxColorValue = 255)) # left, right, interhemisphere

# In case I want to visualize stuff in the brain
labelCoords_vertex <- read.csv2('labelCoords_vertex.csv', sep = ",")[,2:6]
labelCoords_vertex <- transform(labelCoords_vertex, x = as.numeric(as.character(x)), y = as.numeric(as.character(y)), z = as.numeric(as.character(z)))

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
library(corrplot)


## Cohen's D for 2 groups
# for a more flexible approach, make the data input to be a list with entries for n groups, 
# then do length(list) for the number of groups. 
cohenD <- function(group1 = 1, group2 = 2){
  
  # means
  mean1 <- mean(group1, na.rm = T)
  mean2 <- mean(group2, na.rm = T)
  
  # variance
  var1 <- var(group1, na.rm = T)
  var2 <- var(group2, na.rm = T)
  
  # equation
  out <- (mean1 - mean2) / sqrt((var1 + var2)/2)
  
  return(out)
  
}

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
eigenVals <- function(Data = binnedCommunities_p$First) {
  
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
slidingWindow <- function(subjTS = ROI_timeSeries[[1]], mins = 15, jump = 1, Spectral = T, Modularity = F, ROIs = c("_7m_", vmPFC_labels)) {
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
  # Index: VI for variation of information, RI for the adjusted rand index, Cor for a Pearson correlation
  
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
    } else if (Index == "Cor") {
      for (subj in seq(nSubjects)) {
        for (subj2 in seq(nSubjects)) {
          indexMatrix[subj,subj2] <- cor(allVecs[, subj], allVecs[, subj2])
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
    } else if (Index == "Cor") {
      for (subj in seq(nSubjects)) {
        indexMatrix$Index[subj] <- cor(Data[[subj]][, MOI_indx[1]], Data[[subj]][, MOI_indx[2]])
      }
    }
  }
  
  return(indexMatrix)
  
}


###-------- For HCP data --------
setwd('./For_HCP')
temp <- list.files(pattern = '*Binary*')
FB <- lapply(temp, read.table)
Summed <- do.call(cbind, lapply(FB, "[[", 'V1')) - 1
FB_Means <- rowMeans(Summed)

write.table(file = 'propDMN_all.txt', FB_Means, row.names = F, col.names = F, dec = ".")

###-------- Partition comparison: all ------
setwd('../Summary/')
temp <- list.files()
subjList <- sapply(temp, substring, first=1, last=6)
Summaries <- lapply(temp, read.csv)
nSubj <- length(Summaries)

# Cross-method
# Don't worry about plotting, since that's mainly for a general visualization
RI_all_within <- comparePartitions(Data = Summaries, Index = "RI", MOI = c("FiedlerBinary","Membership"), nSubjects = length(Summaries), subjNames = subjList)
plot(density(RI_all_within$Index),
     xlim = c(0,1),
     main = "Distribution of Agreement Values Between Modularity and Spectral Partitioning (per subject)",
     xlab = "Adj. Rand Index")

# Across subjects
RI_all_between <- comparePartitions(Data = Summaries, Index = "RI", nSubjects = length(Summaries), subjNames = subjList)
posteriorIndx <- Summaries[[1]]$y < 0
tempComm_PCC <- lapply(Summaries, "[", i = posteriorIndx, j =)
tempComm_PFC <- lapply(Summaries, "[", i = !posteriorIndx, j =)
RI_PFC <- comparePartitions(Data = tempComm_PFC, Index = "RI", nSubjects = length(Summaries), subjNames = subjList)
RI_PCC <- comparePartitions(Data = tempComm_PCC, Index = "RI", nSubjects = length(Summaries), subjNames = subjList)
corrplot(RI_PCC,
         is.corr = F, 
         method = "color",
         title = "RI-based Agreement of Communities Among Every Pair of Subjects (PCC and PFC)",
         mar = c(0,0,3,0),
         type = "lower", 
         tl.pos = "lt", 
         tl.col = "black", 
         cl.lim = c(0,1))
corrplot(RI_PFC, 
         is.corr = F, 
         method = "color", 
         type = "upper", 
         mar = c(0,0,3,0),
         tl.pos = "lt", 
         tl.col = "black",
         cl.pos = "n", 
         add = T)

# Let's compare stuff
testPCC <- RI_PCC[lower.tri(RI_PCC)]
testPFC <- RI_PFC[lower.tri(RI_PFC)]
permTest <- permute(testPCC, testPFC, statType = mean, paired = T)
ES <- cohenD(testPCC, testPFC) # caveat: maybe not ideal for bounded values, like RI. BUT VI yields similar results.


###-------- Partition comparison: halves ------
# Load
setwd('../Summary_halves/')
temp <- list.files(pattern = "*_finalSummary.csv")
subjList_halves <- sapply(temp, substring, first=1, last=9)
Summaries_halves <- lapply(temp, read.csv)

temp <- list.files(pattern = "*_corHalves.csv")
Correlation_halves <- lapply(temp, read.csv, header = F)

# Compare all halves with RI
RI_all_between_halves <- comparePartitions(Data = Summaries_halves, Index = "RI", nSubjects = length(Summaries_halves), subjNames = subjList_halves)
corrplot(RI_all_between_halves, 
         is.corr = F, 
         title = "RI-based Agreement Between Day 1 and Day 2 (Overall)",
         mar = c(0,0,3,0),
         method = "color",
         tl.pos = "lt", 
         tl.col = "black")

# Compare PCC and PFC among all subjects
tempComm_PCC_halves <- lapply(Summaries_halves, "[", i = posteriorIndx, j =)
tempComm_PFC_halves <- lapply(Summaries_halves, "[", i = !posteriorIndx, j =)
RI_PFC_halves <- comparePartitions(Data = tempComm_PFC_halves, Index = "RI", nSubjects = length(Summaries_halves), subjNames = subjList_halves)
RI_PCC_halves <- comparePartitions(Data = tempComm_PCC_halves, Index = "RI", nSubjects = length(Summaries_halves), subjNames = subjList_halves)
corrplot(RI_PCC_halves, 
         is.corr = F,          
         title = "RI-based Agreement Between Day 1 and Day 2: PCC and mPFC",
         mar = c(0,0,3,0),
         method = "color",
         type = "lower", 
         tl.pos = "lt", 
         tl.col = "black", 
         cl.lim = c(-0.03,1))
corrplot(RI_PFC_halves, 
         is.corr = F, 
         mar = c(0,0,3,0),
         method = "color",
         type = "upper", 
         tl.pos = "lt", 
         tl.col = "black", 
         cl.pos = "n", 
         add = T)



###-------- Partition comparison: session ------
setwd('../Summary_sess/')
temp <- list.files()
subjList_sess <- sapply(temp, substring, first=1, last=9)
Summaries_sess <- lapply(temp, read.csv)

# Across subjects
RI_all_between_sess <- comparePartitions(Data = Summaries_sess, Index = "RI", nSubjects = length(Summaries_sess), subjNames = subjList_sess)
tempComm_PCC_sess <- lapply(Summaries_sess, "[", i = posteriorIndx, j =)
tempComm_PFC_sess <- lapply(Summaries_sess, "[", i = !posteriorIndx, j =)
RI_PFC_sess <- comparePartitions(Data = tempComm_PFC_sess, Index = "RI", nSubjects = length(Summaries_sess), subjNames = subjList_sess)
RI_PCC_sess <- comparePartitions(Data = tempComm_PCC_sess, Index = "RI", nSubjects = length(Summaries_sess), subjNames = subjList_sess)
corrplot(RI_PCC_sess, 
         is.corr = F, 
         mar = c(0,0,2,0),
         method = "color", 
         type = "lower", 
         tl.pos = "n", 
         tl.col = "black", 
         cl.lim = c(-0.03,1))
corrplot(RI_PFC_sess, 
         is.corr = F, 
         mar = c(0,0,2,0),
         title = "RI Values for PCC and PFC across rsfMRI Sessions",
         method = "color", 
         type = "upper", 
         tl.pos = "n", 
         tl.col = "black", 
         cl.pos = "n", 
         add = T)

# Let's compare stuff
testPCC_sess <- RI_PCC_sess[lower.tri(RI_PCC_sess)]
testPFC_sess <- RI_PFC_sess[lower.tri(RI_PFC_sess)]
permTest_sess <- permute(testPCC_sess, testPFC_sess, statType = mean, paired = T)
ES_sess <- cohenD(testPCC_sess, testPFC_sess) # caveat: maybe not ideal for bounded values, like RI. BUT VI yields similar results.


###-------- Stability per subject -------
setwd('../Sliding_window/')
temp <- list.files(pattern = '*Comparisons.csv')
SC <- lapply(temp, read.csv2)
temp <- list.files(pattern = '*Values.csv')
SVals <- lapply(temp, read.csv2)
  
# Overall
SC_all <- do.call(rbind, SC)
SC_plot <- ggplot(data = SC_all)
SC_plot <- SC_plot + geom_line(aes(x = Window, y = Index, group = as.factor(SubjID)), color = "slategray4",alpha = 0.3, show.legend = F) +
              geom_smooth(aes(x = Window, y = Index), color = Cols[2], fill = Cols[2]) +
              ylim(0,1) + 
              theme(text = element_text(size=35)) +
              labs(title = "RI of Each Sliding Window vs Overall Data", y = "Adj. Rand Index") +
              theme_classic()


# Mean RI across windows per subject
SC_summary <- aggregate(SC_all$Index, by = list(SC_all$SubjID), FUN = mean)
SC_summary$SD <- aggregate(SC_all$Index, by = list(SC_all$SubjID), FUN = sd)$x
colnames(SC_summary) <- list("SubjID", "RI", "SD")
SC_sd_plot <- ggplot(data = SC_summary, aes(x=SD)) + 
  geom_histogram(color = "black", fill = Cols[2], bins=50) + 
  xlim(0,1) +
  labs(title = "Standard Deviation of Sliding Window RIs Per Subject", x = "SD") +
  theme(text = element_text(size=35)) +
  theme_classic() 


# Save template
#ggsave(filename = 'testStability.pdf', plot = SC_PCCPFC_sd_plot, device = "pdf", width = 5.5, height = 5)


# Check stability for PCC and mPFC separately
# Get the sliding values for PCC and PFC separately
tempSW_PCC <- lapply(SVals, "[", i = posteriorIndx, j =)
tempSW_PFC <- lapply(SVals, "[", i = !posteriorIndx, j =)

SC_PCC <- data.frame()
SC_PFC <- data.frame()

for (subj in seq(nSubj)) {
  tempPCC <- numeric()
  tempPFC <- numeric()
  nWins <- ncol(tempSW_PCC[[subj]])
  for (Win in seq(nWins)) {
    tempPCC[Win] <- arandi(tempComm_PCC[[subj]]$FiedlerBinary, tempSW_PCC[[subj]][,Win], adjust = T)
    tempPFC[Win] <- arandi(tempComm_PFC[[subj]]$FiedlerBinary, tempSW_PFC[[subj]][,Win], adjust = T)
  }
  dfPCC <- data.frame(SubjID = rep(subjList[subj], nWins),
                      Window = seq(nWins),
                      Index = tempPCC)
  dfPFC <- data.frame(SubjID = rep(subjList[subj], nWins),
                      Window = seq(nWins),
                      Index = tempPFC)
  SC_PCC <- rbind(SC_PCC, dfPCC)
  SC_PFC <- rbind(SC_PFC, dfPFC)
}

SC_PCC_summary <- aggregate(SC_PCC$Index, by = list(SC_PCC$SubjID), FUN = mean)
SC_PCC_summary$SD <- aggregate(SC_PCC$Index, by = list(SC_PCC$SubjID), FUN = sd)$x
SC_PCC_summary$Region <- rep("PCC", nSubj)
SC_PFC_summary <- aggregate(SC_PFC$Index, by = list(SC_PFC$SubjID), FUN = mean)
SC_PFC_summary$SD <- aggregate(SC_PFC$Index, by = list(SC_PFC$SubjID), FUN = sd)$x
SC_PFC_summary$Region <- rep("mPFC", nSubj)
SC_PCCPFC_summary <- rbind(SC_PFC_summary, SC_PCC_summary)
colnames(SC_PCCPFC_summary) <- list("SubjID", "Index", "SD", "Region")
rm(SC_PCC_summary, SC_PFC_summary)

# plots
SC_PCCPFC_plot <- ggplot(data = SC_PCC)
SC_PCCPFC_plot <- SC_PCCPFC_plot + geom_line(aes(x = Window, y = Index, group = as.factor(SubjID)), color = Cols[2], alpha = 0.3, show.legend = F) +
  geom_smooth(aes(x = Window, y = Index), color = Cols[2], fill = Cols[2]) +
  ylim(0,1) + 
  labs(title = "RI of Each Sliding Window vs Overall Data: PCC vs mPFC", y = "Adj. Rand Index") +
  theme(text = element_text(size=35)) +
  geom_line(data = SC_PFC, aes(x = Window, y = Index, group = as.factor(SubjID)), color = Cols[1], alpha = 0.3, show.legend = F) +
  geom_smooth(data = SC_PFC, aes(x = Window, y = Index), color = Cols[1], fill = Cols[1]) +
  theme_classic()
SC_PCCPFC_plot

SC_PCCPFC_mean_plot <- ggplot(data = SC_PCCPFC_summary, aes(x=as.factor(Region), y=Index, fill = as.factor(Region))) + 
  geom_boxplot(show.legend = F) + 
  scale_fill_manual(values = Cols[seq(2)]) +
  ylim(0,1) +
  labs(title = "Mean of Sliding Window RIs Per Subject: PCC vs mPFC", x = "", y = "Adj. Rand Index") +
  theme(text = element_text(size=35)) +
  theme_classic() 
SC_PCCPFC_mean_plot

SC_PCCPFC_sd_plot <- ggplot(data = SC_PCCPFC_summary, aes(x=SD, fill = as.factor(Region))) + 
  scale_fill_manual(values = Cols[seq(2)]) +
  geom_density(color = "black", alpha = 0.6, show.legend = F) + 
  xlim(0,1) +
  labs(title = "Standard Deviation of Sliding Window RIs Per Subject: PCC vs mPFC", x = "SD") +
  theme(text = element_text(size=35)) +
  theme_classic() 
SC_PCCPFC_sd_plot

# Putative steps to analyze the SNR of the full time series (PCC vs mPFC)
# Flexibility vs SNR
# Are vertices with lower SNR switching communities more often?
SNR_vs_Flex <- sapply(Summaries, function(data) cor(data$varCoeff, data$slideFlexibility))

# Compare SNR between PCC and mPFC per subject (uneven sample sizes, perm better?)
SNR_PCCvsPFC_perm <- sapply(seq_along(tempComm_PCC), function(x) permute(tempComm_PCC[[x]]$varCoeff, tempComm_PFC[[x]]$varCoeff)$Pval)
SNR_PCCvsPFC_ES <- sapply(seq_along(tempComm_PCC), function(x) cohenD(tempComm_PCC[[x]]$varCoeff, tempComm_PFC[[x]]$varCoeff))

# Get the mean SNR for each region and the effect size for the comparison 
SNR_meanPCC <- sapply(tempComm_PCC, function(data) mean(data$varCoeff))
SNR_meanPFC <- sapply(tempComm_PFC, function(data) mean(data$varCoeff))
SNR_meanComp <- data.frame(Region = c(rep("PCC", nSubj), rep("mPFC", nSubj)),
                           SNR = c(SNR_meanPCC, SNR_meanPFC))

# Plot the mPFC SNR vs PCC SNR for each subject (one figure containing all subjects)
snrmeanComp_plot <- ggplot(data = SNR_meanComp, aes(Region, SNR)) +
  geom_boxplot() +
  labs(title = "Mean Suject-wise Coefficient of Variation For PCC and mPFC", x = "", y = "SNR") +
  theme(text = element_text(size=35)) +
  theme_classic()
snrmeanComp_plot


###-------- Correlation vs community detection (what do we gain from this?)------
# Compare correlation with community approach on mPFC
tempComm_PFC_halves <- lapply(seq_along(tempComm_PFC_halves), function(x) cbind(tempComm_PFC_halves[[x]], Correlation_halves[[x]]))
tempComm <- sapply(seq(1, 20, by = 2), function(x) cor(tempComm_PFC_halves[[x]]$FiedlerVec, tempComm_PFC_halves[[x+1]]$FiedlerVec))
tempCor <- sapply(seq(1, 20, by = 2), function(x) cor(Correlation_halves[[x]]$V1, Correlation_halves[[x+1]]$V1))
prepostComp <- data.frame(Community = tempComm, Correlation = tempCor)
prepostComp <- melt(prepostComp)
colnames(prepostComp) <- c("Approach", "Correlation")
prepostComp_test <- permute(tempComm, tempCor, paired = T)


## Between halves analysis for correlation and communities (mPFC only)
D1 <- tempComm_PFC_halves
D2 <- lapply(tempComm_PFC_halves, function(data) transform(data, V1 = FiedlerVec)) # the comparison function can only take 1 column name for this stuff
Data <- c(D1,D2)
allHalfCombs <- comparePartitions(Data = Data, MOI = "V1", nSubjects = 40, Index = "Cor", subjNames = c(subjList_halves, subjList_halves))
corrplot(allHalfCombs, 
         title = "Comparing All Community and Correlation mPFC Partitions Across All Subjects",
         method = "color", 
         tl.col = "black", 
         mar = c(0,0,3,0))
segments(0.5,0.5,0.5,40.5, lwd = 2)
segments(0.5,0.5,40.5,0.5, lwd = 2)
segments(0.5,40.5,40.5,40.5, lwd = 2)
segments(40.5,0.5,40.5,40.5, lwd = 2)
segments(0.5,20.5,40.5,20.5, lwd = 2)
segments(20.5,0.5,20.5,40.5, lwd = 2)

# Day 1 vs Day 2 Spatial Correlation Comparison
prepostComp_plot <- ggplot(data = prepostComp, aes(Approach, Correlation)) +
  geom_boxplot() +
  labs(title = "mPFC Spatial Correlation Between Day 1 and Day 2: Correlation vs Community Approaches", x = "", y = "Correlation Coefficient") +
  ylim(0,1) +
  theme(text = element_text(size=35)) +
  theme_classic()
prepostComp_plot
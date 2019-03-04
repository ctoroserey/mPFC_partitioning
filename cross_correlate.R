#!/usr/bin/env Rscript

# cross correlate the mean activity of a subject's two networks

# libraries & functions
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

# simple function to get the mean time series of a community 
mean.activity <- function(data = NA, mask = NA) {
  
  #if (is.na(data)) {stop("You need to input some data")}
  #if (is.na(mask)) {print("No mask provided, will get the overall mean of all vertices")}
  
  temp <- data[mask, ]
  temp <- colMeans(temp)
  
  return(temp)
  
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


# assign a subject
SubjID <- commandArgs(trailingOnly=TRUE)
write(paste("Cross-correlating communities for subject", SubjID), stdout())

# load session data
setwd('./time_series')
temp <- list.files(pattern = paste(SubjID, "_Session*", sep=""))
timeSeries_sess <- lapply(temp, fread, header=F)
timeSeries_sess <- lapply(timeSeries_sess, data.matrix)

# load the resulting mapping and produce indices for both nets (thresholded at the group mean of |0.01|)
# mComm2 is our DMN
map <- read_csv(paste("../For_HCP/", SubjID, "_FiedlerVec_dataforCifti.txt", sep=""), col_names = F, col_types = cols())[[1]]
map[map == -1] <- 0
comm1 <- which(map < -0.01)
comm2 <- which(map > 0.01)

# get the mean of the communities per session
mComm1 <- lapply(timeSeries_sess, mean.activity, mask = comm1)
mComm2 <- lapply(timeSeries_sess, mean.activity, mask = comm2)

# figure out how to regress out the autocorrelated data to avoid inflated estimates (residuals from ARX?)

# cross-correlate and save into a tibble to plot
sessLbls <- c("One", "Two", "Three", "Four")
Cols <- c("grey20", "grey40", "grey60", "grey80")
t1 <- sapply(seq_along(mComm1), function(x) {ccf(mComm2[[x]], mComm1[[x]], plot = F)$acf})
t2 <- sapply(seq_along(mComm1), function(x) {ccf(mComm2[[x]], mComm1[[x]], plot = F)$lag})
ccfData <- tibble(Session = rep(factor(sessLbls, levels = sessLbls), each = nrow(t1)),
                  Lag = matrix(t2, ncol = 1)[, 1],
                  CCF = matrix(t1, ncol = 1)[, 1])

# plot all dists
ccfCI <- quantile(bootstrap(ccfData$CCF), 0.975)
p <- ggplot(aes(Lag, CCF, group = Session, color = Session), data = ccfData) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 0, size = 1, color = "#D9541A") +
  geom_hline(yintercept = abs(ccfCI), linetype = "dashed", color = "#D9541A") +
  geom_hline(yintercept = abs(ccfCI) * -1, linetype = "dashed", color = "#D9541A") +
  scale_color_manual(values = Cols) +
  theme_classic()

# save plot
ggsave(paste("../CCF_", SubjID, ".png", sep = ""),
       height = 4,
       width = 6,
       p,
       device = "png")

print("Cross-correlation plots for each session produced.")




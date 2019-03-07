# produce histograms for the peaks of each ccf across subjects

# libraries & functions
suppressMessages(library(tidyverse))

# load data
setwd('./cross_correlate/')
temp <- list.files(pattern = ".csv")
timeSeries_sess <- lapply(temp, read.csv, header = F)
nSubj <- length(timeSeries_sess)

# save the lag positions and trim mDFs
lag <- timeSeries_sess[[1]]$V1
timeSeries_sess <- lapply(timeSeries_sess, "[", i = , j = 2:ncol(timeSeries_sess[[1]]))


## mean correlation per peak
## this is a very rough estimate, because it disregards the variance within subject
# average within subject and then between subjects + SE
mCorr_subj <- do.call(rbind, lapply(timeSeries_sess, rowMeans))
mCorr_all <- tibble(Lag = lag,
                    meanCorrelation = colMeans(mCorr_subj),
                    SE = apply(mCorr_subj, 2, sd) / sqrt(nSubj)) %>%
             mutate(lower = meanCorrelation - SE,
                    upper = meanCorrelation + SE)

# plot
(p <- ggplot(aes(Lag, meanCorrelation), data = mCorr_all) +
        geom_line(color = "grey50") +
        geom_point(size = 1.5, color = "grey50") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#D9541A") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#D9541A") +
        theme_classic())


## peak counts
# get the absolute peak per subject
peaksAbs <- lapply(timeSeries_sess, function(data) {
    as.numeric(apply(data, 2, function(x) {which(abs(x) == max(abs(x)))}))
  }) %>% 
  lapply(function(x) {lag[x]}) %>%
  do.call(cbind, .) 
  

# get the positive peak per subject
peaksPos <- lapply(timeSeries_sess, function(data) {
    as.numeric(apply(data, 2, function(x) {which(x == max(x))}))
  }) %>% 
  lapply(function(x) {lag[x]}) %>%
  do.call(cbind, .) 

# get the negative peak per subject
peaksNeg <- lapply(timeSeries_sess, function(data) {
    as.numeric(apply(data, 2, function(x) {which(x == min(x))}))
  }) %>% 
  lapply(function(x) {lag[x]}) %>%
  do.call(cbind, .) 


# simple histograms
par(mfrow = c(1,3))
nPeaks <- ncol(peaksAbs) * nrow(peaksAbs)
hist(peaksNeg, breaks = nPeaks)
abline(v = 0, col = "red")
hist(peaksAbs, breaks = nPeaks)
abline(v = 0, col = "red")
hist(peaksPos, breaks = nPeaks)
abline(v = 0, col = "red")

# simple script to count the number of pos/neg correlations that go into SP
# this can be used to compute densities later on

library(data.table)

# fisher r to z transformation, plus p-value adjustment with FDR
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
  transformed$pvals <- 2 * pnorm(abs(z.vec), 0 , sqrt(1 / (n - 3)), lower.tail = F) # no lower tail, but absolute Zs = two tailed
  
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

# get filenames
files <- dir(pattern = "_timeSeries.csv")
subjList <- sapply(files, substring, first = 1, last = 6)

# get sample summary file to get vertex indices
vertexIndex <- read.csv('./Analyzed/Summary/100307_finalSummary.csv')
vertexIndex <- vertexIndex$Vertex

# vectors to store correlations and other stuff Joe asked for later on
nTR <- as.numeric()
neg <- as.numeric()
pos <- as.numeric()
propNonsig <- as.numeric()
propSigpos <- as.numeric()
propSigneg <- as.numeric()
  
for (subj in files) {
  
  
  write(paste("Analyzing data for file", subj), stdout())
  
  # load data for 1 subject
  tempTS <- fread(subj, header = F)
  
  # select only relevant vertices
  tempTS <- tempTS[vertexIndex, ]
  
  # get the pairwise correlations
  r <- cor(t(tempTS))
  rUpper <- r[upper.tri(r)] # we only want to count across pairs once
  
  # threshold to get the n of significant/non-sig correlations remaining
  transfMat <- fisherTanh(Data = r)
  transfMat$tanhZ[transfMat$adjustPvals > 0.05] <- 0
  corrMat <- transfMat$tanhZ
  corrMat <- corrMat[upper.tri(corrMat)] 

  # store nrelevant data
  nTR <- c(nTR, ncol(tempTS))
  neg <- c(neg, sum(rUpper < 0)) 
  pos <- c(pos, sum(rUpper > 0))
  propNonsig <- c(propNonsig, sum(corrMat == 0))
  propSigpos <- c(propSigpos, sum(corrMat > 0))
  propSigneg <- c(propSigneg, sum(corrMat < 0)) 
  
} 

# store in data frame and write to csv
df <- data.frame(SubjID = subjList,
                 nTR = nTR,
                 Positive = pos,
                 Negative = neg,
                 propNonsignificant = propNonsig,
                 propSigpositive = propSigpos,
                 propSignegative = propSigneg)

write.table(file = "corrCounts.csv", df, row.names = F)
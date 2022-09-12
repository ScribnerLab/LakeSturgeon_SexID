#!/usr/bin/env bash
#!/usr/bin/env R

#setwd("Documents/Projects/Lake_Sturgeon_SexID/")

# Load Data
# Columns are:
# 1-WellNumber (1,2,3, etc), 2-WellID (A1,A2,etc), 3-Reading, 4-Temperature, 5-Fluorescence, and 6-Derivative
# These are the standard output columns in the raw melt curve data output by QuantStudio
# You'll just need to strip away the header and rename the columns

infile <- "PO2.txt"
data1 <- read.table(file = infile, header = T)
data1$WellID <- as.character(data1$WellID)

# Load Model
model <- readRDS("SexAssignmentModel.RDS")

# Define function for retrieving summary statistics 
getPeaks <- function(x){
  samps <- unique(x$WellID)
  out <- matrix(nrow = length(samps), ncol = 4)
  for(i in 1:length(samps)){
    out[i,1] <- as.character(samps[i])
    sub <- x[which(x$WellID == samps[i]),]
    t1 <- sub[which(sub$Temperature >73 & sub$Temperature < 77),]
    t2 <- sub[which(sub$Temperature >77 & sub$Temperature < 81),]
    t3 <- sub[which(sub$Temperature >60 & sub$Temperature < 70),]
    out[i,2] <- max(t1$Derivative)
    out[i,3] <- max(t2$Derivative)
    out[i,4] <- mean(t3$Fluorescence)
  }
  prod <- data.frame(WELL = as.character(out[,1]),FEMALE_73_77 = as.numeric(out[,2]),MALE_77_81 = as.numeric(out[,3]),BASELINE_FLUOR_60_70 = as.numeric(out[,4]) )
  return(prod)
}

# Get the summary statisitcs for each well
data2 <- getPeaks(data1)

# Apply the model and posterior threshold
threshold = 0.99
pred <- predict(model, data2)
pred <- data.frame(POSTERIOR = pmax(pred$posterior[,2],pred$posterior[,1]), GENETIC_SEX = pred$class)
pred[which(pred[,1] < threshold),2] <- NA

# Piece together the output. Set samples with low baseline intensity to NAs.
out <- data.frame(WELL = data2$WELL, FEMALE_PEAK = data2$FEMALE_73_77, MALE_PEAK = data2$MALE_77_81, BASELINE = data2$BASELINE_FLUOR_60_70, ASSIGNMENT = pred$GENETIC_SEX)
out$ASSIGNMENT[which(out$BASELINE < 1)] <- NA

# Write the output
write.table(out, file = "SexAssignment_Out.txt",sep = "\t", quote = F, row.names = F, col.names = T)

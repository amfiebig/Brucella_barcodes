#S3_Code_counts_by_LN

# This code merges codes in each LN and each animal.  Based on order of samples in 'sample_names.txt'.  
# User would need to modify sampleNames variable or sample_names.txt for other sets of data. 
# This code will calculate some summary values for each BC in the Parotid LN samples,
# merge all and tally how many nLN contain each BC and how many codes are identified in ALL LNs,
# and add counts data or fractional data.

library(tidyverse)
rm(list = ls()) #clear environment
dir.create("LN_lists")

LNnames <- c("LN1L", "LN1R", "LN2L", "LN2R", "LN3L", "LN3R", "LN4L", "LN4R", "LN5L", "LN5R", "LN6L", "LN6R")
CattleNames <- c("C1", "C2", "C3", "C4", "C5", "C6", "ALL")
groups <- c(LNnames, CattleNames)
sampleNames <- read_table2("sample_names.txt", col_names = TRUE)

#make table to fill
LNData <- tibble(groups, nCodes =0,  nCodesCnt1 = 0, fracCodesCnt1 = 0, R_L_tot =0, R_L_shared = 0)


###__Animal_1

# 1.set up with first join for LN 1 L
name <- sampleNames[1,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN1L <- samp_BC

# 2.merge samples from parotid LNs
Slist<- c(2:4)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN1L <- full_join(LN1L, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
#how many reads in all samples in the group
MeanTotReads <- mean(c((colSums(LN1L[,2], na.rm = TRUE)), (colSums(LN1L[,3], na.rm = TRUE)), (colSums(LN1L[,4], na.rm = TRUE)), (colSums(LN1L[,5], na.rm = TRUE))))
#how many samples contain each barcode
LN1L$PTcount <- apply(LN1L[2:5], 1, function(x) length(which(x>0)))
#add summary stats to the end of the table: fraction of samples with code, total average and median counts, mean and median fractional counts
LN1L <- mutate(LN1L, "fracPT" = PTcount/4, 
                     "LnTotCounts" = rowSums(LN1L[2:5], na.rm = TRUE), 
                     "avgLNCounts" = LnTotCounts/4)    
LN1L$median <- apply(LN1L[2:5], 1, function(x) median(x, na.rm = TRUE))
LN1L <- mutate(LN1L, "median_fracLN" = median/MeanTotReads, 
                    "mean_fracLN" = avgLNCounts/MeanTotReads)  
                     
# 4 write LN list and fill in cells in the LNData table
write_tsv(LN1L, "LN_lists/LN1L.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[1,2] <- nrow(LN1L)
LNData[1,3] <- LN1L %>% filter(LnTotCounts == 1) %>% nrow()
LNData[1,4] <- LNData[1,3]/LNData[1,2]

##___Right

# 1.set up with first join for LN 1 R
name <- sampleNames[6,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN1R <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(7:9)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN1R <- full_join(LN1R, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN1R[,2], na.rm = TRUE)), (colSums(LN1R[,3], na.rm = TRUE)), (colSums(LN1R[,4], na.rm = TRUE)), (colSums(LN1R[,5], na.rm = TRUE))))
LN1R$PTcount <- apply(LN1R[2:5], 1, function(x) length(which(x>0)))
LN1R <- mutate(LN1R, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN1R[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN1R$median <- apply(LN1R[2:5], 1, function(x) median(x, na.rm = TRUE))
LN1R <- mutate(LN1R, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN1R, "LN_lists/LN1R.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[2,2] <- nrow(LN1R)
LNData[2,3] <- LN1R %>% filter(LnTotCounts == 1) %>% nrow()
LNData[2,4] <- LNData[2,3]/LNData[2,2]

#5 join R and L for Animal_1

full <- full_join(LN1L, LN1R, by = "barcode")
shared <- inner_join(LN1L, LN1R, by = "barcode")

LNData[13,5] <- nrow(full)
LNData[13,6] <- nrow(shared)


###___Animal2

# 1.set up with first join for LN 2 L
name <- sampleNames[10,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN2L <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(11:12)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN2L <- full_join(LN2L, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN2L[,2], na.rm = TRUE)), (colSums(LN2L[,3], na.rm = TRUE)), (colSums(LN2L[,4], na.rm = TRUE))))
LN2L$PTcount <- apply(LN2L[2:4], 1, function(x) length(which(x>0)))
LN2L <- mutate(LN2L, "fracPT" = PTcount/3, 
               "LnTotCounts" = rowSums(LN2L[2:4], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/3)    
LN2L$median <- apply(LN2L[2:4], 1, function(x) median(x, na.rm = TRUE))
LN2L <- mutate(LN2L, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  


# 4 write LN list
write_tsv(LN2L, "LN_lists/LN2L.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[3,2] <- nrow(LN2L)
LNData[3,3] <- LN2L %>% filter(LnTotCounts == 1) %>% nrow()
LNData[3,4] <- LNData[3,3]/LNData[3,2]


##___Right

# 1.set up with first join for LN 2 R
name <- sampleNames[13,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN2R <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(14:16)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN2R <- full_join(LN2R, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN2R[,2], na.rm = TRUE)), (colSums(LN2R[,3], na.rm = TRUE)), (colSums(LN2R[,4], na.rm = TRUE)), (colSums(LN2R[,5], na.rm = TRUE))))
LN2R$PTcount <- apply(LN2R[2:5], 1, function(x) length(which(x>0)))
LN2R <- mutate(LN2R, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN2R[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN2R$median <- apply(LN2R[2:5], 1, function(x) median(x, na.rm = TRUE))
LN2R <- mutate(LN2R, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN2R, "LN_lists/LN2R.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[4,2] <- nrow(LN2R)
LNData[4,3] <- LN2R %>% filter(LnTotCounts == 1) %>% nrow()
LNData[4,4] <- LNData[4,3]/LNData[4,2]

#5 join R and L for Animal2

full <- full_join(LN2L, LN2R, by = "barcode")
shared <- inner_join(LN2L, LN2R, by = "barcode")

LNData[14,5] <- nrow(full)
LNData[14,6] <- nrow(shared)


###___Animal_3

# 1.set up with first join for LN 3 L
name <- sampleNames[17,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN3L <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(18:20)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN3L <- full_join(LN3L, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN3L[,2], na.rm = TRUE)), (colSums(LN3L[,3], na.rm = TRUE)), (colSums(LN3L[,4], na.rm = TRUE)), (colSums(LN3L[,5], na.rm = TRUE))))
LN3L$PTcount <- apply(LN3L[2:5], 1, function(x) length(which(x>0)))
LN3L <- mutate(LN3L, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN3L[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN3L$median <- apply(LN3L[2:5], 1, function(x) median(x, na.rm = TRUE))
LN3L <- mutate(LN3L, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN3L, "LN_lists/LN3L.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[5,2] <- nrow(LN3L)
LNData[5,3] <- LN3L %>% filter(LnTotCounts == 1) %>% nrow()
LNData[5,4] <- LNData[5,3]/LNData[5,2]

##___Right

# 1.set up with first join for LN 3 R
name <- sampleNames[22,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN3R <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(23:25)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN3R <- full_join(LN3R, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN3R[,2], na.rm = TRUE)), (colSums(LN3R[,3], na.rm = TRUE)), (colSums(LN3R[,4], na.rm = TRUE)), (colSums(LN3R[,5], na.rm = TRUE))))
LN3R$PTcount <- apply(LN3R[2:5], 1, function(x) length(which(x>0)))
LN3R <- mutate(LN3R, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN3R[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN3R$median <- apply(LN3R[2:5], 1, function(x) median(x, na.rm = TRUE))
LN3R <- mutate(LN3R, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN3R, "LN_lists/LN3R.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[6,2] <- nrow(LN3R)
LNData[6,3] <- LN3R %>% filter(LnTotCounts == 1) %>% nrow()
LNData[6,4] <- LNData[6,3]/LNData[6,2]

#5 join R and L for Animal3

full <- full_join(LN3L, LN3R, by = "barcode")
shared <- inner_join(LN3L, LN3R, by = "barcode")

LNData[15,5] <- nrow(full)
LNData[15,6] <- nrow(shared)

###___Animal__4

# 1.set up with first join for LN 4 L
name <- sampleNames[26,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN4L <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(27:29)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN4L <- full_join(LN4L, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN4L[,2], na.rm = TRUE)), (colSums(LN4L[,3], na.rm = TRUE)), (colSums(LN4L[,4], na.rm = TRUE)), (colSums(LN4L[,5], na.rm = TRUE))))
LN4L$PTcount <- apply(LN4L[2:5], 1, function(x) length(which(x>0)))
LN4L <- mutate(LN4L, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN4L[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN4L$median <- apply(LN4L[2:5], 1, function(x) median(x, na.rm = TRUE))
LN4L <- mutate(LN4L, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  


# 4 write LN list
write_tsv(LN4L, "LN_lists/LN4L.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[7,2] <- nrow(LN4L)
LNData[7,3] <- LN4L %>% filter(LnTotCounts == 1) %>% nrow()
LNData[7,4] <- LNData[7,3]/LNData[7,2]

##___Right

# 1.set up with first join for LN 4 R
name <- sampleNames[30,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN4R <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(31:33)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN4R <- full_join(LN4R, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN4R[,2], na.rm = TRUE)), (colSums(LN4R[,3], na.rm = TRUE)), (colSums(LN4R[,4], na.rm = TRUE)), (colSums(LN4R[,5], na.rm = TRUE))))
LN4R$PTcount <- apply(LN4R[2:5], 1, function(x) length(which(x>0)))
LN4R <- mutate(LN4R, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN4R[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN4R$median <- apply(LN4R[2:5], 1, function(x) median(x, na.rm = TRUE))
LN4R <- mutate(LN4R, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN4R, "LN_lists/LN4R.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[8,2] <- nrow(LN4R)
LNData[8,3] <- LN4R %>% filter(LnTotCounts == 1) %>% nrow()
LNData[8,4] <- LNData[8,3]/LNData[8,2]

#5 join R and L for Animal_4

full <- full_join(LN4L, LN4R, by = "barcode")
shared <- inner_join(LN4L, LN4R, by = "barcode")

LNData[16,5] <- nrow(full)
LNData[16,6] <- nrow(shared)


###___Animal_5

# 1.set up with first join for LN 5 L
name <- sampleNames[34,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN5L <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(35:37)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN5L <- full_join(LN5L, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN5L[,2], na.rm = TRUE)), (colSums(LN5L[,3], na.rm = TRUE)), (colSums(LN5L[,4], na.rm = TRUE)), (colSums(LN5L[,5], na.rm = TRUE))))
LN5L$PTcount <- apply(LN5L[2:5], 1, function(x) length(which(x>0)))
LN5L <- mutate(LN5L, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN5L[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN5L$median <- apply(LN5L[2:5], 1, function(x) median(x, na.rm = TRUE))
LN5L <- mutate(LN5L, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN5L, "LN_lists/LN5L.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
LNData[9,2] <- nrow(LN5L)
LNData[9,3] <- LN5L %>% filter(LnTotCounts == 1) %>% nrow()
LNData[9,4] <- LNData[9,3]/LNData[9,2]

##___Right

# 1.set up with first join for LN 5 R
name <- sampleNames[38,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN5R <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(39:41)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN5R <- full_join(LN5R, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN5R[,2], na.rm = TRUE)), (colSums(LN5R[,3], na.rm = TRUE)), (colSums(LN5R[,4], na.rm = TRUE)), (colSums(LN5R[,5], na.rm = TRUE))))
LN5R$PTcount <- apply(LN5R[2:5], 1, function(x) length(which(x>0)))
LN5R <- mutate(LN5R, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN5R[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN5R$median <- apply(LN5R[2:5], 1, function(x) median(x, na.rm = TRUE))
LN5R <- mutate(LN5R, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN5R, "LN_lists/LN5R.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[10,2] <- nrow(LN5R)
LNData[10,3] <- LN5R %>% filter(LnTotCounts == 1) %>% nrow()
LNData[10,4] <- LNData[10,3]/LNData[10,2]

#5 join R and L for Animal_5

full <- full_join(LN5L, LN5R, by = "barcode")
shared <- inner_join(LN5L, LN5R, by = "barcode")

LNData[17,5] <- nrow(full)
LNData[17,6] <- nrow(shared)

###___Animal_6

# 1.set up with first join for LN 6 L
name <- sampleNames[42,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN6L <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(43:45)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN6L <- full_join(LN6L, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN6L[,2], na.rm = TRUE)), (colSums(LN6L[,3], na.rm = TRUE)), (colSums(LN6L[,4], na.rm = TRUE)), (colSums(LN6L[,5], na.rm = TRUE))))
LN6L$PTcount <- apply(LN6L[2:5], 1, function(x) length(which(x>0)))
LN6L <- mutate(LN6L, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN6L[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN6L$median <- apply(LN6L[2:5], 1, function(x) median(x, na.rm = TRUE))
LN6L <- mutate(LN6L, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN6L, "LN_lists/LN6L.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
LNData[11,2] <- nrow(LN6L)
LNData[11,3] <- LN6L %>% filter(LnTotCounts == 1) %>% nrow()
LNData[11,4] <- LNData[11,3]/LNData[11,2]

##___Right

# 1.set up with first join for LN 6 R
name <- sampleNames[47,1]
samp_BC <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
sample <- paste(name)
samp_BC <- samp_BC %>% rename(!!sample := "counts")
LN6R <- samp_BC

# 2.merge on samples from parotid LNs (exclude retropharyngeal samples)
Slist<- c(48:50)
for (i in Slist) {
  name_i <- sampleNames[i,1]
  samp_BC <- read_table2(paste("filtered_data/", name_i, ".codesCln", sep=""), col_names = TRUE)
  sample <- paste(name_i)
  samp_BC <- samp_BC %>% rename(!!sample := "counts")
  LN6R <- full_join(LN6R, samp_BC, by="barcode") 
}     

# 3.calculations on counts across the LN
MeanTotReads <- mean(c((colSums(LN6R[,2], na.rm = TRUE)), (colSums(LN6R[,3], na.rm = TRUE)), (colSums(LN6R[,4], na.rm = TRUE)), (colSums(LN6R[,5], na.rm = TRUE))))
LN6R$PTcount <- apply(LN6R[2:5], 1, function(x) length(which(x>0)))
LN6R <- mutate(LN6R, "fracPT" = PTcount/4, 
               "LnTotCounts" = rowSums(LN6R[2:5], na.rm = TRUE), 
               "avgLNCounts" = LnTotCounts/4)    
LN6R$median <- apply(LN6R[2:5], 1, function(x) median(x, na.rm = TRUE))
LN6R <- mutate(LN6R, "median_fracLN" = median/MeanTotReads, 
               "mean_fracLN" = avgLNCounts/MeanTotReads)  

# 4 write LN list
write_tsv(LN6R, "LN_lists/LN6R.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

LNData[12,2] <- nrow(LN6R)
LNData[12,3] <- LN6R %>% filter(LnTotCounts == 1) %>% nrow()
LNData[12,4] <- LNData[12,3]/LNData[12,2]

#5 join R and L for Animal_6

full <- full_join(LN6L, LN6R, by = "barcode")
shared <- inner_join(LN6L, LN6R, by = "barcode")

LNData[18,5] <- nrow(full)
LNData[18,6] <- nrow(shared)


####____trim and merge LN codes and summarise

LN_1L <- LN1L[,c(1,7:12)]
LN_1R <- LN1R[,c(1,7:12)]
LN_2L <- LN2L[,c(1,6:11)]
LN_2R <- LN2R[,c(1,7:12)]
LN_3L <- LN3L[,c(1,7:12)]
LN_3R <- LN3R[,c(1,7:12)]
LN_4L <- LN4L[,c(1,7:12)]
LN_4R <- LN4R[,c(1,7:12)]
LN_5L <- LN5L[,c(1,7:12)]
LN_5R <- LN5R[,c(1,7:12)]
LN_6L <- LN6L[,c(1,7:12)]
LN_6R <- LN6R[,c(1,7:12)]



#_______merge all LN wide with total counts for each LN

LN_BCs <- LN_1L[,c("barcode", "LnTotCounts")]
LN_BCs <- LN_BCs %>% rename(LN_1L = "LnTotCounts")
LN_BCs <- full_join(LN_BCs, LN_1R[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_1R = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_2L[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_2L = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_2R[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_2R = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_3L[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_3L = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_3R[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_3R = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_4L[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_4L = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_4R[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_4R = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_5L[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_5L = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_5R[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_5R = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_6L[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_6L = "LnTotCounts")
LN_BCs <- full_join(LN_BCs , LN_6R[,c("barcode", "LnTotCounts")], by = "barcode")
LN_BCs <- LN_BCs %>% rename(LN_6R = "LnTotCounts")

write_tsv(LN_BCs, "R_output/CodeCounts_by_lymphnode.codes", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")


#__________Inner merge all LN wide with total counts for each LN _ to count num codes in ALL LN samples

All_PT_BCs <- LN_1L[,c("barcode", "LnTotCounts")]
All_PT_BCs <- All_PT_BCs %>% rename(LN_1L = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs, LN_1R[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_1R = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_2L[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_2L = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_2R[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_2R = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_3L[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_3L = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_3R[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_3R = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_4L[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_4L = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_4R[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_4R = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_5L[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_5L = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_5R[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_5R = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_6L[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_6L = "LnTotCounts")
All_PT_BCs <- inner_join(All_PT_BCs , LN_6R[,c("barcode", "LnTotCounts")], by = "barcode")
All_PT_BCs <- All_PT_BCs %>% rename(LN_6R = "LnTotCounts")


#___________

LNData[19,5] <- nrow(LN_BCs)  # add tot num codes
LNData[19,6] <- nrow(All_PT_BCs)  # add num chared codes

write_tsv(LNData, "R_output/Table_S1b_LN_Summary_Stats.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

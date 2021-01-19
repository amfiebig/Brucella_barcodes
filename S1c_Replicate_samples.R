#S1c_Replicate_samples.R

# For this code, the MultiCodes.pl output for the replicate samples should be in a "raw_data/replicates" folder in the working directory.

# Start by making a list of replicate samples (copy1, copy2 and samples).  
# Read in inoculum list for filtering.
# Open each raw '.codes' file from the "raw_data/replicates" folder.  Filter each using same approach as "S1_initial_processing"
# Full_join replicates.  Write files with individual filtered replicate counts and codes, and joined codes files. 

# Used the joined files in "replicates/pairs" to make graphs in prism for Supplemental Figure 4.

#~~~~~~~~~~~~~~~

library(tidyverse)
dir.create("replicates")
dir.create("replicates/filtered")
dir.create("replicates/pairs")

#list of replicate samples names 
copy1 <-c('10_L_D1_41', '1705_L_A1_44', '1707_R_A1_16', '1707_R_B1_19', '1708_R_C1_30', '1709_R_D1_58')
copy2 <-c('10_L_D2_42', '1705_L_A2_45', '1707_R_A2_17', '1707_R_B2_20', '1708_R_C2_32', '1709_R_D2_59')
samples <- c('10_L_D', "1705_L_A", "1707_R_A", "1707_R_B", '1708_R_C', "1709_R_D")

numsamp <- length(samples)  #how many samples in set

#read in inoc list for filtering
InocFull <- read_table2("R_output/Inoc_barcode_list.txt", col_names = TRUE)

threshold <- 99

sampleindex <- 1:(numsamp)   #make array of sampleindex numbers corresponding to repliate samples

#filter both replicates and write filtered .codesCln and .countsCln files
for (x in sampleindex) {
  samp1 <- paste(copy1[x]) #original file name replicate 1
  samp2 <- paste(copy2[x]) #original file name replicate 2
  name <- paste(samples[x]) #sample name
  
  #read in rep1 file by original name
  tible1 <- read_table2(paste("raw_data/replicates_raw/", samp1, ".codes", sep=""), col_names = TRUE)
  tible1 <- rename(tible1, counts = 2)
  tibleSort <- arrange(tible1, desc(counts))
  bigcodes <- tibleSort %>% filter(counts > threshold)
  numBigcodes <- nrow(bigcodes)
  
  # make list of off-by-one from biggest code, and remove to generate tibleClean
  i=1
  rmrows <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
  rmrows <- rmrows[rmrows > i]
  
  #loop 1.2 (in loop 1) - go through rest of big codes and remove off-by-one codes
  for (i in c(2:numBigcodes)) {
    rmrows_i <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
    rmrows_i <- rmrows_i[rmrows_i > i]
    rmrows <- c(rmrows,rmrows_i)
  }
  
  rmrows <- rmrows[!duplicated(rmrows)] #remove duplicates
  rmrows <- sort(rmrows, decreasing = FALSE)  #sort in order
  
  codes1 <- tibleSort[-rmrows,] #remove rows containg "off-by-ones"
  codes1 <- inner_join(codes1, InocFull, by = "barcode") # filter out codes not in Inoc
 
  counts1 <- codes1 %>% group_by(counts) %>% summarise(nCodes = length(barcode))
 
  # write new counts and codes files with new name
  write_tsv(counts1, paste("replicates/filtered/", name, "_1.countsCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(codes1, paste("replicates/filtered/", name, "_1.codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  #read in **rep2** file by original name
  tible2 <- read_table2(paste("raw_data/replicates_raw/", samp2, ".codes", sep=""), col_names = TRUE)
  tible2 <- rename(tible2, counts = 2)
  tibleSort <- arrange(tible2, desc(counts))
  bigcodes <- tibleSort %>% filter(counts > threshold)
  numBigcodes <- nrow(bigcodes)
  
  # make list of off-by-one from biggest code, and remove to generate tibleClean
  i=1
  rmrows <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
  rmrows <- rmrows[rmrows > i]
  
  #loop 2.2 (in loop 1) - go through rest of big codes and remove off-by-one codes
  for (i in c(2:numBigcodes)) {
    rmrows_i <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
    rmrows_i <- rmrows_i[rmrows_i > i]
    rmrows <- c(rmrows,rmrows_i)
  }
  
  rmrows <- rmrows[!duplicated(rmrows)] #remove duplicates
  rmrows <- sort(rmrows, decreasing = FALSE)  #sort in order
  
  codes2 <- tibleSort[-rmrows,]
  codes2 <- inner_join(codes2, InocFull, by = "barcode") # filter out codes not in Inoc
   
  counts2 <- codes2 %>% group_by(counts) %>% summarise(nCodes = length(barcode))

  # write new counts and codes files with new name
  write_tsv(counts2, paste("replicates/filtered/", name, "_2.countsCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(codes2, paste("replicates/filtered/", name, "_2.codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  #make join replicate pairs
  codes1 <- rename(codes1, rep_1 = 2)
  codes2 <- rename(codes2, rep_2 = 2)
  
  codes_full <- full_join(codes1, codes2, by = "barcode")
  codes_inner <- inner_join(codes1, codes2, by = "barcode")
  
  write_tsv(codes_full, paste("replicates/pairs/", name, ".codesCln", sep=""), na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
  
 }


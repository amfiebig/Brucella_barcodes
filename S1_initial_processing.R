# S1_initial_processing   

# This code is to open .counts, sort by count, search for and remove off-by-one codes, iteratively for all codes 
# that are counted *100* times or more in the inoculum (loop 1).  Then makes a merged list of inoculum codes.  
# Then in loop 2, filters all lymph node samples as above to remove off-by-ones of abundant (>99 reads) codes, 
# and remove codes not present in the inoculum.  After filtering, a data sheet containing simple information on the codes 
# and counts in each sample is generated.  In addition, new filtered codes and counts files are generated for each sample.
# A summary stats table is generated for all samples.
# Then in the final section of this code, samples are assembled (joined), to make a file with code counts for each LN or Inoc sample, 
# to make a ranked list of similar file, with counts ranked by abundance, to make assemble codes from all LN pieces, 
# from all PT LN samples, and finally from the PT LN samples before filtering out codes not detected in the inoc. 


# To run this code, one needs (1) a text file with sample names (called 'sample_names.txt') with at least two columns new_name and sample_name, corresponding to output and raw_input names, 
# in the working directory, and (2) the .codes and .counts files output from MultiCodes.pl (Wetmore et al) in a sub directory called "raw_data"

# _________________________________________

library(tidyverse)

getwd()
dir.create("R_output")  #to contain summary files
dir.create("filtered_data") #to contain filtered data lists from each sample
dir.create("filtered_data/off-by-one/") #to contain LN code lists after off-by-one filter, but without the inoculum filter
dir.create("threshold_lists") #to contain list of abundant codes in each sample


#read sample names with the base of the original file name in column 2 and new simple sample name in column 1
SampleNameTable <- read_table2("sample_names.txt", col_names = TRUE)

#make list of new names
SampleName <- SampleNameTable %>% pull(new_name)
numsamp <- length(SampleName)  #how many samples in set

threshold <- 99

#build AllData summary table
AllData <- tibble(SampleName, nCounts = 0, nCodes =0, nBigCodes = 0, 
                  nCountsClean = 0, nCodesClean = 0, nBigCodesCln = 0,
                  maxCode = 0, fracMaxCode = 0, nCodes2prcnt = 0, nCodesCnt1 = 0, 
                  fracCodesCnt1 = 0, fracCodesRemov = 0, fracCountsRemov =0)

##### loop 1 - open and filter INOC raw files

#make array of sampleindex numbers corresponding INOC samples, in this case the last 3 samples in the list
sampleindex <- (numsamp-2):numsamp   

for (x in sampleindex) {
  samp <- paste(SampleNameTable[x,2]) #original file name
  name <- paste(SampleNameTable[x,1]) #revised simple name
  
  #read in each file by original name
  tible1 <- read_table2(paste("raw_data/", samp, ".codes", sep=""), col_names = TRUE)
  tible1 <- rename(tible1, counts = 2)   #change column 2 to 'counts'
  tibleSort <- arrange(tible1, desc(counts)) #sort by counts, highest to lowest

  bigcodes <- tibleSort %>% filter(counts > threshold)  #make list of codes counted more than the threshold number of times
  numBigcodes <- nrow(bigcodes)   
  
  # make list of rows corresponding to off-by-one from biggest code
  i=1
  rmrows <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE) #find codes off-by-one from code i
  rmrows <- rmrows[rmrows > i]    # exclude row i from rmrow list by keeping only rows greater than i
  
  #loop 1.2 (in loop 1) - go through rest of big codes and add to list of rows to remove
      for (i in c(2:numBigcodes)) {
        rmrows_i <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
        rmrows_i <- rmrows_i[rmrows_i > i]  # keep only rows greater than i on the rmrows list
        rmrows <- c(rmrows,rmrows_i)  #merge list for this code i with the full list
      }
  
  rmrows <- rmrows[!duplicated(rmrows)] #remove duplicates
  rmrows <- sort(rmrows, decreasing = FALSE)  #sort in order
  
  tibleClean <- tibleSort[-rmrows,] #make new table without off-by-one rows
    bigcodesCln <- tibleClean %>% filter(counts > threshold)
  counts_file <- tibleClean %>% group_by(counts) %>% summarise(nCodes = length(barcode)) #generate new .counts file

  # write new counts and codes files with new name
  write_tsv(counts_file, paste("filtered_data/", name, ".countsCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(tibleClean, paste("filtered_data/", name, ".codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  # add some basic parameters to the codes files
  tibleClean <- tibleClean %>%  
    mutate(fracCounts = counts/(sum(tibleClean$counts))) #adds fracCounts
  tibleClean <- tibleClean %>%
    mutate(cumFrac = cumsum(tibleClean$fracCounts)) #adds cumulative fracCounts
  tibleClean <- tibleClean %>%
    mutate(bcRank = seq.int(nrow(tibleClean))) #adds rank of barcode in sample
  
  tibleClean <- tibleClean%>%
    mutate(sample = name) # appends new sample name to the last column
  tibleClean <- tibleClean[c(6,1,2,3,4,5)] #moves sample name to first column
   
  tibleClean <- tibleClean %>% separate(sample, into = c("animal","side", "piece"), sep = "_")

  # write new codes file with additional paramaters (namePlus.codesCln)
  write_tsv(tibleClean, paste("filtered_data/", name, "_plus.codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  morethan2prcnt <- tibleClean %>% filter(tibleClean$fracCounts > 0.02)
  BigCodes <- tibleClean %>% filter(tibleClean$counts > threshold)
 
  write_tsv(morethan2prcnt, paste("threshold_lists/morethan2prcnt_", name, ".codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")

  rownum <- which(SampleName == name)
      AllData[rownum,2] <- sum(tibleSort$counts) #nCounts original
      AllData[rownum,3] <- nrow(tibleSort) #nCodes original
      AllData[rownum,4] <- numBigcodes #num codes >threshold in original
      AllData[rownum,5] <- sum(tibleClean$counts) #nCounts clean
      AllData[rownum,6] <- nrow(tibleClean) #nCodes clean
      AllData[rownum,7] <- nrow(bigcodesCln) #num codes >threshold in clean
      AllData[rownum,8] <- max(tibleClean$counts) #counts for maximum code
      AllData[rownum,9] <- max(tibleClean$counts)/sum(tibleClean$counts) #maxcode count/total clean read count
      AllData[rownum,10] <- nrow(morethan2prcnt) #num codes at least 2 percent
      AllData[rownum,11] <- counts_file[1,2]  #num codes counted only once
 }

##### Assemble codes for INOC samples

InocA <- read_table2(file="filtered_data/Inoc_x_A.codesCln", col_names = TRUE)
InocA <- InocA %>% rename(inocA = "counts") 
InocB <- read_table2(file="filtered_data/Inoc_x_B.codesCln", col_names = TRUE)
InocB <- InocB %>% rename(inocB = "counts") 
InocC <- read_table2(file="filtered_data/Inoc_x_C.codesCln", col_names = TRUE)
InocC <- InocC %>% rename(inocC = "counts") 

ABfull <- full_join(InocA, InocB, by = "barcode")
ABCfull <- full_join(ABfull, InocC, by = "barcode")
inocFull <- ABCfull[,1]

write_tsv(ABCfull, "R_output/All_inoc_samples.codes", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
write_tsv(inocFull, "R_output/Inoc_barcode_list.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

##### open and filter LN samples using both off-by-one of big codes criteria AND present in inoc criteria 

# loop 2 - LN - open and process raw files for LN filtering off-by_ones and removing codes not in inoc
sampleindex <- 1:(numsamp-3)   #make array of sampleindex numbers corresponding to LN samples (in this case all but last three samples)

for (x in sampleindex) {
  samp <- paste(SampleNameTable[x,2]) #original file name
  name <- paste(SampleNameTable[x,1]) #revised simple name
  
  #read in each file by original name
  tible1 <- read_table2(paste("raw_data/", samp, ".codes", sep=""), col_names = TRUE)
  tible1 <- rename(tible1, counts = 2)
  tibleSort <- arrange(tible1, desc(counts))
  
  bigcodes <- tibleSort %>% filter(counts > threshold)
  numBigcodes <- nrow(bigcodes)
  
  # make list of off-by-one from biggest code, and remove to generate tibleClean
  i=1
  rmrows <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
  rmrows <- rmrows[rmrows > i]
 
  #loop 2.2 (in loop 2) - go through rest of big codes and remove off-by-one codes
  for (i in c(2:numBigcodes)) {
    rmrows_i <- agrep(tibleSort$barcode[i], tibleSort$barcode, 1, value = FALSE)
    rmrows_i <- rmrows_i[rmrows_i > i]
    rmrows <- c(rmrows,rmrows_i)
  }
  
  rmrows<-rmrows[!duplicated(rmrows)] #remove duplicates
  rmrows <- sort(rmrows, decreasing = FALSE)  #sort in order
  
  tible_rm_OBO <- tibleSort[-rmrows,] #remove off-by-one rows
  tibleClean <- inner_join(tible_rm_OBO, inocFull, by = "barcode") # filter out codes not in Inoc
  
  bigcodesCln <- tibleClean %>% filter(counts > threshold)
  counts_file <- tibleClean %>% group_by(counts) %>% summarise(nCodes = length(barcode))
  
  # write new counts and codes files with new name, and a code list without the inoc filter (i.e. tible_rm_OBO)
  write_tsv(counts_file, paste("filtered_data/", name, ".countsCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(tibleClean, paste("filtered_data/", name, ".codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(tible_rm_OBO, paste("filtered_data/off-by-one/", name, ".codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  # add some basic parameters to the codes files
  tibleClean <- tibleClean %>%  
    mutate(fracCounts = counts/(sum(tibleClean$counts))) #adds fracCounts
  tibleClean <- tibleClean %>%
    mutate(cumFrac = cumsum(tibleClean$fracCounts)) #adds cumulative fracCounts
  tibleClean <- tibleClean %>%
    mutate(bcRank = seq.int(nrow(tibleClean))) #adds rank of barcode in sample
  
  tibleClean <- tibleClean%>%
    mutate(sample = name) # appends new sample name to the last column
  tibleClean <- tibleClean[c(6,1,2,3,4,5)] #moves sample name to first column
  
  tibleClean <- tibleClean %>% separate(sample, into = c("animal","side", "piece"), sep = "_")
  
  # write new codes file with additional paramaters (namePlus.codesCln)
  write_tsv(tibleClean, paste("filtered_data/", name, "_plus.codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  morethan2prcnt <- tibleClean %>% filter(tibleClean$fracCounts > 0.02)
  BigCodes <- tibleClean %>% filter(tibleClean$counts > threshold)
  
  write_tsv(morethan2prcnt, paste("threshold_lists/morethan2prcnt_", name, ".codesCln", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  
  #add data to AllData summary table
  rownum <- which(SampleName == name)
  AllData[rownum,2] <- sum(tibleSort$counts) #nCounts original
  AllData[rownum,3] <- nrow(tibleSort) #nCodes original
  AllData[rownum,4] <- numBigcodes #num codes >threshold in original
  AllData[rownum,5] <- sum(tibleClean$counts) #nCounts clean
  AllData[rownum,6] <- nrow(tibleClean) #nCodes clean
  AllData[rownum,7] <- nrow(bigcodesCln) #num codes >threshold in clean
  AllData[rownum,8] <- max(tibleClean$counts) #counts for maximum code
  AllData[rownum,9] <- max(tibleClean$counts)/sum(tibleClean$counts) #maxcode count/total clean read count
  AllData[rownum,10] <- nrow(morethan2prcnt) #num codes at least 2 percent
  AllData[rownum,11] <- counts_file[1,2]  #num codes counted one time
}

AllData$fracCodesCnt1 = (AllData$nCodesCnt1/AllData$nCodesClean)

AllData$fracCodesRemov = ((AllData$nCodes-AllData$nCodesClean)/AllData$nCodes)
AllData$fracCountsRemov = ((AllData$nCounts-AllData$nCountsClean)/AllData$nCounts)

AllData2 <- AllData %>% separate(SampleName, into = c("animal","side", "piece"), sep = "_")

write_tsv(AllData2, "R_output/Table_S1a_Summary_Stats.txt", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")


########### Assemble data for Figure S6.  Counts data for part A, and ranked counts for part B. 

###Assemble COUNTS data from LN samples to make frequency distribution of # codes counted 1, 2, 3, 4,...times (Figure S6a)

#make list of new names
SampleName <- SampleNameTable %>% pull(new_name)

#Set up first cycle
name <- paste(SampleNameTable[1,1]) #revised simple name
countsCln <- read_table2(paste("filtered_data/", name, ".countsCln", sep=""), col_names = TRUE)
sample <- paste(name)
countsClnLN <- countsCln %>% rename(!!sample := "nCodes")

# Add on remaining samples
  numsamp <- length(SampleName)  #how many samples in set
  sampleindex <- 2:(numsamp-3)   #make array of sampleindex numbers corresponding to number of samples

  for (x in sampleindex) {
    name <- paste(SampleNameTable[x,1]) #revised simple name
    countsCln <- read_table2(paste("filtered_data/", name, ".countsCln", sep=""), col_names = TRUE)
    sample <- paste(name)
    countsCln <- countsCln %>% rename(!!sample := "nCodes")
    countsClnLN <- full_join(countsClnLN, countsCln, by="counts") 
}

countsClnLN <- countsClnLN   %>% arrange(counts)

write_tsv(countsClnLN, "R_output/All_LN_quadrants.counts", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
    #plotted data in above tsv file in prism for figure S6A

### Assemble counts data for INOC samples to make frequency distribution of # codes counted 1, 2, 3, 4,...times (Figure S6a)

countsCln_A <- read_table2("filtered_data/inoc_x_A.countsCln", col_names = TRUE)
countsCln_A <- countsCln_A %>% rename(inocA = "nCodes") 
countsCln_B <- read_table2("filtered_data/inoc_x_B.countsCln", col_names = TRUE)
countsCln_B <- countsCln_B %>% rename(inocB = "nCodes") 
countsCln_C <- read_table2("filtered_data/inoc_x_C.countsCln", col_names = TRUE)
countsCln_C <- countsCln_C %>% rename(inocC = "nCodes") 

countsClnInoc <- full_join(countsCln_A, countsCln_B, by="counts") 
countsClnInoc <- full_join(countsClnInoc, countsCln_C, by="counts") 
countsClnInoc <- countsClnInoc %>% arrange(counts)

write_tsv(countsClnInoc, "R_output/All_inoc_samples.counts", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
#plotted data in above tsv file in prism for figure S6A

### Assemble codes for LN samples in ranked order (Fig S6B)
    
    #Set up first cycle
    name <- paste(SampleNameTable[1,1]) #new name
    codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
    sample <- paste(name)
    codesCln <- codesCln[,2]
    codesCln <- mutate(codesCln, row_number())
    codesCln <- codesCln %>% rename(!!sample := "counts", rank = "row_number()")
    codesLNranked <- codesCln[, c('rank','1_L_A')] #change column order
    
    # Add on remaining samples
    numsamp <- length(SampleName)  #how many samples in set
    sampleindex <- 2:(numsamp-3)   #make array of sampleindex numbers corresponding to number of LN samples (in this case all but the last 3)
    
    for (x in sampleindex) {
      name <- paste(SampleNameTable[x,1]) #revised simple name
      codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
      sample <- paste(name)
      codesCln <- codesCln[,2]
      codesCln <- mutate(codesCln, row_number())
      codesCln <- codesCln %>% rename(!!sample := "counts", rank = "row_number()")
      codesLNranked <- full_join(codesLNranked, codesCln, by="rank") 
    }

### Assemble codes for Inoc samples in ranked order

    name <- paste(SampleNameTable[51,1]) #revised simple name
    codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
    sample <- paste(name)
    codesCln <- codesCln[,2]
    codesCln <- mutate(codesCln, row_number())
    codesCln <- codesCln %>% rename(!!sample := "counts", rank = "row_number()")
    codesInocRanked <- codesCln[, c(2,1)] #change column order
    
    # Add on remaining samples
    numsamp <- length(SampleName)  #how many samples in set
    sampleindex <- 52:53   #make array of sampleindex numbers corresponding to number of samples
    
    for (x in sampleindex) {
      name <- paste(SampleNameTable[x,1]) #revised simple name
      codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
      sample <- paste(name)
      codesCln <- codesCln[,2]
      codesCln <- mutate(codesCln, row_number())
      codesCln <- codesCln %>% rename(!!sample := "counts", rank = "row_number()")
      codesInocRanked <- full_join(codesInocRanked, codesCln, by="rank") 
    }

### filter to sample more sparsely as the rank becomes bigger and counts lower.  This makes plotting easier b/c fewer points.
    rank <- c(1:1000,seq(1050, 10000, 50),seq(10500, 100000, 500), seq(105000, 1030000, 5000)) 
    keeprows <- as.data.frame(rank)
    
    ### make trimmed data sets and save.  Make plots in prism.
    codesInocRankedTrim <- inner_join(keeprows, codesInocRanked, by = "rank")
    write_tsv(codesInocRankedTrim, "R_output/Ranked_Barcode_counts_Inoc_trimmed.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
    codesLNrankedTrim <- inner_join(keeprows, codesLNranked, by = "rank")
    write_tsv(codesLNrankedTrim, "R_output/Ranked_Barcode_counts_LN_trimmed.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
    
    
    
### Assemble CODES for LN samples - for use in diversity analysis
    
    #Set up first cycle
    name <- paste(SampleNameTable[1,1]) #revised simple name
    codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
    sample <- paste(name)
    codesClnLN <- codesCln %>% rename(!!sample := "counts")
    
    # Add on remaining samples
    numsamp <- length(SampleName)  #how many samples in set
    sampleindex <- 2:(numsamp-3)   #make array of sampleindex numbers corresponding to number of samples
    
    for (x in sampleindex) {
      name <- paste(SampleNameTable[x,1]) #revised simple name
      codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
      sample <- paste(name)
      codesCln <- codesCln %>% rename(!!sample := "counts")
      codesClnLN <- full_join(codesClnLN, codesCln, by="barcode") 
    }
    
    write_tsv(codesClnLN, "R_output/All_LN_quadrants.codes", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
    
    #--- repeat but exclude retro samples
    sampleindex <- c(2:4,6:20,22:45,47:50)   #make array of sampleindex numbers corresponding to number of samples
    
    name <- paste(SampleNameTable[1,1]) #revised simple name
    codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
    sample <- paste(name)
    codesClnPT <- codesCln %>% rename(!!sample := "counts")
    
    for (x in sampleindex) {
      name <- paste(SampleNameTable[x,1]) #revised simple name
      codesCln <- read_table2(paste("filtered_data/", name, ".codesCln", sep=""), col_names = TRUE)
      sample <- paste(name)
      codesCln <- codesCln %>% rename(!!sample := "counts")
      codesClnPT <- full_join(codesClnPT, codesCln, by="barcode") 
    }
    
    write_tsv(codesClnPT, "R_output/All_PT_LN_quadrants.codes", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
    
    

#### Assemble list of PT LN barcodes filtered by off_by_one only (i.e. includes codes NOT in inoc) for Venn diagrams in Figure S5

    sampleindex <- c(2:4,6:20,22:45,47:50)   #make array of sampleindex numbers corresponding to appropriate samples
    
    name <- paste(SampleNameTable[1,1]) #revised simple name
    codesCln <- read_table2(paste("filtered_data/off-by-one/", name, ".codesCln", sep=""), col_names = TRUE)
    sample <- paste(name)
    PT_wI <- codesCln %>% rename(!!sample := "counts")
    
    for (x in sampleindex) {
      name <- paste(SampleNameTable[x,1]) #revised simple name
      codesCln <- read_table2(paste("filtered_data/off-by-one/", name, ".codesCln", sep=""), col_names = TRUE)
      sample <- paste(name)
      codesCln <- codesCln %>% rename(!!sample := "counts")
      PT_wI <- full_join(PT_wI, codesCln, by="barcode") 
    }

write_tsv(PT_wI, "R_output/All_PT_LN_samples_without_inoc_filter.codes", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")



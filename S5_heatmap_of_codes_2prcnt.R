#S5_heatmap_of_codes>2prcnt

#merge lists of all codes >2% in each sample (in this case, 160 codes)
#then extract the fractional counts from all samples for those barcodes in a wide format
#log10 transform the fractional counts and save file
#then open in excel to manually sort samp1 >2%, then samp2 >2%, then samp3 >2prcnt... 
  #then add a barcode number to preserve sort order
#then make heat map and export as .eps to label in Illustrator


library(tidyverse)
library(reshape2)
rm(list = ls()) #clear environment

SampleNameTable <- read_table2("sample_names.txt", col_names = TRUE)

#make list of new names
SampleName <- SampleNameTable %>% pull(new_name)

numsamp <- length(SampleName)  #how many samples in set

### (a) make a table of codes with abundances greater than 2% in a samples
# and then a list of the unique codes within in this group

  #start with first sample
  All2prcnt <- read_table2(paste("threshold_lists/morethan2prcnt_", SampleName[1], ".codesCln", sep=""), col_names = TRUE)

  #then append all the other samples  
  for (samp in SampleName[2:numsamp]) {
    part2prcnt <- read_table2(paste("threshold_lists/morethan2prcnt_", samp, ".codesCln", sep=""), col_names = TRUE)
    All2prcnt <- bind_rows(All2prcnt, part2prcnt)
    }
  write_tsv(All2prcnt, paste("R_output/BC_2prcnt_Occurances.codes", sep=""), na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")

  ## make non-redundant list of BC 
  BC_2prcnt <- All2prcnt %>% group_by(barcode) %>% summarise(occurances = n())
  BC_2prcnt <- BC_2prcnt %>% arrange(-occurances) #sort by occurances
  BC_2prcnt <- BC_2prcnt[,-2] #remove occurances in column 2 to leave just BC list
  write_tsv(BC_2prcnt, "R_output/BC_2prcnt_uniquecodes.txt", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  

### build a wide table to join counts for each code in the unique subset lists from each sample
  
  #build a table to join to , set up with first join
    samp <- SampleName[1]
    
    samp_BC_list <- read_table2(paste("filtered_data/", samp, "_plus.codesCln", sep=""), col_names = TRUE)
    samp_BC_list <- samp_BC_list[,c(4,6)] #keeps only columns with BC, fracCount
 
    BC_all_2prcnt <- left_join(BC_2prcnt, samp_BC_list, by="barcode")  
    BC_all_2prcnt <- rename(BC_all_2prcnt, !!samp := "fracCounts")
   
  #loop through remaining samples and collect fractional counts for codes on the '>2% list'
    for (samp in SampleName[2:numsamp]) {
      
      samp_BC_list <- read_table2(paste("filtered_data/", samp, "_plus.codesCln", sep=""), col_names = TRUE)
      samp_BC_list <- samp_BC_list[,c(4,6)] #keeps only columns with BC and fracCount
      
      BC_all_2prcnt <- left_join(BC_all_2prcnt, samp_BC_list, by="barcode")  #select codes from samp_BC_list that are in the >2% list
      BC_all_2prcnt <- rename(BC_all_2prcnt, !!samp := "fracCounts")
      }  
   
    BC_2prcnt_log10Frac <- BC_all_2prcnt
    BC_2prcnt_log10Frac[,2:54] <- log10(BC_2prcnt_log10Frac[,2:54])
    
  write_tsv(BC_all_2prcnt, "R_output/BC_2prcnt_gathered_fracCounts.txt", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
  write_tsv(BC_2prcnt_log10Frac, "R_output/BC_2prcnt_gathered_Log10fracCounts.txt", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")

##### here I opened the Log10fracCounts file in excel and sorted manually (because I didn't know how to do this in R)
  #I sorted the barcodes with the BC >2% in samp1 at top, then BC>2% in samp2, then BC>2% in samp3...
  #then added a column with a BC number (BC001-BC160) to preserve to sort order and resaved as BC_2prcnt_gathered_Log10fracCounts_Sorted.txt
     
###make heat map  
  abundance1 <- read_table2('R_output/BC_2prcnt_gathered_Log10fracCounts_Sorted.txt', col_names = TRUE)
  abundance1 <- abundance1[,-1] #removes the column with the BC sequence, leaves column with BC number
  melted_abundance1 = melt(abundance1, na.rm = TRUE) #converts wide to long data format
  
  ggheatmap1 = ggplot(data = melted_abundance1, aes(barcodeNum, variable, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "gray90", high = "gray0", mid = "gray55", 
                         midpoint = -3.5, limit = c(-7,0), space = "Lab", name="Log(frac)") +
    theme_minimal()+ 
    coord_fixed()+ 
    theme(
      #axis.text.x = element_text(angle = 90),
      axis.text.y = element_text(angle = 180),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())+
    labs(title="abundant codes")
  ggheatmap1
  
  #then export ggheatmap1 as .eps
  
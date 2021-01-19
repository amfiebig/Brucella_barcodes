#S4_Population_Diversity_Assessment
  #S4a_Diversity assessment with vegan
  #S4b_Bray-Curtis and Jaccard Dissimilarity matrices

#set_up: merge sample lists and transpose so rows are samples and columns are codes
#PartA: assesses some basic metrics for the populations
#PartB: compares populations to makes heat maps

library(tidyverse)
library(vegan)
library(reshape2)
rm(list = ls()) #clear environment

######-Set-up getting tables in correct format
#read table (with full merged barcodes and NA replaced with 0) and make dataframe 'BC_0'
InocBCs <- read_table2("R_output/All_inoc_samples.codes", col_names = TRUE)
LNpiecesBCs <- read_table2("R_output/All_LN_quadrants.codes", col_names = TRUE) #PT quads and retro samples
LNsumsBCs <- read_table2("R_output/CodeCounts_by_lymphnode.codes", col_names = TRUE)

#make a summary list of all nonredundant samples Inoc, Cow (PT and retro) & LNsums
InocBCs$InocT <- rowSums(InocBCs[2:4], na.rm = TRUE)
Sample_Inoc_LNquads <- full_join(InocBCs, LNpiecesBCs, by = "barcode")
SampleBCs <- full_join(Sample_Inoc_LNquads, LNsumsBCs, by = "barcode")

#write with NA=0 and read the tables back in without NAs
write_tsv(SampleBCs, "R_output/Barcode_count_by_sample_groups.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
SampleBCs <- read_table2("R_output/Barcode_count_by_sample_groups.txt", col_names = TRUE) #Inoc A,B,C,T, PT quads AND retro samples

#### transpose dataframe so rows are samples and columns are BCs
  #Sample of inoc, LN pieces and LN sums (for diversity metrics)
    SampleBCs_t <- as.data.frame(t(SampleBCs[,-1]))
    colnames(SampleBCs_t) <- SampleBCs$barcode

#make a small subset frame to view result
Sample_t_partial <- SampleBCs_t[,1:200]

######-PartA
#make output frame
SampleName <- rownames(SampleBCs_t)

#calculate species number(S), then rarefy species number by sampling 10^6 codes
S <- specnumber(SampleBCs_t)
Srar <- rarefy(SampleBCs_t, sample = 1000000) 
# the three lowest samples are outliers with less than 10^6, most have greater than 10^6 reads, 
# so use 10^6 as sampling size.  Will get an error on undersampled sets, which will return S.


#calculate Shannon H, evenness J,and Fisher alpha, 
H <- diversity(SampleBCs_t, index = "shannon", MARGIN = 1, base = exp(1))
EvenH_J <- H/log(S)
alpha <- fisher.alpha(SampleBCs_t)

#calculate Simpson D, inverse simpson, and evenness (not used in this manuscript)
D <- diversity(SampleBCs_t, index = "simpson")
invD <- diversity(SampleBCs_t, index = "invsimpson")
EvenD <- (invD)*(1/S)

#Make sumary table
Samp_rowSums <- rowSums(SampleBCs_t)
DiversityMetrics <- tibble(SampleName, Samp_rowSums, S, Srar, H, EvenH_J, alpha, D, invD, EvenD )

write_tsv(DiversityMetrics, "R_output/DiversityMetrics.txt", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")

  #plotted Diversity metrics from the above file in prism.


######-PartB
#truncate list to remove rows with LN sums
Sample_t_partial_trim <-Sample_t_partial[1:54,]

SampleBCs_t_trim <- SampleBCs_t[1:54,]

#bray-curtis dis inoc and LN quadrants using counts (not counts normalized with deconstand)
bray <- vegdist(SampleBCs_t_trim,'bray')
newdistB <- as.matrix(bray)
melted_cormatB = melt(newdistB, na.rm = TRUE)

ggheatmapB = ggplot(data = melted_cormatB, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "beige", mid = "red", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", name="Density") +
  theme_minimal()+ 
  coord_fixed()+ 
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_text(angle = 180),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())+
  labs(title="Bray-Curtis Dissimilarity")

ggheatmapB
#export plot as EPS

#jacc dis binary
jaccB <- vegdist(SampleBCs_t_trim,'jaccard', binary = TRUE)
newdistJB <- as.matrix(jaccB)
melted_cormatJB = melt(newdistJB, na.rm = TRUE)

ggheatmapJB = ggplot(data = melted_cormatJB, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "beige", mid = "red", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Density") +
  theme_minimal()+ 
  coord_fixed()+ 
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_text(angle = 180),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank())+labs(title="Jaccard dissimilarity")

ggheatmapJB

# export these two heat maps as .eps and combine in Illustrator.
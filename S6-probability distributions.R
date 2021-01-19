#S7-Probability_distributions
#Binomial distribution plot - probability of detection in LN. Then fit to logistic curve  to assess dose response (ID50)
   
# 1) First, read in inoc and make a "Inoc_ranked_summary"
# 2) then open LN lists and merge inoculum list with fractional abundance
# then bin codes in LN by frac abundance in inoc, and emperically determine the fraction not detected
# in the LN in each bin.  # save these distributions in a new folder
# merge the lists and save a sumary list.
# 3) calculate the theoretical fractions not detected at a range of probabilities using a binomial 
# distribution function.  Save this in a sheet too.
# 4) fit the observed probabilites of loss/detection at different abundancs (doses/number of trials) to a three-parameter 
# logistic curve (dose response curve) to determine the ID50.  


library(tidyverse)
library(drc)
rm(list = ls()) #clear environment

##### 1) Generate list of codes in inoculum with fractional abundance and log(frac abundance) in the inoc pool as a whole

#read in list from each sample
InocA <- read_table2(file="filtered_data/Inoc_x_A.codesCln", col_names = TRUE)
InocB <- read_table2(file="filtered_data/Inoc_x_B.codesCln", col_names = TRUE)
InocC <- read_table2(file="filtered_data/Inoc_x_C.codesCln", col_names = TRUE)

# sum all counts from all three samples
nCountsTot <- sum(c(InocA$counts,InocB$counts,InocC$counts))

#calculate frational abundance of each code:
#sum counts from all inoc then divide by total counts of all inoc.- this way 1s from all samples are treated evenly
Inoc_summary <- bind_rows(InocA, InocB, InocC)
Inoc_summary <- Inoc_summary %>% 
  group_by(barcode) %>% 
  summarise(InocTotCounts = sum(counts), nSamples=n()) %>% 
  arrange(-InocTotCounts) %>%
  mutate ("fracInoc" = InocTotCounts/nCountsTot) 

#add a column with log transformation of  fractional abundance
Inoc_summary <-   mutate(Inoc_summary, "log_frac_inoc" = log10(fracInoc))

#save ranked inoc summary
write_tsv(Inoc_summary, "R_output/Inoc_ranked_summary.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

i <- ggplot(Inoc_summary, aes(x=log_frac_inoc)) + geom_histogram(binwidth=0.1) + scale_y_log10()
i

####### 2) set up thresholds for binning codes (these are the log(fractional abundance) values at the edges of bins)
Bin <- c(-3, -4.2, -4.5, -4.7, -4.8, -4.9, -5.0, -5.1, -5.2, -5.3, -5.4, -5.5, -5.6, -5.7, -5.8, -5.9, -6.0, 
         -6.1, -6.2, -6.35, -6.5, -6.7, -6.9, -7.1, -8.0)

#####Read in list from each LN 

LN_list <- c('LN1L', 'LN2L', 'LN3L', 'LN4L', 'LN5L', 'LN6L', 'LN1R', 'LN2R', 'LN3R', 'LN4R', 'LN5R', 'LN6R')

#for each LN, identify INOC codes detected, then in each inoc abundance bin, count barcodes detected and not detected
# and fill values into a table called LN_Dist.  Then write this LN_Dist file for each of the 12 LNs.

for (LN in LN_list) {
  LNtable <- read_table2(paste("LN_lists/", LN, ".txt", sep=""), col_names = TRUE)
  Inoc_LN <- left_join(Inoc_summary[,c("barcode", "log_frac_inoc")], LNtable[,c("barcode", "median_fracLN")], by = "barcode")
  Inoc_LN <- Inoc_LN %>% replace_na(list(median_fracLN = 0.000000001))  #barcodes not detected in LN represented at 10^-9

  LN_Dist <- tibble(Bin, BinCenter=0, nCodes = 0, nZero = 0, n1 =0, fZero = 0, f1 = 0) #make table to fill in with next loop
  
  for (i in c(2:(length(Bin)))) {
    j <- i-1
    upper <- Bin[j]
    lower <- Bin[i]
    BC_Bin <- Inoc_LN %>% filter(log_frac_inoc >= lower) %>% filter(log_frac_inoc < (upper))
    # fill in table
    LN_Dist[i,2] <- (upper+lower)/2   #calculate BinCenter
    LN_Dist[i,3] <- nrow(BC_Bin) #number of barcodes in bin
    LN_Dist[i,4] <- BC_Bin %>%   filter(median_fracLN == 0.000000001) %>%  summarise(nBC=n())  #number of barcodes in bin not detected
    LN_Dist[i,5] <- BC_Bin %>%   filter(median_fracLN > 0.000000001) %>%  summarise(nBC=n())  #number of barcodes detected
    LN_Dist[i,6] <- LN_Dist[i,4] / nrow(BC_Bin) #fraction of barcodes in bin not detected
    LN_Dist[i,7] <- LN_Dist[i,5] / nrow(BC_Bin) #fraction of barcodes in bin detected
  }
  write_tsv(LN_Dist, paste("LN_lists/", LN, "_distribution.txt", sep = ""), na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")
}

#########merge data from all LN together. 
# Start by reading back in the tables with frac 0 or 1 by bin for each LN

LN1L <- read_table2(file= "LN_lists/LN1L_distribution.txt", col_names = TRUE)
LN2L <- read_table2(file= "LN_lists/LN2L_distribution.txt", col_names = TRUE)
LN3L <- read_table2(file= "LN_lists/LN3L_distribution.txt", col_names = TRUE)
LN4L <- read_table2(file= "LN_lists/LN4L_distribution.txt", col_names = TRUE)
LN5L <- read_table2(file= "LN_lists/LN5L_distribution.txt", col_names = TRUE)
LN6L <- read_table2(file= "LN_lists/LN6L_distribution.txt", col_names = TRUE)

LN1R <- read_table2(file= "LN_lists/LN1R_distribution.txt", col_names = TRUE)
LN2R <- read_table2(file= "LN_lists/LN2R_distribution.txt", col_names = TRUE)
LN3R <- read_table2(file= "LN_lists/LN3R_distribution.txt", col_names = TRUE)
LN4R <- read_table2(file= "LN_lists/LN4R_distribution.txt", col_names = TRUE)
LN5R <- read_table2(file= "LN_lists/LN5R_distribution.txt", col_names = TRUE)
LN6R <- read_table2(file= "LN_lists/LN6R_distribution.txt", col_names = TRUE)


#_make table that sumarizes the frequency of 0 (ie code lost in LN) for each inoculum bin across LNs

Bins_by_LN <- LN1L[,c("BinCenter","nCodes", "fZero")]
Bins_by_LN <- Bins_by_LN %>% rename(LN_1L = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN1R[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_1R = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN2L[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_2L = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN2R[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_2R = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN3L[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_3L = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN3R[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_3R = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN4L[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_4L = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN4R[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_4R = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN5L[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_5L = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN5R[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_5R = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN6L[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_6L = "fZero")
Bins_by_LN <- full_join(Bins_by_LN, LN6R[,c("BinCenter", "fZero")], by = "BinCenter")
Bins_by_LN <- Bins_by_LN %>% rename(LN_6R = "fZero")

Bins_by_LN$meanF0 <- apply(Bins_by_LN[3:14], 1, function(x) mean(x, na.rm = TRUE))
Bins_by_LN$SD_F0 <- apply(Bins_by_LN[3:14], 1, function(x) sd(x, na.rm = TRUE))

Bins_by_LN <- Bins_by_LN[2:24,]  #remove first and last rows which were a placeholders for making bins and represent 0 codes
 
write_tsv(Bins_by_LN, "R_output/FracNotInLN_ByAbundanceBin.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

#transform data from wide to long to make tidy for graphing
Bins_by_LN_tidy <- Bins_by_LN %>% gather(LN_1L:meanF0, key = "LN", value = "fraction0")


LN_F0 <- ggplot(Bins_by_LN_tidy, aes(x=BinCenter, y=fraction0)) + 
  geom_point(aes(color = LN), na.rm= TRUE) + 
  geom_line(aes(color = LN), na.rm= TRUE) + 
  xlim(-7.2, -3) 
LN_F0 

meanLN_F0 <- ggplot(Bins_by_LN, aes(x=BinCenter, y=meanF0, y2=nCodes)) + 
  geom_point() + 
  geom_line() +  xlim(-7.2, -3) 
meanLN_F0

nCodes_bin <- ggplot(Bins_by_LN, aes(x=BinCenter, y=nCodes)) + 
  geom_point() + 
  geom_line() +  xlim(-7.2, -3) + scale_y_log10()
nCodes_bin

###########- 3) determine theoretical fraction of undetected for different abundances with different probabilities 

#this set of codes calculates the frequency (p) of choosing a particular object zero times given s trials and a probability rate of P
#then makes a table where each row represents a particular value of s (number of trials) and each column
#represents P (the probability of success) and each cell if filled with the frequency (p) of success under
#these two conditions using the bionomial distribution

#s = number of trials. In this case approximate number of clones in the inoculum.  For the purpose of modeling, make a vector with a series of possible s.
#f_i = fractional abundance in inoculum = (s / 10^8)

s <- c(1, 2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 32, 40, 50, 63, 79, 100, 126, 158, 200, 251, 316, 398, 501, 631, 794, 1000, 1259, 1585, 1995, 2512, 3162, 3981, 5012, 6310, 7943, 10000, 12589, 15849, 19953, 25119, 31623, 39811, 50119, 63096, 79433, 100000, 200000, 300000, 500000, 750000, 1000000)
position <- c(1:length(s))

P_s <- tibble(s, f_i = 0, P0.01 = 0, P0.001 = 0, P0.0001= 0, P0.00001=0, P0.000001=0)

for (i in position) {
  S <- s[i]  #number of trials
  P_s[i,2] <- S/100000000  # fraction in inoculum (f_i) based on 10^8 CFU
  P_s[i,3] <- dbinom(0, size=S, prob=.01)  #probability of 0 successful events, after S tries, with given proability)
  P_s[i,4] <- dbinom(0, size=S, prob=.001)
  P_s[i,5] <- dbinom(0, size=S, prob=.0001)
  P_s[i,6] <- dbinom(0, size=S, prob=.00001)
  P_s[i,7] <- dbinom(0, size=S, prob=.000001)
  
}

write_tsv(P_s, "R_output/Prob_0_theoretical.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

BN_expect <- ggplot(P_s) + 
  geom_point(aes(x=f_i, y=P0.01, color = 'P0.01')) + 
  geom_point(aes(x=f_i, y=P0.001, color = 'P0.001')) + 
  geom_point(aes(x=f_i, y=P0.0001, color = 'P0.0001')) + 
  geom_point(aes(x=f_i, y=P0.00001, color = 'P0.00001')) + 
  geom_point(aes(x=f_i, y=P0.000001, color = 'P0.000001')) + 
  scale_x_log10()
BN_expect

# I plotted these data together with R_output/FracNotInLN_ByAbundanceBin.txt in prism

#~~~~~~~~~~~ 4) DOSE RESPOPNSE / LOG-LOGISTIC FIT

# import data
lndat <- read_table2(file="R_output/FracNotInLN_ByAbundanceBin.txt", col_names = TRUE)
colnames(lndat) <- c("log.frac", "n.bin", colnames(lndat[,3:16]))
names(lndat) <- make.names(casefold(names(lndat)), allow_=FALSE)

# log.frac =log(fractional abundance)
# n.bin=number of codes in the bin

#  f(x) = d/(1+exp(b(log(x)-log(e))))  [note data already in log space, so in this case f(x) = d/(1+exp(b(x-e))) ]

mL <- drm(meanf0 ~ log.frac, data = lndat, fct = L.3(), type = "continuous") 
summary(mL)
plot(mL, type = "all",log = '')

#Model fitted: Logistic (ED50 as parameter) with lower limit fixed at 0 (3 parms)
#Parameter estimates:

#  Estimate Std. Error  t-value   p-value    
#  b:(Intercept)  5.5082006  0.3010893   18.294 5.882e-14 ***
#  d:(Intercept)  0.9867422  0.0046620  211.655 < 2.2e-16 ***
#  e:(Intercept) -4.4531192  0.0098842 -450.529 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# plot function with fit parameters

x <- c(seq(-8.0, -3.0, 0.1)) 
b <- 5.5082006
d <- 0.9867422
e <- -4.4531192

for (i in x) {
  y= d/(1+exp(b*(x-e)))
  }  

line_fit <- tibble(x,y)
plot(line_fit$x, line_fit$y, type="l")

write_tsv(line_fit, "R_output/DoseResponseFit.txt", na = "0", append = FALSE, col_names = TRUE, quote_escape = "double")

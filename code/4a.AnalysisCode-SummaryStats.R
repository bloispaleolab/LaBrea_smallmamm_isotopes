# install/load packages ----
library(tidyverse)

## Calculate basic statistics for paper ----

# currently within section 'Samples, dates, and stable isotope values' - 1st section of results

tempDat <-read.csv(file="data/processed/final_dataset_alltaxa_finalAMSruns.csv", header=T) 
nrow(tempDat[-which(is.na(tempDat$X14C_age_BP)),]) # 84 samples with dates

tempDat <-read.csv(file="data/processed/final_dataset_alltaxa_finalAMSruns_with_isotopes.csv", header=T) 
nrow(tempDat) # 82 samples with isotopes 

tempDat %>%                               # Summary by group using dplyr
  group_by(Taxon) %>% 
  summarize(minC = min(del13C_permil),
            maxC = max(del13C_permil),
            minN = min(del15N_permil),
            maxN = max(del15N_permil))


tempDat <-read.csv(file="data/processed/final_dataset_focaltaxa_finalAMSruns_with_isotopes_with_all_dates.csv", header=T)
nrow(tempDat) # 76 samples with isotopes from focal taxa 

# isotope stats focal taxa
range(tempDat$del15N_permil) # [1] 3.2 9.9
range(tempDat$del13C_permil) # [1] -22.5 -18.5

# trophic adjustment: 
range(tempDat$del13C_permil - 1)
tempDat[which((tempDat$del13C_permil-1) > -20),]


tempDat <- read.csv(file="data/processed/final_dataset_focaltaxa_dates_isotopes.csv", header=T) #the final set of squirrels and rabbits
#  set of statistics (sample sizes) on final rabbit/squirrel dataset
# currently within section 'Samples, dates, and stable isotope values' - 1st section of results
nrow(tempDat[which(tempDat$Taxon == 'Sylvilagus'),]) # 42
nrow(tempDat[which(tempDat$Taxon == 'Otospermophilus'),]) # 27


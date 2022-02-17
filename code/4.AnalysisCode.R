# install/load packages ----
library(tidyverse)

## Calculate basic statistics for paper ----
# currently within section 'Samples, dates, and stable isotope values' - 1st section of results

mergedDat <-read.csv(file="data/processed/final_dataset_alltaxa_finalAMSruns.csv", header=T) # includes taxa in addition to squirrels, rabbits, excluding the repeated AMS samples.

nrow(mergedDat[-which(is.na(mergedDat$X14C_age_BP)),]) # 84 samples with dates 
nrow(mergedDat[-which(is.na(mergedDat$del15N_permil)),]) # 82 samples with isotopes 
length(which(mergedDat[-which(is.na(mergedDat$del15N_permil)),'Taxon'] == "Sylvilagus")) + length(which(mergedDat[-which(is.na(mergedDat$del15N_permil)),'Taxon'] == "Otospermophilus")) # 76 samples with isotopes from the focal taxa
tempFocalDat <- mergedDat[which(mergedDat$Taxon == 'Sylvilagus' | mergedDat$Taxon == "Otospermophilus"),] # remove non-squirrels or sylvilagus
tempFocalDat <- tempFocalDat[-which(is.na(tempFocalDat$del15N_permil)),] # remove samples without isotopes

# isotope stats focal taxa
range(tempFocalDat$del15N_permil) # [1] 3.2 9.9
range(tempFocalDat$del13C_permil) # [1] -22.5 -18.5

# isotope stats other small mammals
range(mergedDat[c(which(mergedDat$Taxon=="Neotoma"), which(mergedDat$Taxon=="Microtus"), which(mergedDat$Taxon=="Thomomys")),'del15N_permil'], na.rm=T)
range(mergedDat[c(which(mergedDat$Taxon=="Neotoma"), which(mergedDat$Taxon=="Microtus"), which(mergedDat$Taxon=="Thomomys")),'del13C_permil'], na.rm=T)

# isotope stats carnivores
range(mergedDat[c(which(mergedDat$Taxon=="Canis"), which(mergedDat$Taxon=="Mustela")),'del15N_permil'])
range(mergedDat[c(which(mergedDat$Taxon=="Canis"), which(mergedDat$Taxon=="Mustela")),'del13C_permil'])

trimmedDat <- read.csv(file="data/processed/final_dataset_focaltaxa.csv", header=T) #the final set of squirrels and rabbits
#  set of statistics (sample sizes) on final rabbit/squirrel dataset
# currently within section 'Samples, dates, and stable isotope values' - 1st section of results
nrow(trimmedDat[which(trimmedDat$Taxon == 'Sylvilagus'),]) # 42
nrow(trimmedDat[which(trimmedDat$Taxon == 'Otospermophilus'),]) # 27


# read in final dataset ----
matchedDF_all <- read.csv(file="data/processed/final_dataset_focaltaxa.csv", header=T)

# Intra- and inter-specific isotopic niches ----

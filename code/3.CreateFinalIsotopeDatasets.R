# Create the base isotope and date data files ----

## Install relevant packages ----
library(tidyverse)

## Section 1: merge dates and isotopes together ----
### Read in original data from google and arrange ----
datesOrig <- read.delim(file="data/processed/master_dates_file.txt", sep="\t", header=T)
isoOrig <- read.delim(file="data/processed/master_isotopes_file.txt", sep="\t", header=T)

isoSub <- isoOrig %>% select(UCIAMS_Number, Museum_Number, del15N_permil, del13C_permil, prelim_taxon_name)
datesSub <- datesOrig %>% select(UCIAMS_Number, Museum_Number, X14C_age_BP, X14C_age_error, prelim_taxon_name)

# find unique museum numbers and lab numbers
specimens <- unique(union(isoSub$Museum_Number, datesSub$Museum_Number))
samples <- unique(union(isoSub$UCIAMS_Number, datesSub$UCIAMS_Number))

# create shell of the final merged data frame
mergedDat <- as.data.frame(matrix(NA, nrow=0, ncol=7))
colnames(mergedDat)<- c("UCIAMS_Number", "Museum_Number", "prelim_taxon_name", "X14C_age_BP", "X14C_age_error", "del15N_permil",  "del13C_permil")

### merge non-repeated samples ----
# create dataframe of non-repeated samples first - specimens with one isotope value and one or no dates 
for (i in 1:length(specimens)){
  id <- specimens[i]
  IsoSamp <- (which(isoSub$Museum_Number == id))
  DatesSamp <- (which(datesSub$Museum_Number == id))
  
  if (length(IsoSamp) == 1 && length(DatesSamp) == 1){ # samples with dates and isotopes
    mergedDat <- rbind(mergedDat, cbind(datesSub[DatesSamp, c("UCIAMS_Number", "Museum_Number", "prelim_taxon_name", "X14C_age_BP", "X14C_age_error")], isoSub[IsoSamp, c("del15N_permil",  "del13C_permil")]))
  }
  if (length(IsoSamp) == 1 && length(DatesSamp) == 0){ # samples with isotopes and no dates
    tempDat <- cbind(isoSub[IsoSamp, c("UCIAMS_Number", "Museum_Number", "prelim_taxon_name")], NA, NA,
                     isoSub[IsoSamp, c("del15N_permil",  "del13C_permil")])
    colnames(tempDat)[4:5] <- c("X14C_age_BP", "X14C_age_error")
    mergedDat <- rbind(mergedDat, tempDat)
  }
}

mergedDat$repeated <- as.vector(rep("N", nrow(mergedDat)))

### figure out which samples had repeated dates or isotopes ----
leftovers <- specimens[which(is.na(match(specimens, mergedDat$Museum_Number)))]
leftoverDF <- as.data.frame(leftovers)
leftoverDF$NumIso <- NULL
leftoverDF$NumDates <- NULL

for (j in 1:length(leftovers)){
  id <- leftovers[j]
  IsoSamp <- (which(isoSub$Museum_Number == id))
  DatesSamp <- (which(datesSub$Museum_Number == id))
  
  leftoverDF$NumIso[j] <- length(IsoSamp)
  leftoverDF$DatesIso[j] <- length(DatesSamp)
}

# which UCIAMS to use?
# LACMP23-33228: 198302 #alt 198206 use 198206 because the C:N ratio is  higher

# LACMHC-142773: 216768 #alt 217076 use 217076 because the C:N ratio is slightly higher (3.0, vs 2.9 for 216768)

# LACMHC-142779: 216770 #alt 217077 use 217077 because the C:N ratio is slightly higher (3.0, vs 2.9 for 216770)

# LACMP23-35541: 223585 #alt 223495  ## use 223495
# note: 223495 used in SIBERraw.csv, 223585 used in Iso.Clim.R script. Sample 223495 has a higher C:N (2.9) vs 2.8 for 223585. We will use sample with higher C:N (223495).

# LACMP23-40642: 223587 #alt 223521  ## use 223587
# note: 223521 used in SIBERraw.csv, 223587 used in Iso.Clim.R script. The two samples have same C:N, isotope values, slightly different dates (26940 +- 320 vs 26710 +- 110). We will use sample with more precise age range (223587). 

leftoverDF$correctUCIAMS<- c(198206, 217076, 217077, 223495, 223587)
leftoverDF$repeatedUCIAMS <- c(198302, 216768, 216770, 223585, 223521)

## CHECK!
leftoverDF[,c(1,4)]
# Should have in correctUCIAMS:
# 1 LACMP23-33228        198206
# 2 LACMHC-142773        217076
# 3 LACMHC-142779        217077
# 4 LACMP23-35541        223495
# 5 LACMP23-40642        223587
  
### merge repeated samples ----
# merge correct and repeated samples in turn
for (i in 1:nrow(leftoverDF)){
  
  correct_ams <- leftoverDF$correctUCIAMS[i]
  repeated_ams <- leftoverDF$repeatedUCIAMS[i]
  
  IsoSamp <- c(which(isoSub$UCIAMS_Number == correct_ams), which(isoSub$UCIAMS_Number == repeated_ams))
  DatesSamp <- c(which(datesSub$UCIAMS_Number == correct_ams), which(datesSub$UCIAMS_Number == repeated_ams))
  
  repeated <- as.vector(c("N", "Y"))
  
  mergedDat <- rbind(mergedDat, cbind(datesSub[DatesSamp, c("UCIAMS_Number", "Museum_Number", "prelim_taxon_name", "X14C_age_BP", "X14C_age_error")], isoSub[IsoSamp, c("del15N_permil",  "del13C_permil")], repeated))
  
}

### create the Species and Taxon columns ----
Species <- vector(length=nrow(mergedDat))
Species[grep("audu", mergedDat$prelim_taxon_name)] <- 'S. audubonii'
Species[grep("bach", mergedDat$prelim_taxon_name)] <- 'S. bachmani'
Species[grep("Otosperm", mergedDat$prelim_taxon_name)] <- 'Otospermophilus'
Species[grep("Neotoma", mergedDat$prelim_taxon_name)] <- 'Neotoma sp'
Species[grep("Thomomys", mergedDat$prelim_taxon_name)] <- 'Thomomys sp'
Species[grep("Microtus", mergedDat$prelim_taxon_name)] <- 'Microtus sp'
Species[grep("Mustela", mergedDat$prelim_taxon_name)] <- 'Mustela frenata'
Species[grep("Canis", mergedDat$prelim_taxon_name)] <- 'Canis latrans'
Species[grep("Mammalia-Lagomorpha", mergedDat$prelim_taxon_name)] <- 'Leporidae' #Sylvilagus
Species[which(Species=="FALSE")] <- 'Sylvilagus sp'

Taxon <- vector(length=nrow(mergedDat))
Taxon[grep("audu", mergedDat$prelim_taxon_name)] <- 'Sylvilagus'
Taxon[grep("bach", mergedDat$prelim_taxon_name)] <- 'Sylvilagus'
Taxon[grep("Otosperm", mergedDat$prelim_taxon_name)] <- 'Otospermophilus'
Taxon[grep("Neotoma", mergedDat$prelim_taxon_name)] <- 'Neotoma'
Taxon[grep("Thomomys", mergedDat$prelim_taxon_name)] <- 'Thomomys'
Taxon[grep("Microtus", mergedDat$prelim_taxon_name)] <- 'Microtus'
Taxon[grep("Mustela", mergedDat$prelim_taxon_name)] <- 'Mustela'
Taxon[grep("Canis", mergedDat$prelim_taxon_name)] <- 'Canis'
Taxon[grep("Mammalia-Lagomorpha", mergedDat$prelim_taxon_name)] <- 'Sylvilagus'
Taxon[which(Species=="Sylvilagus sp")] <- 'Sylvilagus'

mergedDat$Species <- Species
mergedDat$Taxon <- Taxon

# filter and write various datasets ----

# save the most expanded dataset
# includes taxa in addition to squirrels, rabbits, as well as the repeated AMS samples.
write.csv(mergedDat, file="data/processed/final_dataset_alldata.csv", row.names=F) 

# filter and save final datasets for different analyses
# remove the repeated AMS samples
mergedDat <- mergedDat[which(mergedDat$repeated == "N"),]
write.csv(mergedDat, file="data/processed/final_dataset_alltaxa_finalAMSruns.csv", row.names=F) 

# remove samples without isotopes
trimmedDat <- mergedDat[-which(is.na(mergedDat$del15N_permil)),] 
write.csv(trimmedDat, file="data/processed/final_dataset_alltaxa_finalAMSruns_with_isotopes.csv", row.names=F) 

# remove non-squirrels or sylvilagus
trimmedDat <- trimmedDat[which(trimmedDat$Taxon == 'Sylvilagus' | trimmedDat$Taxon == "Otospermophilus"),] 
write.csv(trimmedDat, file="data/processed/final_dataset_focaltaxa_finalAMSruns_with_isotopes_with_all_dates.csv", row.names=F) 

# further trim down dataset for primary analyses
# remove samples without dates 
trimmedDat <- trimmedDat[-which(is.na(trimmedDat$X14C_age_BP)),] 
# remove samples with too old of dates
trimmedDat <- trimmedDat[-which(is.na(trimmedDat$X14C_age_error)),] 

#this is the final set of squirrels and rabbits
write.csv(trimmedDat, file="data/processed/final_dataset_focaltaxa_dates_isotopes.csv", row.names = F) 

## Section 2: Match calibrated ages to samples ----

### create sample names for matching isotope file to ages file ----
samples <- paste("UCIAMS", trimmedDat$UCIAMS_Number) # this is the final set of samples with both dates and isotope data

### read in calibrated ages ----
allAges<- read.csv('data/OxCal/final oxcal models/AllAges_ages_probs.csv', header=T) # this file stores the raw probabilities for each age across the distribution
all_calibrated_ages <- read.csv('data/OxCal/final oxcal models/AllAges_forinput.csv', header=T) # this file stores the calibrated age statistics for each age
all_calibrated_ages$trimmedName <- unlist(lapply(strsplit(all_calibrated_ages$Name, " R_"), '[[', 1))
sample_median_ages <- as.data.frame(cbind(samples, all_calibrated_ages[match(samples, all_calibrated_ages$trimmedName), 'Unmodelled._BP_median']))
colnames(sample_median_ages)[2] <- 'median_age'  
sample_median_ages$median_age <- as.numeric(sample_median_ages$median_age)

### data cleaning on allAges file----
allAges<- allAges[allAges$probability != 0, ] # remove age estimates with 0 probability

#remove repetitive information and duplicate samples from allAges file
sampsToDelete <- unique(allAges$name)[which(is.na(match(unique(allAges$name), samples)))]
for (i in 1:length(sampsToDelete)){
  allAges<- allAges[-which(allAges$name==sampsToDelete[i]),]
}

## Section 3: Match d180 to all age estimates ----
# (both median_age as well as the full distribution of ages) for every sample 

if (all(allAges$value<0)==T){ # this should be true, then take the absolute values to match to climate
  xout <- abs(allAges$value)
}

# first match hendy data to age data
hendyDat<- read.delim("data/raw/climate/hendy2002data.txt")

# adding d18O to allAges to match to the full distribution of age estimates for each sample
Hendy_extracted <- approx(x=hendyDat$HendyAge, y=hendyDat$pach.d18O, method="linear", xout=xout)
d18O_hendy <-Hendy_extracted$y
allAges_d18O<-cbind(allAges, d18O_hendy) # this matches d180 to all the age estimates, for all probabilities

# matching d18O to the median ages
d18O_hendy_medianage <- approx(x=hendyDat$HendyAge, y=hendyDat$pach.d18O, 
                               method="linear", xout=sample_median_ages$median_age)$y

# next match ngrip data to age data
ngripDat<- read.delim("data/raw/climate/ngrip.txt")
ngripDat$Age <- ngripDat$Age-50 #ngrip dates are "before 2000", so need to subtract 50 years

# adding d18O to allAges to match to the full distribution of age estimates for each sample
ngrip_extracted <- approx(x=ngripDat$Age, y=ngripDat$d18O, method="linear", xout=xout)
d18O_ngrip <-ngrip_extracted$y
allAges_d18O <- cbind(allAges_d18O, d18O_ngrip) # this matches d180 to all the age estimates, for all probabilities

# matching d18O to the median ages
d18O_ngrip_medianage <- approx(x=ngripDat$Age, y=ngripDat$d18O, 
                               method="linear", xout=sample_median_ages$median_age)$y

# clean up matched d18O data
# d18O_hendy # not needed - bound to allAges_d18O
# d18O_ngrip # not needed - bound to allAges_d18O
rm(d18O_hendy)
rm(d18O_ngrip)

# d18O_hendy_medianage - matched to median age within sample_median_ages
# d18O_ngrip_medianage - matched to median age within sample_median_ages
sample_median_ages <- cbind(sample_median_ages, d18O_hendy_medianage, d18O_ngrip_medianage)
colnames(sample_median_ages)[3:4] <- c("d18O_hendy", "d18O_ngrip")
rm(d18O_hendy_medianage)
rm(d18O_ngrip_medianage)

### at end of Section 3, have 2 main files for later use: ###
# allAges_d18O (stores both ngrip and hendy)
# sample_median_ages (stores both ngrip and hendy)

## Section 4: additional data cleaning ----
# remove any unmatched specimen from the files
matches <- match(samples, allAges_d18O$name)

# added ifelse statement to catch instances where there are or are not matches
if (any(is.na(matches))){
  isoDat2_final <- trimmedDat[-which(is.na(matches)),]
  samples_final <- samples[-which(is.na(matches))]
}else{
  isoDat2_final <- trimmedDat
  samples_final <- samples
}
# end result dataframes:
# isoDat2_final - final C/N isotope data
# samples_final - final set of samples
# allAges_d18O - ages, probabilities, and d18O for the full distribution of calibrated age estimates associated with each sample
# sample_median_ages - median ages and d180 estimates at those median ages for each sample

## Section 5: Age uncertainty: calculate weighted age and weighted mean d18O for each specimen, add time gorup ----
# For each specimen, calculate the weighted d180 and weighted mean age, based on the probability of each calibrated age from the age models

specimen_wage <- vector(mode="numeric", length=length(samples_final))
specimen_wd18O_hendy <- vector(mode="numeric", length=length(samples_final))
specimen_wd18O_ngrip <- vector(mode="numeric", length=length(samples_final))

for (i in 1:length(samples_final)){
  # pull out specimen ID
  spec <- samples_final[i] 
  # pull out all calibrated ages and probabilities for that single specimen
  specd18O <- allAges_d18O[which(allAges_d18O$name == spec),]
  # calculate weighted mean d18O
  specimen_wd18O_hendy[i] <- weighted.mean(specd18O$d18O_hendy, specd18O$probability)
  specimen_wd18O_ngrip[i] <- weighted.mean(specd18O$d18O_ngrip, specd18O$probability)
  # calculate weighted mean age
  specimen_wage[i] <- weighted.mean(specd18O$value, specd18O$probability)
}

# match all age and d18O estimates with the original C&N isotope data
matchedDF_weighted <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Species', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP', 'X14C_age_error')], specimen_wage, specimen_wd18O_hendy, specimen_wd18O_ngrip)

# make sure still in same order (ie, no nas, numbers increase in sequence)
match(gsub("UCIAMS ", "", sample_median_ages$samples), matchedDF_weighted$UCIAMS_Number)

matchedDF_all <- cbind(matchedDF_weighted, sample_median_ages[,2:4]) 

# add time_group onto the dataframe
time_group <- vector(length=nrow(matchedDF_all))
time_group[which(matchedDF_all$median_age <= 11700)] <- "Holocene"
time_group[which(matchedDF_all$median_age > 11700)] <- "Pleistocene"

matchedDF_all$time_group <- time_group

# Notes about the final dataset! ----
# matchedDF_all thus is the full final dataset
# UCIAMS_Number
# Species
# Taxon
# del15N_permil      
# del13C_permil
# X14C_age_BP
# X14C_age_error
# specimen_wage # sensitivity age
# specimen_wd18O_hendy # sensitivity d18O estimate - hendy
# specimen_wd18O_ngrip # sensitivity d18O estimate - ngrip
# median_age # primary age
# d18O_hendy #primary d18O estimate at median age - hendy
# d18O_ngrip #primary d18O estimate at median age - hendy
# time_group # pleistocene or holocene

write.csv(matchedDF_all, file="data/processed/final_dataset_focaltaxa_with_calages_climate.csv", row.names = F) #just the final set of squirrels and rabbits




# this should work for anything that uses SIBERraw.
# for the analyses in the iso.clim.R script, need to reconcile the column names. 
# Here are the original column names
# 
# iso_dat:
# UCIAMS_Number	 ok, same
# Original_sample_name	# deleted, not used in iso.clim
# Species	ok, same
# Taxon	 ok, same
# d15N	
# d13C	
# Calibrated_mean_age	
# Sigma	
# Cal_median_age	
# X14C_age_BP	 ok, same
# X14C_age_error ok, same	
# pach.d18O_mean	
# pach.d18O_median	
# P3_L

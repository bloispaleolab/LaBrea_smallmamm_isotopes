# Sensitivity analysis - age sensitivity ----

# Read in data ----
matchedDF_all <- read.csv(file="data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T) # 

# sensitivity analysis #1 ----
# compare weighted and median age models, using the final models from the primary analysis

C_model_weighted <- lm(del13C_permil ~ specimen_wd18O_hendy + Taxon, data=matchedDF_all)
summary(C_model_weighted)
C_model_median <- lm(del13C_permil ~ d18O_hendy + Taxon, data=matchedDF_all)
summary(C_model_median) # this is the final model from the primary analyses

N_model_weighted <- lm(del15N_permil ~ specimen_wd18O_hendy + Taxon, data=matchedDF_all)
summary(N_model_weighted)
N_model_median <- lm(del15N_permil ~ d18O_hendy + Taxon, data=matchedDF_all)
summary(N_model_median) # this is the final model from the primary analyses

# Results: compare weighted d180 vs d18O at median age
# No substantial difference. N still marginlly signif, C still highly signif. 
# Median is more significant than weighted, for what that's worth. Maybe only matters when we go to match isotopes with a specific age? weighted age and median age can be quite different


# sensitivity analysis #2 ----
# How much of a difference does the variation in age make?
# for this, need to go back to original code and a dataset earlier in the cleaning process to match samples to climate and sample climate probabilistically across the full age spectrum of each specimen.

# reconstruct datasets (copied over relevant code from 3.CreateFinalIsotopeDatasets.R)
#this is the final set of squirrels and rabbits
trimmedDat <- read.csv(file="data/processed/final_dataset_focaltaxa_dates_isotopes.csv", header=T)

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

## Section 3: Match d180 to all age estimates - just for hendy ----
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

# clean up matched d18O data
# d18O_hendy # not needed - bound to allAges_d18O
rm(d18O_hendy)

# d18O_hendy_medianage - matched to median age within sample_median_ages
sample_median_ages <- cbind(sample_median_ages, d18O_hendy_medianage)
colnames(sample_median_ages)[3] <- c("d18O_hendy")
rm(d18O_hendy_medianage)

### at end of Section 3, have 2 main files for later use: ###
# allAges_d18O (stores just hendy)
# sample_median_ages (stores just hendy)

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

# ok, now that the appropriate datasets have been recreated, do the sensitivity analysis
N=100 # Note: some specimens do not have 100 age estimates.
N_model_res <- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
C_model_res<- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
colnames(N_model_res) <- colnames(C_model_res) <- c('Fstat', 'lm_pVal', 'coeff', 'AdjR2', 'cor', 'cor_pVal')

for (k in 1:N){
  
  # for each specimen, sample a d18O value
  sampledd18O <- vector(mode="numeric", length=length(samples_final))
  
  for (i in 1:length(samples_final)){
    # pull out specimen ID
    spec <- samples_final[i] 
    # pull out all calibrated ages and probabilities
    specd18O <- allAges_d18O[which(allAges_d18O$name == spec),]
    # sample a single age, with sampling weighted per the probability
    sampledd18O[i] <- sample(specd18O$d18O, 1, prob=specd18O$probability)
  }
  
  # match with isotope data
  matchedDF_sensitivity <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], sampledd18O)
  
  N_model_sensitivity <- lm(del15N_permil ~ sampledd18O, data=matchedDF_sensitivity)
  C_model_sensitivity <- lm(del13C_permil ~ sampledd18O, data=matchedDF_sensitivity)
  
  par(mfrow=c(1,2))
  plot(del15N_permil ~ sampledd18O, data=matchedDF_sensitivity, pch=16)
  abline(N_model_sensitivity)
  plot(del13C_permil ~ sampledd18O, data=matchedDF_sensitivity, pch=16)
  abline(C_model_sensitivity)
  
  N_cor.test_sensitivity <- cor.test(matchedDF_sensitivity$del15N_permil, matchedDF_sensitivity$sampledd18O)
  C_cor.test_sensitivity <- cor.test(matchedDF_sensitivity$del13C_permil, matchedDF_sensitivity$sampledd18O)
  
  N_model_res[k, 1]<- summary(N_model_sensitivity)$fstatistic[1]
  N_model_res[k, 2]<- summary(N_model_sensitivity)$coefficients[8]
  N_model_res[k, 3]<- summary(N_model_sensitivity)$coefficients[2]
  N_model_res[k, 4]<- summary(N_model_sensitivity)$adj.r.squared
  N_model_res[k, 5]<- N_cor.test_sensitivity$estimate
  N_model_res[k, 6]<- N_cor.test_sensitivity$p.value
  
  C_model_res[k, 1]<- summary(C_model_sensitivity)$fstatistic[1]
  C_model_res[k, 2]<- summary(C_model_sensitivity)$coefficients[8]
  C_model_res[k, 3]<- summary(C_model_sensitivity)$coefficients[2]
  C_model_res[k, 4]<- summary(C_model_sensitivity)$adj.r.squared
  C_model_res[k, 5]<- C_cor.test_sensitivity$estimate
  C_model_res[k, 6]<- C_cor.test_sensitivity$p.value
  
}

# Calculate summary statistics on each column
apply(N_model_res, 2, summary)
apply(C_model_res, 2, summary)

# Supplemental figure for appendix - sensitivity analysis ----
## Compare final estimate with sensitivity analysis ----

# new supplemental figure, showing the test results for both Adj. R2 and the fitted coefficient for d18O.
# Note the following models to use:
# Carbon median ages univariate climate: carbon.lm.clim
# Nitrogen median ages univariate climate: nitrogen.lm.clim

# Carbon weighted age univariate climate
carbon.lm.clim.weighted <- lm(del13C_permil ~ specimen_wd18O_hendy, data=matchedDF_all)
# Nitrogen weighted age univariate climate
nitrogen.lm.clim.weighted <- lm(del15N_permil ~ specimen_wd18O_hendy, data=matchedDF_all)

# original models
carbon.lm.clim<-lm(del13C_permil~d18O_hendy, matchedDF_all)
summary(carbon.lm.clim)

nitrogen.lm.clim<-lm(del15N_permil ~ d18O_hendy, data=matchedDF_all)
summary(nitrogen.lm.clim)

# grDevices::pdf("output/SuppFigure1_climate_sensitivity_R2andcoeff.pdf", width=8, height=6)
pdf(file="output/SuppFigure1_climate_sensitivity_R2andcoeff-NF.pdf", height=6, width=8)
#grDevices::cairo_pdf("output/SuppFigure1_climate_sensitivity_R2andcoeff.pdf", width=8, height=6)

par(mfrow=c(2,2)) 

# plot Carbon R2
y1<- hist(C_model_res$AdjR2, xlab=expression('Adjusted R'^2), main=expression('Model:'~{delta}^13*C~'\u2030'~'~'~{delta}^18*O~'\u2030'))
segments(summary(carbon.lm.clim.weighted)$adj.r.squared, 0, 
         summary(carbon.lm.clim.weighted)$adj.r.squared, max(y1$counts), 
         col="red", lwd=1)
segments(summary(carbon.lm.clim)$adj.r.squared, 0, 
         summary(carbon.lm.clim)$adj.r.squared, max(y1$counts), 
         col="blue", lwd=2)
legend(xpd=T, x=0.07, y=24, bty="n", legend=c(expression('median'~{delta}^18*O), expression('weighted'~{delta}^18*O)), 
       col=c("blue", "red"), lwd=c(2,1), cex=0.75) # change x and y to fit with panel locations


# plot Carbon d18O fitted coefficient
y2 <- hist(C_model_res$coeff, xlab=expression({delta}^18*O~'fitted coefficient'), main="")
segments(summary(carbon.lm.clim.weighted)$coefficients[2], 0, 
         summary(carbon.lm.clim.weighted)$coefficients[2], max(y2$counts), 
         col="red", lwd=1)
segments(summary(carbon.lm.clim)$coefficients[2], 0, 
         summary(carbon.lm.clim)$coefficients[2], max(y2$counts), 
         col="blue", lwd=2)


# plot Nitrogen R2
x1<- hist(N_model_res$AdjR2, xlab=expression('Adjusted R'^2), main=expression('Model:'~{delta}^15*N~'\u2030'~'~'~{delta}^18*O~'\u2030'))
segments(summary(nitrogen.lm.clim.weighted)$adj.r.squared, 0, 
         summary(nitrogen.lm.clim.weighted)$adj.r.squared, max(x1$counts), 
         col="red", lwd=1)
segments(summary(nitrogen.lm.clim)$adj.r.squared, 0, 
         summary(nitrogen.lm.clim)$adj.r.squared, max(x1$counts), 
         col="blue", lwd=2)

# plot Nitrogen d18O fitted coefficient
x2 <- hist(N_model_res$coeff, xlab=expression({delta}^18*O~'fitted coefficient'), main="")
segments(summary(nitrogen.lm.clim.weighted)$coefficients[2], 0, 
         summary(nitrogen.lm.clim.weighted)$coefficients[2], max(x2$counts), 
         col="red", lwd=1)
segments(summary(nitrogen.lm.clim)$coefficients[2], 0, 
         summary(nitrogen.lm.clim)$coefficients[2], max(x2$counts), 
         col="blue", lwd=2)
dev.off()


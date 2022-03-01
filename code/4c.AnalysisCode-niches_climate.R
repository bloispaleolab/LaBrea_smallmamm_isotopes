# Set up R session ----
library(dplyr)
library(stringr)
library(ggplot2)

# random function so I can extract p-values from lm
Regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Read in data ----
matchedDF_all <- read.csv(file="data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T) # 


######### Section 5: Final Models and Figures ##############
# using the final dataframe: matchedDF_all
# we should use the median d18O values for final analysis, and present the weighted ones in supplemental, along with sensitivity analysis

# Methods potential text:
# For carbon and nitrogen, we fit a linear model that originally included oxygen, taxon, and the interaction between the two variables as independent variables. We then performed stepwise regression to determine a final model. 
# Results potential text
# Stepwise regression indicated that there was no significant interaction between 18O and taxon for either carbon or nitrogen stable isotope values. For carbon, variation in 13C was significantly associated with both 18O and taxon (stats from summary(carbon.lm.final)). For nitrogen, neither 18O nor taxon explained significant variation in 15N, though taxon as a single variable was marginally significant (stats from summary(nitrogen.lm.taxon)).

final.model.stats <- as.data.frame(matrix(data=NA, nrow=4, ncol=5))
colnames(final.model.stats) <- c("element", "climate", "final_model", "p", "adjR2")

## start with Hendy climate data ----
### models - Carbon ----

# Start simple and build complexity

# climate-only model
carbon.lm.clim<-lm(del13C_permil~d18O_hendy, matchedDF_all)
summary(carbon.lm.clim)
plot(del13C_permil ~ d18O_hendy, data=matchedDF_all, pch=16)
abline(carbon.lm.clim)

# sensitivity test: taxon specific climate-only models (both are still significant)
# Sylvilagus have a stronger relationship with climate than Otospermophilus

carbon.lm.clim.oto <-lm(del13C_permil~d18O_hendy, data=matchedDF_all[which(matchedDF_all$Taxon=="Otospermophilus"),])
summary(carbon.lm.clim.oto)

carbon.lm.clim.syl<-lm(del13C_permil~d18O_hendy, data=matchedDF_all[which(matchedDF_all$Taxon=="Sylvilagus"),])
summary(carbon.lm.clim.syl)

# sensitivity test: Pleistocene models. Not significant!!
# just climate, pleistocene (older than 11500)
carbon.lm.clim.ple <-lm(del13C_permil~d18O_hendy, data=matchedDF_all[which(matchedDF_all$X14C_age_BP>11500),])
summary(carbon.lm.clim.ple)
plot(del13C_permil ~ d18O_hendy, data=matchedDF_all[which(matchedDF_all$X14C_age_BP>11500),], pch=16)
abline(carbon.lm.clim.ple)

# climate + taxon, just pleistocene. Taxon signficant, climate not
carbon.lm.all.additive.ple <-lm(del13C_permil~d18O_hendy + Taxon, data=matchedDF_all[which(matchedDF_all$X14C_age_BP>11500),])
summary(carbon.lm.all.additive.ple)

# otosperm climate only, pleistocene only model. climate not significant
carbon.lm.clim.ple.oto <-lm(del13C_permil~d18O_hendy, data=matchedDF_all[intersect(which(matchedDF_all$X14C_age_BP>11500), which(matchedDF_all$Taxon == "Otospermophilus")),])
summary(carbon.lm.clim.ple.oto)
plot(del13C_permil ~ d18O_hendy, data=matchedDF_all[intersect(which(matchedDF_all$X14C_age_BP>11500), which(matchedDF_all$Taxon == "Otospermophilus")),], pch=16, col="royalblue2")
abline(carbon.lm.clim.ple.oto)

# sylvilagus climate only, pleistocene only model. climate not significant
carbon.lm.clim.ple.syl <-lm(del13C_permil~d18O_hendy, data=matchedDF_all[intersect(which(matchedDF_all$X14C_age_BP>11500), which(matchedDF_all$Taxon == "Sylvilagus")),])
summary(carbon.lm.clim.ple.syl)
plot(del13C_permil ~ d18O_hendy, data=matchedDF_all[intersect(which(matchedDF_all$X14C_age_BP>11500), which(matchedDF_all$Taxon == "Sylvilagus")),], pch=16, col="darkorange")
abline(carbon.lm.clim.ple.syl)


# original models: taxon-only model
carbon.lm.taxon<-lm(del13C_permil~Taxon, data=matchedDF_all)
summary(carbon.lm.taxon)

# model with no interaction term 
# NOTE: this is what the final model is in the end, so should be the same as carbon.lm.final
carbon.lm.all.additive<-lm(del13C_permil~d18O_hendy + Taxon, data=matchedDF_all)
summary(carbon.lm.all.additive)

# model with interaction term included
carbon.lm.all.interaction<-lm(del13C_permil~d18O_hendy * Taxon, 
                              data=matchedDF_all)
summary(carbon.lm.all.interaction)

##NOTE - no significant interactions, total amount of variation explained is lower than with additive model. And, when interaction term is included, climate emerges as only significant variable. generally indicating that there still is a significant amount of variation explained by taxon, but that adding the interaction spreads it out among two variables, neither of which are significant.

#stepwise regression 
carbon.lm.final <- step(lm(del13C_permil~d18O_hendy * Taxon, data=matchedDF_all), direction="both")
summary(carbon.lm.final)
# --> this shows that the d180 + Taxon model (no interaction) is the best final model. 

# t-test of residuals from climate only model
# residuals of Otospermophilus and Sylvilagus are signif different, again indicating effect of climate
c.t.clim <- t.test(carbon.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))])
c.t.clim 

### models - Nitrogen ----

# Start simple and build complexity

# climate-only model
nitrogen.lm.clim<-lm(del15N_permil ~ d18O_hendy, data=matchedDF_all)
summary(nitrogen.lm.clim)

# taxon-only model
nitrogen.lm.taxon<-lm(del15N_permil~Taxon, data=matchedDF_all)
summary(nitrogen.lm.taxon)

# model with no interaction term 
nitrogen.lm.all.additive<-lm(del15N_permil ~ d18O_hendy + Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.additive)

# model with interaction term included
nitrogen.lm.all.interaction<-lm(del15N_permil ~ d18O_hendy * Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.interaction)

#stepwise regression --> this shows that there is not a great final model
nitrogen.lm.final <- step(lm(del15N_permil~d18O_hendy*Taxon, data=matchedDF_all), direction="both")
summary(nitrogen.lm.final)
# --> the best final model is the additive model, which is the same as carbon.

# t-test of residuals from climate-only model
# residuals of Otospermophilus and Sylvilagus are NOT signif different
n.t.clim <- t.test(nitrogen.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))])
n.t.clim


### Do we need to add age as a co-variate in the models? ----

# Does this answer the question: "does isotope vary through time?"
carbon.lm.all.age<-step(lm(del13C_permil~d18O_hendy * Taxon * median_age, data=matchedDF_all), direction="both")
summary(carbon.lm.all.age)
# explain a bit more variation with age included: 0.4877 vs 0.428

nitrogen.lm.all.age<-step(lm(del15N_permil ~ d18O_hendy * Taxon * median_age, data=matchedDF_all), direction="both")
summary(nitrogen.lm.all.age)
# explain a bit more variation with age included: 0.1702 vs 0.08476

final.model.stats[1,] <- c("carbon", "hendy", "carbon.lm.all.additive", 
                           Regressionp(carbon.lm.all.additive), summary(carbon.lm.all.additive)$adj)
final.model.stats[2,] <- c("nitrogen", "hendy", "nitrogen.lm.all.additive", 
                           Regressionp(nitrogen.lm.all.additive), summary(nitrogen.lm.all.additive)$adj)

# rename final models before changing climate data
carbon.lm.final.hendy <- carbon.lm.final
carbon.lm.clim.hendy <- carbon.lm.clim
carbon.lm.all.age.hendy <- carbon.lm.all.age
c.t.clim.hendy <- c.t.clim
nitrogen.lm.final.hendy <- nitrogen.lm.final
nitrogen.lm.clim.hendy <- nitrogen.lm.clim
nitrogen.lm.all.age.hendy <- mitrogen.lm.all.age
n.t.clim.hendy <- n.t.clim

# Save model objects for use in plotting
save(carbon.lm.final.hendy, carbon.lm.clim.hendy, carbon.lm.all.age.hendy, c.t.clim.hendy, 
     nitrogen.lm.final.hendy, nitrogen.lm.clim.hendy, nitrogen.lm.all.age.hendy, n.t.clim.hendy, 
     file="output/hendy_models_for_plotting.RData")


# remove variables that are being generated anew in ngrip section
rm("c.t.clim", "carbon.lm.all.additive", "carbon.lm.all.age", "carbon.lm.all.interaction", "carbon.lm.clim", "carbon.lm.final", "carbon.lm.taxon",
   "n.t.clim", "nitrogen.lm.all.additive", "nitrogen.lm.all.age", "nitrogen.lm.all.interaction", "nitrogen.lm.clim", "nitrogen.lm.final", "nitrogen.lm.taxon")

# Sensitivity analysis - Climate ----

## ngrip climate data ----
# Note that the ngrip data are negatively correlated with hendy data, so I have converted them to negative values so that the relationships sync up with the hendy dataset.
matchedDF_all$d18O_ngrip <- -matchedDF_all$d18O_ngrip

### models - Carbon ----

# Start simple and build complexity

# climate-only model
carbon.lm.clim<-lm(del13C_permil~d18O_ngrip, matchedDF_all)
summary(carbon.lm.clim)
plot(del13C_permil ~ d18O_ngrip, data=matchedDF_all, pch=16)
abline(carbon.lm.clim)

# taxon-only model
carbon.lm.taxon<-lm(del13C_permil~Taxon, data=matchedDF_all)
summary(carbon.lm.taxon)

# model with no interaction term 
# NOTE: this is what the final model is in the end, so should be the same as carbon.lm.final
carbon.lm.all.additive<-lm(del13C_permil~d18O_ngrip + Taxon, data=matchedDF_all)
summary(carbon.lm.all.additive)

# model with interaction term included
carbon.lm.all.interaction<-lm(del13C_permil~d18O_ngrip * Taxon, 
                              data=matchedDF_all)
summary(carbon.lm.all.interaction)

##NOTE - no significant interactions, total amount of variation explained is about same as additive model

#stepwise regression 
carbon.lm.final <- step(lm(del13C_permil~d18O_ngrip * Taxon, data=matchedDF_all), direction="both")
summary(carbon.lm.final)
# --> this shows that the d180 + Taxon model (no interaction, additive model) is the best final model. 

# t-test of residuals from climate only model
# residuals of Otospermophilus and Sylvilagus are signif different
c.t.clim <- t.test(carbon.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_ngrip))])
c.t.clim 

### models - Nitrogen ----

# Start simple and build complexity

# climate-only model
nitrogen.lm.clim<-lm(del15N_permil ~ d18O_ngrip, data=matchedDF_all)
summary(nitrogen.lm.clim)

# taxon-only model
nitrogen.lm.taxon<-lm(del15N_permil~Taxon, data=matchedDF_all)
summary(nitrogen.lm.taxon)

# model with no interaction term 
nitrogen.lm.all.additive<-lm(del15N_permil ~ d18O_ngrip + Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.additive)

# model with interaction term included
nitrogen.lm.all.interaction<-lm(del15N_permil ~ d18O_ngrip * Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.interaction)

#stepwise regression --> this shows that there is not a great final model
nitrogen.lm.final <- step(lm(del15N_permil~d18O_ngrip*Taxon, data=matchedDF_all), direction="both")
summary(nitrogen.lm.final)
# --> the best final model is the additive model, which is the same as carbon.
# --> Not a very strong movel overall - adj R2 is low, no vars sign.

# t-test of residuals from climate-only model
# residuals of Otospermophilus and Sylvilagus are NOT signif different
n.t.clim <- t.test(nitrogen.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_ngrip))])
n.t.clim

### Do we need to add age as a co-variate in the models? ----
# Does this answer the question: "does isotope vary through time?"
carbon.lm.all.age<-step(lm(del13C_permil~d18O_ngrip * Taxon * median_age, data=matchedDF_all), direction="both")
summary(carbon.lm.all.age)
# explain a bit more variation with age included: 0.4554 vs 0.3714

nitrogen.lm.all.age<-step(lm(del15N_permil ~ d18O_ngrip * Taxon * median_age, data=matchedDF_all), direction="both")
summary(nitrogen.lm.all.age)
# explain a bit more variation with age included: 0.08062 vs 0.05773

final.model.stats[3,] <- c("carbon", "ngrip", "carbon.lm.all.additive", 
                           Regressionp(carbon.lm.all.additive), summary(carbon.lm.all.additive)$adj)
final.model.stats[4,] <- c("nitrogen", "ngrip", "nitrogen.lm.all.additive", 
                           Regressionp(nitrogen.lm.all.additive), summary(nitrogen.lm.all.additive)$adj)

# rename final models 
carbon.lm.final.ngrip <- carbon.lm.final
carbon.lm.clim.ngrip <- carbon.lm.clim
carbon.lm.all.age.ngrip <- carbon.lm.all.age
c.t.clim.ngrip <- c.t.clim
nitrogen.lm.final.ngrip <- nitrogen.lm.final
nitrogen.lm.clim.ngrip <- nitrogen.lm.clim
nitrogen.lm.all.age.ngrip <- nitrogen.lm.all.age
n.t.clim.ngrip <- n.t.clim

# Save model objects for use in plotting
save(carbon.lm.final.ngrip, carbon.lm.clim.ngrip, carbon.lm.all.age.ngrip, c.t.clim.ngrip, 
     nitrogen.lm.final.ngrip, nitrogen.lm.clim.ngrip, nitrogen.lm.all.age.ngrip, n.t.clim.ngrip, 
     file="output/ngrip_models_for_plotting.RData")

# remove variables 
rm("c.t.clim", "carbon.lm.all.additive", "carbon.lm.all.age", "carbon.lm.all.interaction", "carbon.lm.clim", "carbon.lm.final", "carbon.lm.taxon",
   "n.t.clim", "nitrogen.lm.all.additive", "nitrogen.lm.all.age", "nitrogen.lm.all.interaction", "nitrogen.lm.clim", "nitrogen.lm.final", "nitrogen.lm.taxon")


## Summarize results ----
print(final.model.stats[c(1,3,2,4),])
summary(carbon.lm.final.ngrip)
summary(carbon.lm.final.hendy)
# for carbon, the hendy data better explains variation in carbon isotope values. this is seen for the full model, as well as the greater partial variation explained by d18O. In the ngrip models, the strength of Taxon relative to d18O increases in explaining variation in carbon isotopes.
summary(nitrogen.lm.final.ngrip)
summary(nitrogen.lm.final.hendy)
# for nitrogen, ngrip and hendy data produce virtually identical and very weak models, both overall as well as in terms of the partial effect of d18O on carbon isotope variation.
temp.table <- final.model.stats[c(1,3,2,4), -3]
temp.table[,3] <- round(as.numeric(temp.table[,3]), 4)
temp.table[,4] <- round(as.numeric(temp.table[,4]), 4)
write.csv(temp.table, file="output/SupplementaryTable2.csv", row.names=F)

# Based on this, move forward with Hendy dataset




# Final additional plot - age of Squirrels vs rabbits ----
age_ttest<- t.test(specimen_medianage~Taxon, data=matchedDF_all)

pdf(file="output/SuppFigY_taxon_age_distribution.pdf", height=6, width=6)
boxplot(specimen_medianage~Taxon, data=matchedDF_all, 
        ylim=c(55000, 0), 
        xlab="", ylab="Median Specimen Age (cal years BP)")
stripchart(specimen_medianage~Taxon, data=matchedDF_all, vertical=TRUE, ylim=c(55000, 0), add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
legend("topright", legend=paste0("t=", round(age_ttest$statistic,2), "; df=", round(age_ttest$parameter,2), "; p=", round(age_ttest$p.value,2)), bty = "n", cex = 0.8)
dev.off()
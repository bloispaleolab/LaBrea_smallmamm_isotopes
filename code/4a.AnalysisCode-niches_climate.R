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
n.t.clim.hendy <- n.t.clim

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
n.t.clim.ngrip <- n.t.clim

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

## Final plots ####

### Figure 3 ----

#grDevices::pdf("output/Figure3_lm_carbon_nitrogen_all_July2021_NF.pdf", width=8, height=8)
grDevices::cairo_pdf("output/Figure3_lm_carbon_nitrogen_all.pdf", width=8, height=8)

layout(matrix(seq(1:6), ncol=2, nrow=3, byrow=F), heights=c(2.5,2.5,1))
par(mar=c(4,4,1,1), cex.axis=1, bty="l")


# carbon
plot(del13C_permil~d18O_hendy, data=matchedDF_all, pch=16, type="n", xlab="", ylab="")
points(del13C_permil~d18O_hendy, 
       data=matchedDF_all[which(matchedDF_all$Taxon == "Sylvilagus"),], 
       pch=16, col="darkorange")
points(del13C_permil~d18O_hendy, 
       data=matchedDF_all[which(matchedDF_all$Taxon == "Otospermophilus"),], 
       pch=16, col="royalblue2")
abline(carbon.lm.final.hendy, lty=2)
abline(carbon.lm.clim.hendy, lty=1)
mtext(expression({delta}^18*O~'value ('~'\u2030'~', PDB)'), side=1, line=2.5, cex=0.75) 
mtext(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)'), side=2, line=2.25, cex=0.75)

boxplot(carbon.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))], 
        xlab="", ylab="")
stripchart(carbon.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
mtext("Taxon", side=1, line=2.25)
mtext(expression('Residuals (Carbon Climate-only Model)'), side=2, line=2.25, cex=0.8)
# mtext(expression('Residuals ('~{delta}^13*C~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
legend("topright", legend=paste0("t=", round(c.t.clim.hendy$statistic,2), "; df=", round(c.t.clim.hendy$parameter,2), "; p=", round(c.t.clim.hendy$p.value,2)), bty = "n", cex = 0.8)

# taxon legend
# Draw an empty plot
plot(5, 5,
     type="n", axes=FALSE, ann=FALSE,
     xlim=c(0, 10), ylim = c(0,10))
legend("left", xpd=T, ncol=2, 
       legend = c("Otospermophilus", "Sylvilagus"),
       title="Taxon", title.adj=0, cex = 1.25,
       col = c("royalblue2","darkorange"), pch = 16,
       bty = "n")

# nitrogen
plot(del15N_permil~d18O_hendy, data=matchedDF_all, pch=16, 
     xlab = "", ylab = "", type="n")
points(del15N_permil~d18O_hendy, 
       data=matchedDF_all[which(matchedDF_all$Taxon == "Sylvilagus"),], 
       pch=16, col="darkorange")
points(del15N_permil~d18O_hendy, 
       data=matchedDF_all[which(matchedDF_all$Taxon == "Otospermophilus"),], 
       pch=16, col="royalblue2")
abline(nitrogen.lm.final.hendy, lty=2)
abline(nitrogen.lm.clim.hendy, lty=1)
mtext(expression({delta}^18*O~'value ('~'\u2030'~', PDB)'), side=1, line=2.5, cex=0.75)
mtext(expression({delta}^15*N~'value ('~'\u2030'~', AIR)'), side=2, line=2.25, cex=0.75)

boxplot(nitrogen.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))],
        xlab="", ylab="")
stripchart(nitrogen.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
mtext("Taxon", side=1, line=2.25)
mtext(expression('Residuals (Nitrogen Climate-only Model)'), side=2, line=2.25, cex=0.8)
# mtext(expression('Residuals ('~{delta}^15*N~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
legend("topright", legend=paste0("t=", round(n.t.clim.hendy$statistic,2), "; df=", round(n.t.clim.hendy$parameter,2), "; p=", round(n.t.clim.hendy$p.value,2)), bty = "n", cex = 0.8)


# model legend
plot(5, 5, 
     type="n", axes=FALSE, ann=FALSE, 
     xlim=c(0, 10), ylim = c(0,10))
legend("left", xpd=T,
       legend = c("Climate+Taxon", "Climate-only"),
       title="Model", title.adj=0, 
       lty = c(2,1), cex=1.25, ncol=2,
       bty = "n")

dev.off()

# I haven't put A-D labels on them yet.
# Figure caption text.
# Figure 3. The relationship between isotope niche and climate. A) & C) show the relationship between 13C or 15N, respectively, and 18O.  In both panels, the dashed line indicates the fitted relationship between 13C or 15N and 18O from the final model which  includes Taxon as an independent variable. The solid line indicates the fitted relationship between 13C or 15N and 18O from a linear model that just includes 18O as the independent variable. B) & D) indicate the residuals from the climate-only model, plotted by taxon.

### Figure 4 ----
#plot carbon and climate thru time

# order matchedDF_all
matchedDF_all <- matchedDF_all[order(matchedDF_all$median_age),]

#grDevices::pdf(file="output/Figure4_carbon_climate_time_updated_NF.pdf", height=6, width=8)
grDevices::cairo_pdf(file="output/isotope paper final/Figure4_carbon_climate_time_updated_Jan2022_JB.pdf", height=6, width=8)

layout(matrix(seq(1:2), ncol=1, nrow=2), heights=c(0.75,1))
par(mar=c(4,5,0,5), cex.axis=1, bty="l", xpd=F)

# d13C plot
plot(del13C_permil~median_age, data=matchedDF_all, 
     type="n", axes=FALSE, ann=FALSE, xaxs="i", yaxs="r", 
     xlim=c(55000, 0))
axis(4)
mtext(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)'), side=4, line=2.25)
lines(del13C_permil~median_age, data=matchedDF_all, lty=1, col="black")
points(del13C_permil~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Otospermophilus"),], pch=16, col="royalblue2", cex=0.5)
points(del13C_permil~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Sylvilagus"),], pch=16, col="darkorange", cex=0.5)


#d18O plot - Hendy
plot(pach.d18O~HendyAge, dat=hendyDat,
     type="l",
     xlab="Years before present",
     ylab = expression({delta}^18*O~'value ('~'\u2030'~', PDB)'),
     bty="n",
     xlim=c(55000, 0), ylim=c(3, 0),
     lab=c(12, 8, 7), xaxs="i", yaxs="r", col="lightgray")
x1.05<- loess(hendyDat$pach.d18O~hendyDat$HendyAge, span=0.05)
lines(x1.05$fitted~x1.05$x, type="l", col="darkgray")

points(d18O_hendy~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Sylvilagus"),], pch=16, col="darkorange", cex=0.5)
points(d18O_hendy~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Otospermophilus"),], pch=16, col="royalblue2", cex=0.5)
dev.off()




# Final additional plot - age of Squirrels vs rabbits ----
age_ttest<- t.test(specimen_medianage~Taxon, data=matchedDF_all)

pdf(file="output/SuppFigY_taxon_age_distribution.pdf", height=6, width=6)
boxplot(specimen_medianage~Taxon, data=matchedDF_all, 
        ylim=c(55000, 0), 
        xlab="", ylab="Median Specimen Age (cal years BP)")
stripchart(specimen_medianage~Taxon, data=matchedDF_all, vertical=TRUE, ylim=c(55000, 0), add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
legend("topright", legend=paste0("t=", round(age_ttest$statistic,2), "; df=", round(age_ttest$parameter,2), "; p=", round(age_ttest$p.value,2)), bty = "n", cex = 0.8)
dev.off()
## Final plots ####

# Figure 1 plot - in script 4b
# Fiure 2 plot

# Read in data ----
matchedDF_all <- read.csv(file="data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T) # 
load(file="output/hendy_models_for_plotting.RData") # saved models

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


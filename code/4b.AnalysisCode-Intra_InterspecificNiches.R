options(scipen=999)
library(ggplot2)
library(tidyverse)
library(ggpubr)


# library(scales)
# library(RColorBrewer)
# library(rstatix)

# Read in file
iso_dat <- read.csv(file="data/processed/final_dataset_focaltaxa_finalAMSruns_with_isotopes_with_all_dates.csv", header=T) # All rabbits and Squirrels with isotopes - some have no dates, or dates too old to calibrate 

## Inter-specific - Compare squirrels to rabbits ----

# examine data
#carbon
wilcox.c.taxon <- wilcox.test(del13C_permil~Taxon, data=iso_dat, correct = TRUE)
wilcox.c.taxon

p3 <- ggboxplot(iso_dat, x = "Taxon", y = "del13C_permil",
                color = "Taxon", palette = c("royalblue2", "darkorange"),
                add = "jitter") +
  labs(y = expression({delta}^13*C~'\u2030'), x="")+
  #  labs(y = expression(italic(delta)^13*C~("\211")))+
  theme(axis.title.x=element_blank(),
        legend.position = "none")
#  Add p-value
genus_C <- p3 + stat_compare_means(method = "wilcox.test", label.x = 2)
genus_C

#nitrogen
wilcox.n.taxon <- wilcox.test(del15N_permil~Taxon, data=iso_dat, correct = TRUE)
wilcox.n.taxon
p4 <- ggboxplot(iso_dat, x = "Taxon", y = "del15N_permil",
                color = "Taxon", palette = c("royalblue2", "darkorange"),
                add = "jitter") +
  #  labs(y = expression(italic(delta)^15*N~("\211")))+
  labs(y = expression(italic(delta)^15*N~("\u2030")))+
  theme(axis.title.x=element_blank(),
        legend.position = "none")
#  Add p-value
genus_N <- p4 + stat_compare_means(method = "wilcox.test", label.x = 2)
genus_N

## Intra - specific diet ---- 

## Compare desert cottontails to brush rabbits

iso_rabbits<- iso_dat[grep("S. audubonii|S. bachmani", iso_dat$Species), c('del13C_permil', 'del15N_permil', 'Species', 'Taxon')]

#carbon
wilcox.c.rabbits <- wilcox.test(del13C_permil~Species, data=iso_rabbits, correct = TRUE, alternative = "two.sided")
wilcox.c.rabbits
#t.test(del13C_permil~Taxon, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F, data=iso_rabbits)

p1 <- ggboxplot(iso_rabbits, x = "Species", y = "del13C_permil",
                color = "Species", palette = c("orange", "darkorange3"),
                add = "jitter") +
  labs(y = expression({delta}^13*C~'\u2030'))+
  theme(axis.title.x=element_blank(),
        legend.position = "none")
#  Add p-value
rabbits_C <- p1 + stat_compare_means(method = "wilcox.test", label.x = 2)
#p-values appear on figure as a single plot, but not on multiplots
rabbits_C

#nitrogen
wilcox.n.rabbits <- wilcox.test(del15N_permil~Species, data=iso_rabbits, correct = TRUE)
wilcox.n.rabbits

p2 <- ggboxplot(iso_rabbits, x = "Species", y = "del15N_permil",
                color = "Species", palette = c("orange", "darkorange3"),
                add = "jitter") +
  #  labs(y = expression(italic(delta)^15*N~("\211")))+
  labs(y = expression(italic(delta)^15*N~("\u2030")))+
  theme(axis.title.x=element_blank(),
        legend.position = "none")
#  Add p-value
rabbits_N <- p2 + stat_compare_means(method = "wilcox.test", label.x = 2)
rabbits_N

figure <- ggarrange(genus_C, rabbits_C, genus_N, rabbits_N,
                    labels = c("A", "C", "B", "D"),
                    ncol = 2, nrow = 2)
figure

# summarize basic stats for the paper ---- 
# TAXON
wilcox.c.taxon
wilcox.n.taxon

iso_dat %>% count(Taxon)

# carbon results by Taxon
C.taxon.res <- iso_dat %>%      # Summary by group using dplyr
  group_by(Taxon) %>% 
  summarize(min = min(del13C_permil),
            median = median(del13C_permil),
            max = max(del13C_permil),
            mean = mean(del13C_permil))
C.taxon.res
C.taxon.res$mean[1] - C.taxon.res$mean[2]

  
# nitrogen results by Taxon
N.taxon.res <- iso_dat %>%                               # Summary by group using dplyr
  group_by(Taxon) %>% 
  summarize(min = min(del15N_permil),
            median = median(del15N_permil),
            mean = mean(del15N_permil),
            max = max(del15N_permil))
N.taxon.res
N.taxon.res$mean[1] - N.taxon.res$mean[2]


# SPECIES
wilcox.c.rabbits
wilcox.n.rabbits

iso_rabbits %>% count(Species)

# carbon results by Taxon
C.rabbits.res <- iso_rabbits %>%                               # Summary by group using dplyr
  group_by(Species) %>% 
  summarize(min = min(del13C_permil),
            median = median(del13C_permil),
            mean = mean(del13C_permil),
            max = max(del13C_permil))
C.rabbits.res
C.rabbits.res$mean[1] - C.rabbits.res$mean[2]

# nitrogen results by Taxon
N.rabbits.res <- iso_rabbits %>%                               # Summary by group using dplyr
  group_by(Species) %>% 
  summarize(min = min(del15N_permil),
            median = median(del15N_permil),
            mean = mean(del15N_permil),
            max = max(del15N_permil))
N.rabbits.res
N.rabbits.res$mean[1] - N.rabbits.res$mean[2]


# Figure 1 Final ----
 # to sync up Sylvilagus color scheme across panels

# create new dataframes and colors
#data
Cdatalist <- list(iso_dat$del13C_permil[which(iso_dat$Taxon=="Otospermophilus")],
                  iso_dat$del13C_permil[which(iso_dat$Taxon=="Sylvilagus")])

Ndatalist <- list(iso_dat$del15N_permil[which(iso_dat$Taxon=="Otospermophilus")],
                  iso_dat$del15N_permil[which(iso_dat$Taxon=="Sylvilagus")])

# colors
oto.col<- rep("royalblue2", length(Ndatalist[[1]]))
syl.col <- iso_dat$Species[which(iso_dat$Taxon=="Sylvilagus")]
syl.col[which(syl.col=="Leporidae")] <- "darkorange" # assuming this is Sylvilagus sp.
syl.col[which(syl.col=="Sylvilagus sp")] <- "darkorange"
syl.col[which(syl.col=="S. bachmani")] <- "red1"
syl.col[which(syl.col=="S. audubonii")] <- "gold1"
colorlist <- list(oto.col, syl.col)

grDevices::cairo_pdf("output/Figure1_niche_boxplots_Feb2022_JB.pdf", width=8, height=8)
par(bty="l", mfcol=c(2,2), mar=c(4,4,1,1)+0.1)

## Inter-specific ----
# Compare squirrels to rabbits

# carbon, interspecific
# create basic boxplot
boxplot(del13C_permil ~ Taxon, data=iso_dat, col = "white", ylab = "", cex.axis=0.75)

# add individual data points
for (i in 1:2) { 
  stripchart(na.omit(Cdatalist[[i]]), at = i, add = T, 
             pch =21, bg = colorlist[[i]], lwd=0.5, vertical = T, method = "jitter") 
}
mtext(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)'), side=2, line = 2.25)
legend("top", legend=paste0("W = ", wilcox.c.taxon$statistic, "; p = ", round(wilcox.c.taxon$p.value, 3)), bty="n", cex=0.75)

legend("topright", legend=c('S. bachmani', 'Sylvilagus sp.', 'S. audubonii'), pch=21, col="black", pt.lwd=0.5, pt.bg=c("red1", "darkorange", "gold1"), bty="n", cex=0.5)

# nitrogen, interspecific
# create basic boxplot
boxplot(del15N_permil ~ Taxon, data=iso_dat, col = "white", ylab = "", cex.axis=0.75)

# add individual data points
for (i in 1:2) { 
  stripchart(na.omit(Ndatalist[[i]]), at = i, add = T, pch = 21, 
             bg = colorlist[[i]], lwd=0.5, vertical = T, method = "jitter") 
}
mtext(expression({delta}^15*N~'value ('~'\u2030'~', AIR)'), side=2, line = 2.25)
legend("top", legend=paste0("W = ", wilcox.n.taxon$statistic, "; p = ", round(wilcox.n.taxon$p.value, 3)), bty="n", cex=0.75)

## Intra-specific diet ---- 

## Compare desert cottontails to brush rabbits

# carbon, introspecific
# create basic boxplot
boxplot(del13C_permil ~ Species, data=iso_rabbits, col = "white", ylab = "", cex.axis=0.75)
# add individual data points
stripchart(del13C_permil ~ Species, data=iso_rabbits, vertical=TRUE, add=TRUE, method="jitter", pch=21, col="black", bg=c("gold1","red1"), lwd=0.5)
mtext(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)'), side=2, line = 2.25)
legend("top", legend=paste0("W = ", wilcox.c.rabbits$statistic, "; p = ", round(wilcox.c.rabbits$p.value, 3)), bty="n", cex=0.75)

# nitrogen, intraspecific
# create basic boxplot
boxplot(del15N_permil ~ Species, data=iso_rabbits, col = "white", ylab = "", cex.axis=0.75)
# add individual data points
stripchart(del15N_permil ~ Species, data=iso_rabbits, vertical=TRUE, add=TRUE, method="stack", pch=21, col="black", bg=c("gold1","red1"), lwd=0.5)
mtext(expression({delta}^15*N~'value ('~'\u2030'~', AIR)'), side=2, line = 2.25)
legend("top", legend=paste0("W = ", wilcox.n.rabbits$statistic, "; p = ", round(wilcox.n.rabbits$p.value, 3)), bty="n", cex=0.75)

dev.off()


# In Illustrator: Add A,B,C,D labels, italicize taxon names, move jittered points over to outliers.

# MANOVA models ----
data1 <- read.csv(file="data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T) # 
# Summary stats ----

cn.manova <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon, data1)
summary(cn.manova)

# factor two of the grouping variables the run t-tests
data2 <- data1
data2$Taxon <- as.factor(data2$Taxon)
data2$time_group <- as.factor(data2$time_group)

# t.test
c.t <- t.test(del13C_permil~Taxon, data=data2)
c.t
n.t <- t.test(del15N_permil~Taxon, data=data2)
n.t

# Just pleistocene data
data3 <- data2[which(data2$time_group=="Pleistocene"),]
cn.manova.p <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon, data3)
summary(cn.manova.p)

# t.test
c.t.p <- t.test(del13C_permil~Taxon, data=data3)
c.t.p
n.t.p <- t.test(del15N_permil~Taxon, data=data3)
n.t.p

# try adding sample age to the manova
cnt.manova <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon + median_age, data1)
summary(cnt.manova)

cnt.manova.p <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon + median_age, data3)
summary(cnt.manova.p)

# try adding sample age to the manova
cnt.manova.interact <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon * median_age, data1)
summary(cnt.manova.interact)

cnt.manova.p.interact <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon * median_age, data3)
summary(cnt.manova.p.interact)

## Table 2 ----
manova.res.PH <- summary(cnt.manova)$stats
manova.res.P <- summary(cnt.manova.p)$stats
save(manova.res.PH, manova.res.P, file = "output/manovaResults.RData")
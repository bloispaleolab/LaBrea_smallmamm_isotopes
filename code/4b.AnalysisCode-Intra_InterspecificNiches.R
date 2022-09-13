options(scipen=999)
library(ggplot2)
library(tidyverse)
library(ggpubr)


# library(scales)
# library(RColorBrewer)
# library(rstatix)

# Read in file
iso_dat <- read.csv(file="data/processed/final_dataset_focaltaxa_finalAMSruns_with_isotopes_with_all_dates.csv", header=T) # All rabbits and Squirrels with isotopes - some have no dates, or dates too old to calibrate 

# determine normality, Shapiro-Wilk tests ----

#carbon
CS <- iso_dat %>%
  filter(Taxon == 'Sylvilagus') %>%
  select(del13C_permil)
shapiro.test(CS$del13C_permil)
# not normal

CO <- iso_dat %>%
  filter(Taxon == 'Otospermophilus') %>%
  select(del13C_permil)
shapiro.test(CO$del13C_permil)
# not normal

NS <- iso_dat %>%
  filter(Taxon == 'Sylvilagus') %>%
  select(del15N_permil)
shapiro.test(NS$del15N_permil)
# normal

NO <- iso_dat %>%
  filter(Taxon == 'Otospermophilus') %>%
  select(del15N_permil)
shapiro.test(NO$del15N_permil)
# normal

CSa <- iso_dat %>%
  filter(Species == 'S. audubonii') %>%
  select(del13C_permil)
shapiro.test(CSa$del13C_permil)
# normal

CSb <- iso_dat %>%
  filter(Species == 'S. bachmani') %>%
  select(del13C_permil)
shapiro.test(CSb$del13C_permil)
# not normal

NSa <- iso_dat %>%
  filter(Species == 'S. audubonii') %>%
  select(del15N_permil)
shapiro.test(NSa$del15N_permil)
# not normal

NSb <- iso_dat %>%
  filter(Species == 'S. bachmani') %>%
  select(del15N_permil)
shapiro.test(NSb$del15N_permil)
# normal


## Inter-specific - Compare squirrels to rabbits ----

# examine data
#carbon - both taxa non-normal, use wilcon.test
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

#nitrogen - could use normal t-test for this. t.test marginally not significant, just use wilcox test for consistency with other comparisons.
t.test(del15N_permil~Taxon, data=iso_dat)
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

#carbon - one taxon normal, other non-normal. Use wilcox.test 
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

#nitrogen - one taxon normal, other non-normal. Use wilcox.test
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
#data - genus
Cdatalist <- list(iso_dat$del13C_permil[which(iso_dat$Taxon=="Otospermophilus")],
                  iso_dat$del13C_permil[which(iso_dat$Taxon=="Sylvilagus")])

Ndatalist <- list(iso_dat$del15N_permil[which(iso_dat$Taxon=="Otospermophilus")],
                  iso_dat$del15N_permil[which(iso_dat$Taxon=="Sylvilagus")])

# colors - genus
oto.col<- rep("royalblue2", length(Ndatalist[[1]]))
syl.col <- iso_dat$Species[which(iso_dat$Taxon=="Sylvilagus")]
syl.col[which(syl.col=="Leporidae")] <- "darkorange" # assuming this is Sylvilagus sp.
syl.col[which(syl.col=="Sylvilagus sp")] <- "darkorange"
syl.col[which(syl.col=="S. bachmani")] <- "red1"
syl.col[which(syl.col=="S. audubonii")] <- "gold1"
colorlist <- list(oto.col, syl.col)

#data - rabbits
S_Cdatalist <- list(iso_rabbits$del13C_permil[which(iso_rabbits$Species=="S. audubonii")],
                    iso_rabbits$del13C_permil[which(iso_rabbits$Species=="S. bachmani")])

S_Ndatalist <- list(iso_rabbits$del15N_permil[which(iso_rabbits$Species=="S. audubonii")],
                    iso_rabbits$del15N_permil[which(iso_rabbits$Species=="S. bachmani")])

# colors - rabbits
sa.col <- rep("gold1", length(which(iso_rabbits$Species == "S. audubonii"))) 
sb.col <- rep("red1", length(which(iso_rabbits$Species == "S. bachmani")))
rabbit_colorlist <- list(sa.col, sb.col)

#pdf(file="output/Figure1_niche_boxplots_Mar2022_NF.pdf", height=8, width=8)
grDevices::cairo_pdf("output/Figure1_niche_boxplots_Sep2022_JB.pdf", width=8, height=8)
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
for (i in 1:2) { 
  stripchart(na.omit(S_Cdatalist[[i]]), at = i, add = T, pch = 21, 
             bg = rabbit_colorlist[[i]], lwd=0.5, vertical = T, method = "jitter") 
}
mtext(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)'), side=2, line = 2.25)
legend("top", legend=paste0("W = ", wilcox.c.rabbits$statistic, "; p = ", round(wilcox.c.rabbits$p.value, 3)), bty="n", cex=0.75)

# nitrogen, intraspecific
# create basic boxplot
boxplot(del15N_permil ~ Species, data=iso_rabbits, col = "white", ylab = "", cex.axis=0.75)
# add individual data points
for (i in 1:2) { 
  stripchart(na.omit(S_Ndatalist[[i]]), at = i, add = T, pch = 21, 
             bg = rabbit_colorlist[[i]], lwd=0.5, vertical = T, method = "jitter") 
}
mtext(expression({delta}^15*N~'value ('~'\u2030'~', AIR)'), side=2, line = 2.25)
legend("top", legend=paste0("W = ", wilcox.n.rabbits$statistic, "; p = ", round(wilcox.n.rabbits$p.value, 3)), bty="n", cex=0.75)

dev.off()


# stripchart(del13C_permil ~ Species, data=iso_rabbits, vertical=TRUE, add=TRUE, method="jitter", pch=21, col="black", bg=c("gold1","red1"), lwd=0.5)
# stripchart(del15N_permil ~ Species, data=iso_rabbits, vertical=TRUE, add=TRUE, method="stack", pch=21, col="black", bg=c("gold1","red1"), lwd=0.5)

# In Illustrator: Add A,B,C,D labels, italicize taxon names, move jittered points over to outliers.

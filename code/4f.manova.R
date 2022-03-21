# MANOVA stats - originally in SIBER script----

# load in the isotope dataset
dat<- read.csv("data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T)

# factor two of the grouping variables the run t-tests
dat$Taxon <- as.factor(dat$Taxon)
dat$time_group <- as.factor(dat$time_group)

# does Taxon significantly explain joint variation in C and N? Yes
cn.manova.all <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon, dat)
summary(cn.manova.all)

# t.test - does taxon explain C and N separately?
c.t <- t.test(del13C_permil~Taxon, data=dat)
c.t # taxa are sig different in Carbon
n.t <- t.test(del15N_permil~Taxon, data=dat)
n.t # taxa are NOT sig different in Nitrogen

# Just pleistocene data
# does Taxon significantly explain joint variation in C and N during the Pleistocene? Yes
dat_p <- dat[which(dat$time_group=="Pleistocene"),]
cn.manova.p <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon, dat_p)
summary(cn.manova.p)

# t.test
c.t.p <- t.test(del13C_permil~Taxon, data=dat_p)
c.t.p # taxa are sig different in Carbon in the Pleisto
n.t.p <- t.test(del15N_permil~Taxon, data=dat_p)
n.t.p # taxa are NOT sig different in Nitrogen in the Pleisto

# try adding sample age to the manova
# both taxon and age are significant in explaining joint variation in C and N across Quaternary
cnt.manova.all <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon + median_age, dat)
summary(cnt.manova.all)

# only taxon is significant in explaining joint variation in C and N across Quaternary. 
# age is close (0.0515)
cnt.manova.p <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon + median_age, dat_p)
summary(cnt.manova.p)

# try adding sample age to the manova as an interaction
# interaction term doesn't change anything. age still only sig with full Quaternary data, only marginal non-sig during Pleisto. 
cnt.manova.all.interact <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon * median_age, dat)
summary(cnt.manova.all.interact)

cnt.manova.p.interact <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon * median_age, dat_p)
summary(cnt.manova.p.interact)



# MANOVA models - originally in intra-interspecific script----
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

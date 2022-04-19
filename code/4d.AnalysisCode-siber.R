library(tidyverse)
library(SIBER)

# SIBER Stats ----
source('code/maxLikOverlap.PH.R')

# load in the isotope dataset
dat<- read.csv("data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T)

## Create data frames for SIBER analysis ----
# First create a df with groups (Sylv and Oto) and communities (Holo vs Pleisto) separated. 
siber.data.all <- dat %>% 
  rename(iso1 = del13C_permil,
        iso2 = del15N_permil,
        community = time_group,
        group = Taxon) %>%
    select(iso1, iso2, group, community)

siber.data.all <- siber.data.all %>%
  mutate(community = factor(community),
         group = factor(group))

# create object that breaks out data into 4 separate groups
siber.data.all.obj <- createSiberObject(siber.data.all)
# Note: The warning is because sample sizes are too low for Holo for the individual taxa

# create new object that groups Pleisto and Holo into one Quaternary community, taxa separate
siber.data.q <- siber.data.all
siber.data.q$community <- "Quaternary"
siber.data.q.obj <- createSiberObject(siber.data.q)

# Then create separate Pleistocene dataframe, with taxa still separate
siber.data.p <- siber.data.all[which(siber.data.all$community == "Pleistocene"),]
siber.data.p$community <- droplevels(siber.data.p$community)
siber.data.p.obj <- createSiberObject(siber.data.p)

# Compare all taxa for Pleistocene and Holocene separately
# For this, we need to switch groups and communities
siber.data.PH <- dat %>% 
  rename(iso1 = del13C_permil,
         iso2 = del15N_permil,
         community = Taxon,
         group = time_group) %>%
  select(iso1, iso2, group, community)

siber.data.PH$community <- "AllTaxa"

siber.data.PH <- siber.data.PH %>%
  mutate(community = factor(community),
         group = factor(group))

siber.data.PH.obj <- createSiberObject(siber.data.PH)
# Note: The warning is because sample sizes are too low for Holo for the individual taxa

# create one huge group - AllTaxa and Quaternary
siber.data.Q.AllTaxa <- siber.data.PH
siber.data.Q.AllTaxa$group <- "Quaternary"
siber.data.Q.AllTaxa.obj <- createSiberObject(siber.data.Q.AllTaxa)

# calculate SIBER metrics ----
group.ML.siber.all.obj <- groupMetricsML(siber.data.all.obj)
print(group.ML.siber.all.obj)

group.ML.siber.data.q.obj<- groupMetricsML(siber.data.q.obj)
print(group.ML.siber.data.q.obj)

group.ML.siber.data.p.obj <- groupMetricsML(siber.data.p.obj)
print(group.ML.siber.data.p.obj)

group.ML.siber.data.PH.obj<- groupMetricsML(siber.data.PH.obj)
print(group.ML.siber.data.PH.obj)

group.ML.siber.data.Q.AllTaxa.obj<- groupMetricsML(siber.data.Q.AllTaxa.obj)
print(group.ML.siber.data.Q.AllTaxa.obj)

# create summary table of group statistics
save(group.ML.siber.all.obj,  # data with all four groups
     group.ML.siber.data.p.obj, # data with taxa separate, just for Pleisto
     group.ML.siber.data.q.obj, # data with taxa separate, for Pleisto + Holo
     group.ML.siber.data.PH.obj, # data with taxa merged, for Pleisto and Holo separate
     file="output/siber.summary.stats.RData")

# start building Summary Stats Table - not used, just for comparison with final siber_output table
# Quaternary Sylvilagus
# Quaternary Otospermophilus
# Pleistocene Sylvilagus
# Pleistocene Otospermphilus
# Pleistocene small mammals
# Holocene small mammals
# Quaternary small mammals
summary.stats <- as.data.frame(matrix(data=NA, nrow=7, ncol=4)) 

summary.stats[1:2,2:4] <- t(group.ML.siber.data.q.obj)
summary.stats[1:2,1] <- colnames(group.ML.siber.data.q.obj)

summary.stats[3:4,2:4] <- t(group.ML.siber.data.p.obj)
summary.stats[3:4,1] <- colnames(group.ML.siber.data.p.obj)

summary.stats[5:6,2:4] <- t(group.ML.siber.data.PH.obj)
summary.stats[5:6,1] <- colnames(group.ML.siber.data.PH.obj)

summary.stats[7,2:4] <- t(group.ML.siber.data.Q.AllTaxa.obj)
summary.stats[7,1] <- colnames(group.ML.siber.data.Q.AllTaxa.obj)

colnames(summary.stats) <- c('Group', rownames(group.ML.siber.data.Q.AllTaxa.obj))

# plot data with all four separate groups and communities ----
# Note - don't need to plot with SIBER for the paper, just doing this to visualize as I go
p.ell <- 0.68 #set ellipse area for plotting

community.hulls.args <- list(col = "black", lty = 1, lwd = 1) # time group
group.ellipses.args  <- list(n = 100, p.interval = p.ell, lty = 1, lwd = 2) # Sylv and Oto
group.hulls.args     <- list(col="gray", lty = 2) #taxon

rlbPalette <- palette(c("royalblue2","darkorange"))
par(mfrow=c(1,1))
plotSiberObject(siber.data.all.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
print(group.ML.siber.all.obj)

# plot data with two taxon groups and 1 Quaternary community
par(mfrow=c(1,1))
plotSiberObject(siber.data.q.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = list(n = 100, p.interval = p.ell, lty = 1, lwd = 2),
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
print(group.ML.siber.data.q.obj)

# plot data with two groups and 1 Pleisto community
par(mfrow=c(1,1))
plotSiberObject(siber.data.p.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = list(n = 100, p.interval = p.ell, lty = 1, lwd = 2),
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
print(group.ML.siber.data.p.obj)

# plot data with two time groups and 1 taxonomic community
par(mfrow=c(1,1))
plotSiberObject(siber.data.PH.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args, # Pleisto and Holo
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
print(group.ML.siber.data.PH.obj)


# Calculate SEA statistics ----

# data check: compare SIBER output with the standard ellipse stats to make sure plotting is ok
  # temp<- rlb_data[which(rlb_data$time_group == "Holocene"), c('d13C', 'd15N')]
  # Sigma <- cov(temp)
  # evalues <- eigen(Sigma, symmetric = T, only.values = TRUE)$values
  # p <- p.ell # make sure you're calculate the same area as plotting
  # r <- 2 * qf(p, 2, nrow(temp)-1)
  # a <- sqrt(r * evalues[1])
  # b <- sqrt(r * evalues[2]) # SIBER does not include r values in their calculation
  # ea_H <- pi*a*b
  # # compare with:
  # sigmaSEA(Sigma)
  # # Note, same, except that SIBER does not include r values in their calculation of a, b, and SEA. So I think we're fine plotting the ellipses, but using the groupML function to quantify SEA and SEAc.

# create data frames to store siber_output
siber_output <- as.data.frame(matrix(NA, nrow=0, ncol=8))
colnames(siber_output) <- c("Group", "n", "SEA", "SEAc", "Lower 95 CI", "Upper 95 CI", "Probability", "Proportion Overlap")

# Description of the different comparisons ----
# Note, all this code below is for SEA. This calculates niche size, not niche overlap.
# Comparisons we want to make:
#   4. Are Oto sig diff from Sylv in Quaternary? group.ML.siber.data.q.obj
#   3. Are Oto sig diff from Sylv in Pleistocene? group.ML.siber.data.p.obj
#   2. Are Holocene all taxa significantly different from Pleistocene all taxa ? siber.data.PH.obj
#   0. Are Pleistocene all taxa different from Quaternary all taxa? siber.data.Q.AllTaxa.obj
#   1. Is each subgroup significantly different (by taxon AND time)? siber.data.all.obj 

# For all of these, we will first fit Bayesian multivariate normal distributions to each group in the dataset. Then, we will use siberEllipses to calculate the bayesian Standard Ellipse Area for all groups, then calculate credible intervals and compare groups
# we can compare groups by plotting the groups and the distribution of their ellipses, as well as calculating the probability that one groups posterior distribution is smaller (or larger) than another
# generally following directions here: https://cran.r-project.org/web/packages/SIBER/vignettes/siber-comparing-populations.html

# Set bayesian parameters for all calculations

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-2

# Comparison 4. Are Oto sig diff from Sylv in Quaternary ----
# use siber.data.q.obj and group.ML.siber.data.q.obj

## Note, the code chunk below is written generally. Choose the comparison to make in the next lines 
siber.obj.data <- siber.data.q.obj
group.ML.data <- group.ML.siber.data.q.obj

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.obj.data, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.data), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.data[3,], col="red", pch = "x", lwd = 2)

# Calculate credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)
# convert for plotting
CI_95_lower <-c(SEA.B.credibles$V1[2,1], SEA.B.credibles$V2[2,1]) 
CI_95_upper <- c(SEA.B.credibles$V1[2,2], SEA.B.credibles$V2[2,2])
x<- c(1,2)
segments(x0=x, y0=CI_95_lower, x1=x, y1=CI_95_upper, col="red")

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

######### store output for later
ellipses.posterior.q <- ellipses.posterior
SEA.B.q <- SEA.B

## Comparison 4 niche size statistics ----
# Pr = probability HP=holo and pleis, Q= Quaternary, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_QS_lt_QO <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_QS_lt_QO)

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes)),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper, c(Pr_QS_lt_QO, "NaN")))
siber_output[1:2,1:7] <- siber_output_temp

## Comparison 4 niche overlap statistics ----
#Quaternary
obs.overlap.QS.QO <- maxLikOverlap("Quaternary.Sylvilagus", 
                                   "Quaternary.Otospermophilus", 
                                   siber.obj.data, p = 0.95, n = 100)
# calculate max likelihood proportion overlaps
prop.overlap.both <- as.numeric(obs.overlap.QS.QO["overlap"] / 
                                  (obs.overlap.QS.QO["area.1"] + obs.overlap.QS.QO["area.2"]))
prop.overlap.QS <- as.numeric(obs.overlap.QS.QO["overlap"] / 
                                obs.overlap.QS.QO["area.1"])
prop.overlap.QO <- as.numeric(obs.overlap.QS.QO["overlap"] / 
                                obs.overlap.QS.QO["area.2"])
obs.overlap.QS.QO <- c(obs.overlap.QS.QO, prop.overlap.both, prop.overlap.QS, prop.overlap.QO)
names(obs.overlap.QS.QO)[4:6] <- c('prop.overlap.both', 'prop.overlap.QS', 'prop.overlap.QO')

siber_output[1:2,8] <- c(prop.overlap.QS, prop.overlap.QO)

# calculate bayes overlap distribution
bayes.overlap.QS.QO <- bayesianOverlap("Quaternary.Sylvilagus", 
                                       "Quaternary.Otospermophilus",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)

bayes.prop.overlap.both <- as.numeric(bayes.overlap.QS.QO[,'overlap'] / 
                                        (bayes.overlap.QS.QO[,'area1'] + bayes.overlap.QS.QO[,'area2']))
bayes.prop.overlap.QS <- as.numeric(bayes.overlap.QS.QO[,'overlap'] / bayes.overlap.QS.QO[,'area1'])
bayes.prop.overlap.QO <- as.numeric(bayes.overlap.QS.QO[,'overlap'] / bayes.overlap.QS.QO[,'area2'])
bayes.overlap.QS.QO <- cbind(bayes.overlap.QS.QO, bayes.prop.overlap.both, 
                             bayes.prop.overlap.QS, bayes.prop.overlap.QO)
# save objects
save(obs.overlap.QS.QO, bayes.overlap.QS.QO, file="output/overlap_QS.QO.RData")


# Comparison 3: Are Oto sig diff from Sylv in Pleistocene ---- 
# use siber.data.p.obj and group.ML.siber.data.p.obj

## Note, the code chunk below is written generally. Choose the comparison to make in the next lines 
siber.obj.data <- siber.data.p.obj
group.ML.data <- group.ML.siber.data.p.obj

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.obj.data, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.data), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.data[3,], col="red", pch = "x", lwd = 2)

# Calculate credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)
CI_95_lower <-c(SEA.B.credibles$V1[2,1], SEA.B.credibles$V2[2,1]) 
CI_95_upper <- c(SEA.B.credibles$V1[2,2], SEA.B.credibles$V2[2,2])
x<- c(1,2)
segments(x0=x, y0=CI_95_lower, x1=x, y1=CI_95_upper, col="red")

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

######### store output for later
ellipses.posterior.p <- ellipses.posterior
SEA.B.p <- SEA.B

## Comparison 3 niche size statistics ----
# Pr = probability HP=holo and pleis, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_PS_lt_PO <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_PS_lt_PO)

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes)),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper, c(Pr_PS_lt_PO, "NaN")))
siber_output[3:4,1:7] <- siber_output_temp

## Comparison 3 niche overlap statistics ----
#Pleistocene
obs.overlap.PS.PO <- maxLikOverlap("Pleistocene.Sylvilagus", 
                                   "Pleistocene.Otospermophilus", 
                                   siber.obj.data, p = 0.95, n = 100)
# calculate max likelihood proportion overlaps
prop.overlap.both <- as.numeric(obs.overlap.PS.PO["overlap"] / 
                                  (obs.overlap.PS.PO["area.1"] + obs.overlap.PS.PO["area.2"]))
prop.overlap.PS <- as.numeric(obs.overlap.PS.PO["overlap"] / 
                                obs.overlap.PS.PO["area.1"])
prop.overlap.PO <- as.numeric(obs.overlap.PS.PO["overlap"] / 
                                obs.overlap.PS.PO["area.2"])
obs.overlap.PS.PO <- c(obs.overlap.PS.PO, prop.overlap.both, prop.overlap.PS, prop.overlap.PO)
names(obs.overlap.PS.PO)[4:6] <- c('prop.overlap.both', 'prop.overlap.PS', 'prop.overlap.PO')

siber_output[3:4,8] <- c(prop.overlap.PS, prop.overlap.PO)

# calculate bayes overlap distribution
bayes.overlap.PS.PO <- bayesianOverlap("Pleistocene.Sylvilagus", 
                                       "Pleistocene.Otospermophilus",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)

bayes.prop.overlap.both <- as.numeric(bayes.overlap.PS.PO[,'overlap'] / 
                                        (bayes.overlap.PS.PO[,'area1'] + bayes.overlap.PS.PO[,'area2']))
bayes.prop.overlap.PS <- as.numeric(bayes.overlap.PS.PO[,'overlap'] / bayes.overlap.PS.PO[,'area1'])
bayes.prop.overlap.PO <- as.numeric(bayes.overlap.PS.PO[,'overlap'] / bayes.overlap.PS.PO[,'area2'])
bayes.overlap.PS.PO <- cbind(bayes.overlap.PS.PO, bayes.prop.overlap.both, 
                             bayes.prop.overlap.PS, bayes.prop.overlap.PO)
# save objects
save(obs.overlap.PS.PO, bayes.overlap.PS.PO, file="output/overlap_PS.PO_comparison3.RData")


# Comparison 2: Holo taxa vs Pleisto taxa ----
# use siber.data.PH.obj and group.ML.siber.data.PH.obj

## Note, the code chunk below is written generally. Choose the comparison to make in the next lines 
siber.obj.data <- siber.data.PH.obj
group.ML.data <- group.ML.siber.data.PH.obj

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.obj.data, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.data), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.data[3,], col="red", pch = "x", lwd = 2)

# Calculate credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)
CI_95_lower <-c(SEA.B.credibles$V1[2,1], SEA.B.credibles$V2[2,1]) 
CI_95_upper <- c(SEA.B.credibles$V1[2,2], SEA.B.credibles$V2[2,2])
x<- c(1,2)
segments(x0=x, y0=CI_95_lower, x1=x, y1=CI_95_upper, col="red")

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

######### store output for later
ellipses.posterior.PH <- ellipses.posterior
SEA.B.PH <- SEA.B

## Comparison 2 niche size statistics ----
# Pr = probability HP=holo and pleis, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_P_lt_H <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_P_lt_H)

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes)),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper, c(Pr_P_lt_H, "NaN")))
siber_output[5:6,1:7] <- siber_output_temp

## Comparison 2 niche overlap statistics ----
#Pleistocene
obs.overlap.P.H <- maxLikOverlap.PH("Pleistocene", 
                                   "Holocene", 
                                 siber.obj.data, p = 0.95, n = 100)
# calculate max likelihood proportion overlaps
prop.overlap.both <- as.numeric(obs.overlap.P.H["overlap"] / 
                                  (obs.overlap.P.H["area.1"] + obs.overlap.P.H["area.2"]))
prop.overlap.P <- as.numeric(obs.overlap.P.H["overlap"] / 
                                obs.overlap.P.H["area.1"])
prop.overlap.H <- as.numeric(obs.overlap.P.H["overlap"] / 
                                obs.overlap.P.H["area.2"])
obs.overlap.P.H <- c(obs.overlap.P.H, prop.overlap.both, prop.overlap.P, prop.overlap.H)
names(obs.overlap.P.H)[4:6] <- c('prop.overlap.both', 'prop.overlap.P', 'prop.overlap.H')

siber_output[5:6,8] <- c(prop.overlap.P, prop.overlap.H)

# calculate bayes overlap distribution
bayes.overlap.P.H <- bayesianOverlap("AllTaxa.Pleistocene", 
                                       "AllTaxa.Holocene",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)

bayes.prop.overlap.both <- as.numeric(bayes.overlap.P.H[,'overlap'] / 
                                        (bayes.overlap.P.H[,'area1'] + bayes.overlap.P.H[,'area2']))
bayes.prop.overlap.P <- as.numeric(bayes.overlap.P.H[,'overlap'] / bayes.overlap.P.H[,'area1'])
bayes.prop.overlap.H <- as.numeric(bayes.overlap.P.H[,'overlap'] / bayes.overlap.P.H[,'area2'])
bayes.overlap.P.H <- cbind(bayes.overlap.P.H, bayes.prop.overlap.both, 
                             bayes.prop.overlap.P, bayes.prop.overlap.H)
# save objects
save(obs.overlap.P.H, bayes.overlap.P.H, file="output/overlap_P.H.RData")


# Comparison 0. All taxa, Quaternary ----
# use siber.data.Q.AllTaxa.obj and group.ML.siber.data.Q.AllTaxa.obj

## Note, the code chunk below is written generally. Choose the comparison to make in the next lines 
siber.obj.data <- siber.data.Q.AllTaxa.obj
group.ML.data <- group.ML.siber.data.Q.AllTaxa.obj

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.obj.data, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.data), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.data[3,], col="red", pch = "x", lwd = 2)

# Calculate credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)
CI_95_lower <-c(SEA.B.credibles$V1[2,1], SEA.B.credibles$V2[2,1]) 
CI_95_upper <- c(SEA.B.credibles$V1[2,2], SEA.B.credibles$V2[2,2])
x<- c(1,2)
segments(x0=x, y0=CI_95_lower, x1=x, y1=CI_95_upper, col="red")

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

######### store output for later
ellipses.posterior.Q.AllTaxa <- ellipses.posterior
SEA.B.Q.AllTaxa <- SEA.B

## Comparison 0 niche size statistics ----
# Pr = probability HP=holo and pleis, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare credible intervals all taxa across time bins
Pr_P_lt_Q <- sum( SEA.B.PH[,1] < SEA.B.Q.AllTaxa[,1] ) / nrow(SEA.B.Q.AllTaxa) # this works because number of iterations the same  across comparisons
print(Pr_P_lt_Q)

## Comparison 0 niche overlap statistics ----
# No comparisons to make here, simply describing the niche size credible intervals

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes)),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper, c(Pr_P_lt_Q)))
siber_output[7,1:7] <- siber_output_temp

write.csv(siber_output, file="output/siber_stats_table.csv", row.names = F)

# Comparison 1. Each subgroup ----
# use siber.data.all.obj and group.ML.siber.data.all.obj

# This is the code for the initial comparison, which I didn't end up using because the Holocene sample sizes are too small when broken out for each group


## Note, the code chunk below is written generally. Choose the comparison to make in the next lines 
siber.obj.data <- siber.data.all.obj
group.ML.data <- group.ML.siber.all.obj

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.obj.data, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.data), 
                 ylim=c(0, 100),
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.data[3,], col="red", pch = "x", lwd = 2)
# no red x for first becauase not enough samples

# Calculate credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)
CI_95_lower <-c(SEA.B.credibles$V1[2,1], SEA.B.credibles$V2[2,1], SEA.B.credibles$V3[2,1], SEA.B.credibles$V4[2,1]) 
CI_95_upper <- c(SEA.B.credibles$V1[2,2], SEA.B.credibles$V2[2,2], SEA.B.credibles$V3[2,2], SEA.B.credibles$V4[2,2])
x<- c(1,2)
segments(x0=x, y0=CI_95_lower, x1=x, y1=CI_95_upper, col="blue")

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

#print(SEA.B.modes)


## Comparison 1 niche size statistics ----
# Pr = probability HP=holo and pleis, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_HS_lt_HO <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_HS_lt_HO)
Pr_PS_lt_PO <- sum( SEA.B[,3] < SEA.B[,4] ) / nrow(SEA.B)
print(Pr_PS_lt_PO)

# compare same taxa across time bins
Pr_HS_lt_PS <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B)
print(Pr_HS_lt_PS)
Pr_HO_lt_PO <- sum( SEA.B[,2] < SEA.B[,4] ) / nrow(SEA.B)
print(Pr_HO_lt_PO)

# siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
#                                          as.vector(siber.obj.data$sample.sizes)[c(3,1,4,2)],
#                                          as.vector(group.ML.data[2,]), 
#                                          as.vector(group.ML.data[3,]), 
#                                          CI_95_lower, CI_95_upper))
# colnames(siber_output_temp) <- colnames(siber_output)
# siber_output <- rbind(siber_output, siber_output_temp)

## Comparison 1 niche overlap statistics ----
#Pleistocene
obs.overlap.PS.PO <- maxLikOverlap("Pleistocene.Sylvilagus", 
                                   "Pleistocene.Otospermophilus", 
                                   siber.obj.data, p = 0.95, n = 100)
# calculate max likelihood proportion overlaps
prop.overlap.both <- as.numeric(obs.overlap.PS.PO["overlap"] / 
                                  (obs.overlap.PS.PO["area.1"] + obs.overlap.PS.PO["area.2"]))
prop.overlap.PS <- as.numeric(obs.overlap.PS.PO["overlap"] / 
                                obs.overlap.PS.PO["area.1"])
prop.overlap.PO <- as.numeric(obs.overlap.PS.PO["overlap"] / 
                                obs.overlap.PS.PO["area.2"])
obs.overlap.PS.PO <- c(obs.overlap.PS.PO, prop.overlap.both, prop.overlap.PS, prop.overlap.PO)
names(obs.overlap.PS.PO)[4:6] <- c('prop.overlap.both', 'prop.overlap.PS', 'prop.overlap.PO')

# calculate bayes overlap distribution
bayes.overlap.PS.PO <- bayesianOverlap("Pleistocene.Sylvilagus", 
                                       "Pleistocene.Otospermophilus",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)

bayes.prop.overlap.both <- as.numeric(bayes.overlap.PS.PO[,'overlap'] / 
                                        (bayes.overlap.PS.PO[,'area1'] + bayes.overlap.PS.PO[,'area2']))
bayes.prop.overlap.PS <- as.numeric(bayes.overlap.PS.PO[,'overlap'] / bayes.overlap.PS.PO[,'area1'])
bayes.prop.overlap.PO <- as.numeric(bayes.overlap.PS.PO[,'overlap'] / bayes.overlap.PS.PO[,'area2'])
bayes.overlap.PS.PO <- cbind(bayes.overlap.PS.PO, bayes.prop.overlap.both, 
                             bayes.prop.overlap.PS, bayes.prop.overlap.PO)
# save objects
save(obs.overlap.PS.PO, bayes.overlap.PS.PO, file="output/overlap_PS.PO.RData")

#holocene
obs.overlap.HS.HO <- maxLikOverlap("Holocene.Sylvilagus", 
                                   "Holocene.Otospermophilus", 
                                   siber.data.all.obj, p = 0.95, n = 100)
#'Error in verify.xypolygon(P) : x and y coordinates must not contain NA values'
#Nate comments - ran is.na(siber.data.all.obj), all list objects returned FALSE. Not
#sure what the issue is here 
# JLB: Yes, I get the same error - I think it's because you can't make a polygon with only 2 points. Doesn't affect Table 2 (which is written on line 605)

# # calculate max likelihood proportion overlaps
# prop.overlap.both <- as.numeric(obs.overlap.HS.HO["overlap"] / 
#                                   (obs.overlap.HS.HO["area.1"] + obs.overlap.HS.HO["area.2"]))
# prop.overlap.HS <- as.numeric(obs.overlap.HS.HO["overlap"] / 
#                                 obs.overlap.HS.HO["area.1"])
# prop.overlap.HO <- as.numeric(obs.overlap.HS.HO["overlap"] / 
#                                 obs.overlap.HS.HO["area.2"])
# obs.overlap.HS.HO <- c(obs.overlap.HS.HO, prop.overlap.both, prop.overlap.HS, prop.overlap.HO)
# names(obs.overlap.HS.HO)[4:6] <- c('prop.overlap.both', 'prop.overlap.HS', 'prop.overlap.HO')

# calculate bayes overlap distribution
bayes.overlap.HS.HO <- bayesianOverlap("Holocene.Sylvilagus", 
                                       "Holocene.Otospermophilus",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)
bayes.prop.overlap.both <- as.numeric(bayes.overlap.HS.HO[,'overlap'] / 
                                        (bayes.overlap.HS.HO[,'area1'] + bayes.overlap.HS.HO[,'area2']))
bayes.prop.overlap.HS <- as.numeric(bayes.overlap.HS.HO[,'overlap'] / bayes.overlap.HS.HO[,'area1'])
bayes.prop.overlap.HO <- as.numeric(bayes.overlap.HS.HO[,'overlap'] / bayes.overlap.HS.HO[,'area2'])
bayes.overlap.HS.HO <- cbind(bayes.overlap.HS.HO, bayes.prop.overlap.both, 
                             bayes.prop.overlap.HS, bayes.prop.overlap.HO)
# save objects
save(obs.overlap.HS.HO, bayes.overlap.HS.HO, file="output/overlap_HS.HO.RData")

# #sylvilagus
# obs.overlap.PS.HS <- maxLikOverlap("Pleistocene.Sylvilagus", 
#                                    "Holocene.Sylvilagus", 
#                                    siber.data.all.obj, p = 0.95, n = 100)

# # calculate max likelihood proportion overlaps
# prop.overlap.both <- as.numeric(obs.overlap.PS.HS["overlap"] / 
#                                   (obs.overlap.PS.HS["area.1"] + obs.overlap.PS.HS["area.2"]))
# prop.overlap.PS <- as.numeric(obs.overlap.PS.HS["overlap"] / 
#                                 obs.overlap.PS.HS["area.1"])
# prop.overlap.HS <- as.numeric(obs.overlap.PS.HS["overlap"] / 
#                                 obs.overlap.PS.HS["area.2"])
# obs.overlap.PS.HS <- c(obs.overlap.PS.HS, prop.overlap.both, prop.overlap.PS, prop.overlap.HS)
# names(obs.overlap.PS.HS)[4:6] <- c('prop.overlap.both', 'prop.overlap.PS', 'prop.overlap.HS')

# calculate bayes overlap distribution
bayes.overlap.PS.HS <- bayesianOverlap("Pleistocene.Sylvilagus", 
                                       "Holocene.Sylvilagus",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)
bayes.prop.overlap.both <- as.numeric(bayes.overlap.PS.HS[,'overlap'] / 
                                        (bayes.overlap.PS.HS[,'area1'] + bayes.overlap.PS.HS[,'area2']))
bayes.prop.overlap.PS <- as.numeric(bayes.overlap.PS.HS[,'overlap'] / bayes.overlap.PS.HS[,'area1'])
bayes.prop.overlap.HS <- as.numeric(bayes.overlap.PS.HS[,'overlap'] / bayes.overlap.PS.HS[,'area2'])
bayes.overlap.PS.HS <- cbind(bayes.overlap.PS.HS, bayes.prop.overlap.both, 
                             bayes.prop.overlap.PS, bayes.prop.overlap.HS)
# save objects
save(obs.overlap.PS.HS, bayes.overlap.PS.HS, file="output/overlap_PS.HS.RData")


#otospermophilus
obs.overlap.PO.HO <- maxLikOverlap("Pleistocene.Otospermophilus", 
                                   "Holocene.Otospermophilus", 
                                   siber.data.all.obj, p = 0.95, n = 100)

# calculate max likelihood proportion overlaps
prop.overlap.both <- as.numeric(obs.overlap.PO.HO["overlap"] / 
                                  (obs.overlap.PO.HO["area.1"] + obs.overlap.PO.HO["area.2"]))
prop.overlap.PO <- as.numeric(obs.overlap.PO.HO["overlap"] / 
                                obs.overlap.PO.HO["area.1"])
prop.overlap.HO <- as.numeric(obs.overlap.PO.HO["overlap"] / 
                                obs.overlap.PO.HO["area.2"])
obs.overlap.PO.HO <- c(obs.overlap.PO.HO, prop.overlap.both, prop.overlap.PO, prop.overlap.HO)
names(obs.overlap.PO.HO)[4:6] <- c('prop.overlap.both', 'prop.overlap.PO', 'prop.overlap.HO')

# calculate bayes overlap distribution
bayes.overlap.PO.HO <- bayesianOverlap("Pleistocene.Otospermophilus", 
                                       "Holocene.Otospermophilus",
                                       ellipses.posterior, 
                                       draws = 100, p.interval = 0.95,
                                       n = 360)
bayes.prop.overlap.both <- as.numeric(bayes.overlap.PO.HO[,'overlap'] / 
                                        (bayes.overlap.PO.HO[,'area1'] + bayes.overlap.PO.HO[,'area2']))
bayes.prop.overlap.PO <- as.numeric(bayes.overlap.PO.HO[,'overlap'] / bayes.overlap.PO.HO[,'area1'])
bayes.prop.overlap.HO <- as.numeric(bayes.overlap.PO.HO[,'overlap'] / bayes.overlap.PO.HO[,'area2'])
bayes.overlap.PO.HO <- cbind(bayes.overlap.PO.HO, bayes.prop.overlap.both, 
                             bayes.prop.overlap.PO, bayes.prop.overlap.HO)
# save objects
save(obs.overlap.PO.HO, bayes.overlap.PO.HO, file="output/overlap_PO.HO.RData")





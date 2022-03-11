library(tidyverse)
#library(gridExtra)
library(SIBER)
library(ggpubr)

rlbPalette <- palette(c("royalblue2","darkorange"))

# ggplot ellipses ----
# load in the isotope dataset
dat<- read.csv("data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T)

# arrange for ggplot
rlb_data <- dat %>% mutate(Taxon = factor(Taxon), 
                             time_group = factor(time_group),
                             d13C = del13C_permil, 
                             d15N = del15N_permil,
                             .keep = "unused") 

## basic plots for the isotope data ####

rlbPalette <- palette(c("royalblue2","darkorange"))
p.ell <- 0.68

first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(colour = Taxon, shape = time_group), size = 3) +
  scale_colour_manual(labels = c("Otospermophilus", "Sylvilagus"), 
                      values=rlbPalette) +
  scale_shape_manual(labels = c("Holocene", "Pleistocene"), 
                     values=c(16,17)) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme_classic() +
  theme(text = element_text(size=14),
        axis.ticks.length = unit(0.15, "cm")) + 
  labs(colour = "Taxon", shape="Time Period") 
 
print(first.plot) 

# error bars
sbg <- rlb_data %>% 
  group_by(Taxon, time_group) %>% 
  summarise(count = n(),
            mC = mean(d13C), 
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N))

second.plot <- first.plot +
  geom_errorbar(data = sbg, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN, 
                              ymax = mN + 1.96*sdN), 
                width = 0) +
  geom_errorbarh(data = sbg, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0) + 
  geom_point(data = sbg, aes(x = mC, 
                             y = mN,
                             fill = Taxon), 
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values=rlbPalette)

print(second.plot)

ellipse.plot <- first.plot + 
  stat_ellipse(aes(Taxon = interaction(Taxon, time_group), 
                   fill = Taxon, 
                   color = Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

print(ellipse.plot)

## final plots - individual panels ----

# add probability ellipses for all (Pleistocene and Holocene, solid line) and Pleistocene-only (dashed line) data 
# (taxa differentiated)

alldata.taxa.ellipse.plot <- first.plot + 
  stat_ellipse(data = rlb_data, 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[-which(rlb_data$time_group == "Holocene"),], 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               linetype = 2,
               alpha = 0, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

print(alldata.taxa.ellipse.plot)


# add probability ellipses for Pleistocene and Holocene data separately
# (combining taxa)

alldata.time.ellipse.plot <- first.plot + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Pleistocene"), ], 
               aes(x = d13C, 
                   y = d15N), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Holocene"),], 
               aes(x = d13C, 
                   y = d15N), 
               linetype = 2,
               alpha = 0.5, 
               level = p.ell,
               type = "norm",
               geom = "polygon")  

print(alldata.time.ellipse.plot)

## merged - Final Figure 2 plot ----
# will alter in Illustrator afterwards

new.first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  lims(x = c(-23, -17), y= c(0,11)) +
  geom_point(aes(colour = Taxon, shape = time_group), size = 3) +
  scale_colour_manual(labels = c("Otospermophilus", "Sylvilagus"), 
                      values=rlbPalette) +
  scale_shape_manual(labels = c("Holocene", "Pleistocene"), 
                     values=c(16,17)) +
  ylab(expression({delta}^15*N~'value ('~'\u2030'~', AIR)')) +
  xlab(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)')) + 
  theme_classic() +
  theme(text = element_text(size=14),
        axis.ticks.length = unit(0.15, "cm")) + 
  labs(colour = "Taxon", shape="Time Period") 

new.alldata.taxa.ellipse.plot <- new.first.plot + 
  stat_ellipse(data = rlb_data, 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Pleistocene"),], 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               linetype = 2,
               alpha = 0, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

new.alldata.time.ellipse.plot <- new.first.plot + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Pleistocene"),], 
               aes(x = d13C, 
                   y = d15N), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Holocene"),], 
               aes(x = d13C, 
                   y = d15N), 
               linetype = 2,
               alpha = 0.5, 
               level = p.ell,
               type = "norm",
               geom = "polygon") 

figure <- ggarrange(new.alldata.time.ellipse.plot, 
                    new.alldata.taxa.ellipse.plot, 
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

grDevices::cairo_pdf("output/Figure2_SIBERplots_Mar2022_JLB_both.pdf", width=12, height=4)
  print(figure)
dev.off()


# SIBER Stats ----

## Create data frames for SIBER analysis ----
# Note - don't need to plot with SIBER for the paper, just doing this to visualize as I go.
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

# Note: The warning is because sample sizes are too low for Holo for the individual taxa
siber.data.all.obj <- createSiberObject(siber.data.all)

# create new object that groups Pleisto and Holo into one community
siber.data.q <- siber.data.all
siber.data.q$community <- "Quaternary"
siber.data.q.obj <- createSiberObject(siber.data.q)

# Then create separate Pleistocene dataframe
siber.data.p <- siber.data.all[which(siber.data.all$community == "Pleistocene"),]
siber.data.p$community <- droplevels(siber.data.p$community)
siber.data.p.obj <- createSiberObject(siber.data.p)

community.hulls.args <- list(col = "black", lty = 1, lwd = 1) # time group
group.ellipses.args  <- list(n = 100, p.interval = 0.68, lty = 1, lwd = 2) # Sylv and Oto
group.hulls.args     <- list(col="gray", lty = 2) #taxon

# plot data with all four separate groups and communities
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

group.ML.siber.all.obj <- groupMetricsML(siber.data.all.obj)
print(group.ML.siber.all.obj)

# plot data with two taxon groups and 1 Quaternary community
par(mfrow=c(1,1))
plotSiberObject(siber.data.q.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = list(n = 100, p.interval = 0.95, lty = 1, lwd = 2),
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

group.ML.siber.data.q.obj<- groupMetricsML(siber.data.q.obj)
print(group.ML.siber.data.q.obj)

# plot data with two groups and 1 Pleisto community
par(mfrow=c(1,1))
plotSiberObject(siber.data.p.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = list(n = 100, p.interval = 0.95, lty = 1, lwd = 2),
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

group.ML.siber.data.p.obj <- groupMetricsML(siber.data.p.obj)
print(group.ML.siber.data.p.obj)

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

# Note: The warning is because sample sizes are too low for Holo for the individual taxa
siber.data.PH.obj <- createSiberObject(siber.data.PH)

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

group.ML.siber.data.PH.obj<- groupMetricsML(siber.data.PH.obj)
print(group.ML.siber.data.PH.obj)

# create summary table of group statistics
save(group.ML.siber.all.obj,  # data with all four groups
     group.ML.siber.data.p.obj, # data with taxa separate, just for Pleisto
     group.ML.siber.data.q.obj, # data with taxa separate, for Pleisto + Holo
     group.ML.siber.data.PH.obj, # data with taxa merged, for Pleisto and Holo separate
     file="output/siber.summary.stats.RData")

# Calculate SEA statistics ----

# first, compare SIBER output with the standard ellipse stats to make sure plotting is ok
temp<- rlb_data[which(rlb_data$time_group == "Holocene"), c('d13C', 'd15N')]
Sigma <- cov(temp) 
evalues <- eigen(Sigma, symmetric = T, only.values = TRUE)$values
p <- p.ell # make sure you're calculate the same area as plotting
r <- 2 * qf(p, 2, nrow(temp)-1)
a <- sqrt(r * evalues[1]) 
b <- sqrt(r * evalues[2]) # SIBER does not include r values in their calculation
ea_H <- pi*a*b
# compare with:
sigmaSEA(Sigma)
# Note, same, except that SIBER does not include r values in their calculation of a, b, and SEA. So I think we're fine plotting the ellipses, but using the groupML function to quantify SEA and SEAc.

# create a data frame to store siber_output
siber_output <- as.data.frame(matrix(NA, nrow=0, ncol=6))
colnames(siber_output) <- c("community", "n", "SEA", "SEAc", "lower95CI", "upper95CI")


# Comparisons we want to make:
#   1. Is each subgroup significantly different (by taxon AND time)? siber.data.all.obj
#   2. Are Holocene and Pleistocene groups significantly different? siber.data.PH.obj
#   3. Are Oto sig diff from Sylv in Pleistocene? group.ML.siber.data.p.obj
#   4. Are Oto sig diff from Sylv in Quaternary? group.ML.siber.data.q.obj

# For all of these, we will first fit Bayesian multivariate normal distributions to each group in the dataset. Then, we will use siberEllipses to calculate the bayesian Standard Ellipse Area for all groups, then calculate credible intervals and compare groups
# we can compare groups by plottng the groups and the distribution of their ellipses, as well as calculating the probability that one groups posterior distribution is smaller (or larger) than another
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

# Comparison 1. Each subgroup ----
# use siber.data.all.obj and group.ML.siber.data.all.obj

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

# Make comparisons - comment or uncomment depending on which comparison you are making

## Comparison 1 statistics ----
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

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         as.vector(siber.obj.data$sample.sizes)[c(3,1,4,2)],
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper))
colnames(siber_output_temp) <- colnames(siber_output)
siber_output <- rbind(siber_output, siber_output_temp)

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

## Comparison 2 statistics ----
# Pr = probability HP=holo and pleis, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_P_lt_H <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_P_lt_H)

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes)),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper))
colnames(siber_output_temp) <- colnames(siber_output)
siber_output <- rbind(siber_output, siber_output_temp)


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

## Comparison 3 statistics ----
# Pr = probability HP=holo and pleis, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_PS_lt_PO <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_PS_lt_PO)

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes)),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper))
colnames(siber_output_temp) <- colnames(siber_output)
siber_output <- rbind(siber_output, siber_output_temp)


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


## Comparison 4 statistics ----
# Pr = probability HP=holo and pleis, Q= Quaternary, SO = sylv and oto, lt = lower than
# e.g., Prob holocene sylv standard ellipse area is smaller than the SEA of holocene otospermophilus

# compare taxa within time bins
Pr_QS_lt_QO <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pr_QS_lt_QO)

siber_output_temp <- as.data.frame(cbind(as.vector(dimnames(group.ML.data)[[2]]), 
                                         rev(as.vector(siber.obj.data$sample.sizes[2,1])),
                                         as.vector(group.ML.data[2,]), 
                                         as.vector(group.ML.data[3,]), 
                                         CI_95_lower, CI_95_upper))
colnames(siber_output_temp) <- colnames(siber_output)
siber_output <- rbind(siber_output, siber_output_temp)

write.csv(siber_output, file="output/siber_stats_table.csv", row.names = F)











# old ----

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.data.all.obj, parms, priors)

# # extract the posterior means
# mu.post <- extractPosteriorMeans(siber.data.all.obj, ellipses.posterior)

# 
# ----------------------------------------------------------------
# Plot out some of the data and results
# ----------------------------------------------------------------

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.siber.all.obj), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.siber.all.obj[3,], col="red", pch = "x", lwd = 2)


# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)

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


# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

# --------------------------------------
# Visualise the first community
# --------------------------------------

# drop the 3rd column of the posterior which is TA using -3.
siberDensityPlot(layman.B[[1]][ , -3], 
                 xticklabels = colnames(layman.B[[1]][ , -3]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber.data.all.obj$ML.mu[[1]][1,1,], # corresponds to iso1
                                 siber.data.all.obj$ML.mu[[1]][1,2,] # corresponds to iso2
)

# again drop the 3rd entry which relates to TA
points(1:5, comm1.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)






# --------------------------------------
# Visualise the second community
# --------------------------------------
siberDensityPlot(layman.B[[2]][ , -3], 
                 xticklabels = colnames(layman.B[[2]][ , -3]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(siber.data.all.obj$ML.mu[[2]][1,1,],
                                 siber.data.all.obj$ML.mu[[2]][1,2,]
)
points(1:5, comm2.layman.ml$metrics[-3], 
       col = "red", pch = "x", lwd = 2)

--------------------------------------
  # Alternatively, pull out TA from both and aggregate them into a 
  # single matrix using cbind() and plot them together on one graph.
  # --------------------------------------

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

# Now we only plot the TA data. We could address this as either
# layman.B[[1]][, "TA"]
# or
# layman.B[[1]][, 3]
siberDensityPlot(cbind(layman.B[[1]][ , "TA"], 
                       layman.B[[2]][ , "TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0, 90),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML.siber.data.q.obj), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML.siber.data.q.obj[3,], col="red", pch = "x", lwd = 2)


# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

# extract the posterior means
mu.post <- extractPosteriorMeans(siber.data.q.obj, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

# --------------------------------------
# Visualise the first community
# --------------------------------------
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber.data.q.obj$ML.mu[[1]][1,1,],
                                 siber.data.q.obj$ML.mu[[1]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)



# extract the posterior means
mu.post <- extractPosteriorMeans(siber.example, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)


# --------------------------------------
# Visualise the second community
# --------------------------------------
siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(siber.RLB$ML.mu[[2]][1,1,],
                                 siber.RLB$ML.mu[[2]][1,2,]
)
points(1:6, comm2.layman.ml$metrics, col = "red", pch = "x", lwd = 2)



# --------------------------------------
# Alternatively, pull out TA from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.
# --------------------------------------

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0,20),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")

# default SIBER plots - using SIBER plotting ####

# change ellipse color
palette(c("royalblue2","darkorange"))

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.68, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20") # change color here to be organe vs blue? 

#pdf("output/isotope paper final/Figure2_SIBER.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(5,5,4,1)+0.01)
plotSiberObject(siber.RLB,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
                
)

legend("bottomleft", legend = c("Otospermophilus Pre-LGM", "Sylvilagus Pre-LGM"),
       col = palette(c("royalblue2","darkorange")), pch = 1, 
       bty = "n", cex = 0.8)

legend("topright", legend = c("Otospermophilus Post-LGM", "Sylvilagus Post-LGM"),
       col = palette(c("royalblue2", "darkorange")), pch = 2,
       bty = "n", cex = 0.8)
dev.off()


# MANOVA stats ----

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


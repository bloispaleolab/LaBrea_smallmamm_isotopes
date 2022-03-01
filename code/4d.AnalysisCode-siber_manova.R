library(tidyverse)
#library(gridExtra)
library(SIBER)
library(ggpubr)

rlbPalette <- palette(c("royalblue2","darkorange"))

# load in the isotope dataset
dat<- read.csv("data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T)

# arrange for ggplot
rlb_data <- dat %>% mutate(Taxon = factor(Taxon), 
                             time_group = factor(time_group),
                             d13C = del13C_permil, 
                             d15N = del15N_permil,
                             .keep = "unused") 

# basic plots for the isotope data ####

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

# final plots - individual panels ----

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

# merged - Final Figure 2 plot ----
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


#SIBER Stats ----

## Create data frames for SIBER analysis
# First create a df with both groups (Sylv and Oto) and communities (Holo vs Pleisto). 
siber.data.all <- dat %>% 
  rename(iso1 = del13C_permil,
        iso2 = del15N_permil,
        community = time_group,
        group = Taxon) %>%
    select(iso1, iso2, group, community)

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
group.hulls.args     <- list(col="gray", lty = 2) #taxon
group.ellipses.args  <- list(n = 100, p.interval = 0.68, lty = 1, lwd = 2) # Sylv and Oto

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

group.ML.siber.obj <- groupMetricsML(siber.data.all.obj)
print(group.ML.siber.obj)

# plot data with two groups and 1 Quaternary community
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


## JESSICA -finished code through here.

# fit Bayesian multivariate normal distributions to each group in the dataset. do this for all four separately ----
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

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.obj, parms, priors)


# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


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
mu.post <- extractPosteriorMeans(siber.RLB, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

# --------------------------------------
# Visualise the first community
# --------------------------------------
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber.RLB$ML.mu[[1]][1,1,],
                                 siber.RLB$ML.mu[[1]][1,2,]
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


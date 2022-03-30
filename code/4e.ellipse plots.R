library(tidyverse)
#library(gridExtra)
library(ggpubr)

# ggplot ellipses ----
# load in the isotope dataset
dat<- read.csv("data/processed/final_dataset_focaltaxa_with_calages_climate.csv", header=T)

# arrange for ggplot
rlb_data <- dat %>% mutate(Taxon = factor(Taxon), 
                           time_group = factor(time_group),
                           d13C = del13C_permil, 
                           d15N = del15N_permil,
                           .keep = "unused") 

## Final Figure 2 plots ----
# will alter in Illustrator afterwards

rlbPalette <- palette(c("royalblue2","darkorange"))
p.ell <- 0.68

# lay out basic plot structure
base.plot <- ggplot(data = rlb_data, 
                         aes(x = d13C, 
                             y = d15N)) + 
  lims(x = c(-23, -17), y= c(0.5,11)) +
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

# Figure 2a
# add probability ellipses for all (Quaternary, solid line, light ellipse) and Pleistocene-only (dashed line) data 
# (taxa differentiated)

new.alldata.taxa.ellipse.plot <- base.plot + 
  stat_ellipse(data = rlb_data, 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               linetype = 1,
               alpha = 0.1, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Pleistocene"),], 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               linetype = 2,
               alpha = 0.3, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

# 2b
# add probability ellipses for Quaternary, Pleistocene and Holocene data separately
# (combining taxa)

new.alldata.time.ellipse.plot <- base.plot + 
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Holocene"),], 
               color="black",
               linetype = 3,
               alpha = 0.1, 
               level = p.ell,
               type = "norm",
               geom = "polygon") +
  stat_ellipse(data = rlb_data, 
               color="black", 
               linetype = 1,
               alpha = 0.3, 
               level = p.ell,
               type = "norm",
               geom = "polygon") +
  stat_ellipse(data = rlb_data[which(rlb_data$time_group == "Pleistocene"),],
               color="black",
               linetype = 2,
               alpha = 0.5, 
               level = p.ell,
               type = "norm",
               geom = "polygon")

figure <- ggarrange(new.alldata.taxa.ellipse.plot,
                    new.alldata.time.ellipse.plot, 
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
pdf(file="output/Figure2_SIBERplots_Mar2022_NF.pdf", height=4, width=12)
#grDevices::cairo_pdf("output/Figure2_SIBERplots_Mar2022_JLB_both.pdf", width=12, height=4)
print(figure)
dev.off()



# other plots ----
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
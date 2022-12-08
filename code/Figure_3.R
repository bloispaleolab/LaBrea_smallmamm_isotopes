library(tidyverse)
library(readr)
library(tidyverse)
library(readr)

### Final Figure 3 code ####

# plot the isotope data - using ggplot2 ####
fig3Dat<- read_csv(file="data/processed/final_dataset_Figure3.csv", col_names=T, trim_ws=T)

fig3Dat$Taxon <- as.factor(fig3Dat$Taxon)
fig3Dat$Time <- as.factor(fig3Dat$Time)
fig3Dat$Group <- as.factor(fig3Dat$Group)

fig3Dat<- fig3Dat %>%
  arrange(Group) %>%
  mutate(Taxon = fct_relevel(Taxon, 
                             "Bison", "Camelops", "Equus", "Mammut", "Paramylodon", "Neotoma", "Otospermophilus", "Sylvilagus", "Thomomys"), Group = fct_relevel(Group, "Megaherbivore", "Herbivore"), Time = fct_relevel(Time, "Holocene", "Pleistocene"))

fig3dat_mega <-  filter(fig3Dat, Group =="Megaherbivore")
fig3dat_micro <-  filter(fig3Dat, Group =="Herbivore")
fig3dat_micro <- droplevels(fig3dat_micro)

p.ell <- 0.68

labelColors <- c("Bison" = "gray30",
                 "Camelops" = "gray70",
                 "Equus" = "gray50", 
                 "Mammut" = "gray10",
                 "Paramylodon" = "gray90",
                 "Neotoma" = "forestgreen",
                 "Otospermophilus" = "royalblue2",
                 "Sylvilagus" = "darkorange",
                 "Thomomys" = "red3")

lineTypes <- c("Bison" = 2,
                 "Camelops" = 3,
                 "Equus" = 4, 
                 "Mammut" = 5,
                 "Paramylodon" = 6,
                 "Neotoma" = 1,
                 "Otospermophilus" = 1,
                 "Sylvilagus" = 1,
                 "Thomomys" = 1)

# error bars
sbg <- fig3dat_micro %>% 
  group_by(Taxon, Time) %>% 
  summarise(count = n(),
            mC = mean(del13C_permil), 
            sdC = sd(del13C_permil), 
            mN = mean(del15N_permil), 
            sdN = sd(del15N_permil))
sbg.shape.code <- c(24, 21)
sbg.color.code <- c("forestgreen", "royalblue2", "royalblue2", "darkorange", "darkorange", "red3")


first.plot_new <- ggplot(data=fig3Dat, aes(x=del13C_permil, y=del15N_permil, color=Taxon, shape=Time)) +
  ylab(expression({delta}^15*N~'value ('~'\u2030'~', AIR)')) +
  xlab(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)')) + 
  theme_classic() +
  theme(text = element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size=12), 
        legend.key.size = unit(0.4, 'cm')) 
print(first.plot_new)

second.plot_new <- first.plot_new + 
  stat_ellipse(data = fig3dat_mega, 
             aes(x = del13C_permil, 
                 y = del15N_permil, color=Taxon, linetype=Taxon, fill=Taxon),
             alpha = 0.8, level = p.ell, type = "norm", geom = "polygon") +
  geom_point(data = fig3dat_micro,
             aes(colour = Taxon, shape = Time), alpha=0.5) +
  geom_errorbar(data = sbg, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN, 
                              ymax = mN + 1.96*sdN), 
                width = 0, alpha=1) +
  geom_errorbarh(data = sbg, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0, alpha=1) + 
  geom_point(data = sbg, aes(x = mC, 
                             y = mN,
                             fill = Taxon, shape = Time), 
             color = sbg.color.code, size = 4,
             alpha = 1, show.legend = TRUE)+
  scale_fill_manual(values=labelColors) +
  scale_linetype_manual(values = lineTypes ) +
  scale_color_manual(values=labelColors)
#grDevices::cairo_pdf(file="output/Figure3_v4_megaherbivore_ellipses.pdf", height=4, width=6)
pdf(file="output/Figure3_v4_megaherbivore_ellipses.pdf", height=4, width=6)
print(second.plot_new) 
dev.off()

second.plot_new2 <- first.plot_new + 
  stat_ellipse(data = fig3dat_mega, 
               aes(x = del13C_permil, 
                   y = del15N_permil, color=Taxon, linetype=Taxon, fill=Taxon),
               alpha = 0.8, level = p.ell, type = "norm", geom = "polygon") +
  stat_ellipse(data = fig3dat_micro, 
               aes(x = del13C_permil, 
                   y = del15N_permil, color=Taxon, linetype=Taxon, fill=Taxon),
               alpha = 0.5, level = p.ell, type = "norm", geom = "polygon")+
  scale_fill_manual(values=labelColors) +
  scale_linetype_manual(values = lineTypes ) +
  scale_color_manual(values=labelColors)+
  geom_point(data = fig3dat_micro,
             aes(colour = Taxon, shape = Time), alpha=1, size=3) 
 
grDevices::cairo_pdf(file="output/Figure3_v4_megaherbivore_ellipses.pdf", height=4, width=6)
#pdf(file="output/Figure3_v4_megaherbivore_ALLellipses.pdf", height=4, width=6)
print(second.plot_new2) 
dev.off()

# note: cairo_pdf gets the axes labels right, but the linetype symbology is not as differentiated. Used cairo_pdf to generate axes to copy over, but used pdf for main plot.
# changed the Taxon legend in Illustrator to remove circles for megafauna and make the colors more clear, plus a few other legend changes

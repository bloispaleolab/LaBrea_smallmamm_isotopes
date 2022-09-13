library(tidyverse)
library(readr)
library(tidyverse)
library(readr)

# plot the isotope data - using ggplot2 ####
fig3Dat<- read_csv(file="data/processed/final_dataset_Figure3.csv", col_names=T, trim_ws=T)

fig3Dat$Taxon <- as.factor(fig3Dat$Taxon)
fig3Dat$Time <- as.factor(fig3Dat$Time)
fig3Dat$Group <- as.factor(fig3Dat$Group)

fig3Dat<- fig3Dat %>%
  arrange(Group) %>%
  mutate(Taxon = fct_relevel(Taxon, 
                             "Bison", "Camelops", "Equus", "Mammut", "Paramylodon", "Neotoma", "Otospermophilus", "Sylvilagus", "Thomomys"), Group = fct_relevel(Group, "Megaherbivore", "Herbivore"), Time = fct_relevel(Time, "Holocene", "Pleistocene"))

labelColors <- c("Bison" = "gray20",
                 "Camelops" = "gray30",
                 "Equus" = "gray40", 
                 "Mammut" = "gray50",
                 "Paramylodon" = "gray60",
                 "Neotoma" = "forestgreen",
                 "Otospermophilus" = "royalblue2",
                 "Sylvilagus" = "darkorange",
                 "Thomomys" = "red3")
# error bars
sbg <- fig3Dat %>% 
  group_by(Taxon, Time) %>% 
  summarise(count = n(),
            mC = mean(del13C_permil), 
            sdC = sd(del13C_permil), 
            mN = mean(del15N_permil), 
            sdN = sd(del15N_permil))
sbg.shape.code <- c(22,22,22,22,22,22, 21, 22, 21, 22, 22)
sbg.color.code <- c("gray20", "gray30", "gray40", "gray50", "gray60", "forestgreen", "royalblue2", "royalblue2", "darkorange", "darkorange", "red3")


first.plot <- ggplot(data = fig3Dat, 
                     aes(x = del13C_permil, y = del15N_permil)) +
                       geom_point(aes(colour = Taxon, 
                                      shape = Time), alpha=0.5) +
                       scale_size_manual(values = c(1, 2)) +
                       scale_color_manual(values=labelColors) +
  ylab(expression({delta}^15*N~'value ('~'\u2030'~', AIR)')) +
  xlab(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)')) + 
  theme_classic() +
  theme(text = element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size=12), 
        legend.key.size = unit(0.4, 'cm')) + 
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(colour = "Taxon") 

  print(first.plot) 

second.plot <- first.plot+
  geom_errorbar(data = sbg, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN, 
                              ymax = mN + 1.96*sdN), 
                width = 0, color=sbg.color.code, alpha=1) +
  geom_errorbarh(data = sbg, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0, color=sbg.color.code, alpha=1) + 
  geom_point(data = sbg, aes(x = mC, 
                             y = mN,
                             fill = Taxon), 
             color = "black", shape = sbg.shape.code, size = 4,
             alpha = 1, show.legend = FALSE) +
  scale_fill_manual(values=labelColors)

grDevices::cairo_pdf("output/Figure3_v3_witherrorbars_Sep2022.pdf", width=6, height=4)
#pdf(file="output/Figure3_v2_witherrorbars-NF.pdf", height=4, width=6)
  print(second.plot) 
dev.off()

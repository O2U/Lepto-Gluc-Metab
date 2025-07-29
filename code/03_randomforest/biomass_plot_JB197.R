library(tidyverse)
library(ggsignif)

# read data
jb197.29 <- read_delim("data/processed/RIPTiDe/jb197_29_20240227_130901/flux_samples.tsv") %>% select(-1)
jb197.37 <- read_delim("data/processed/RIPTiDe/jb197_37_20240227_130922/flux_samples.tsv") %>% select(-1)

# biomass flux
bm.jb197.29 <- jb197.29$biomass_GmNeg
bm.jb197.37 <- jb197.37$biomass_GmNeg

summary(bm.jb197.29)
summary(bm.jb197.37)

df.bm <- data.frame(
  "t29" <- c(bm.jb197.29),
  "t37" <- c(bm.jb197.37)
)
colnames(df.bm) <- c("29°C", "37°C")

head(df.bm)

df.bm.l = pivot_longer(df.bm, cols = c("29°C", "37°C"), names_to = "Temp", values_to = 'biomass')
head(df.bm.l)

g <- ggplot(df.bm.l, aes(x = Temp, y = biomass, fill = Temp)) +
  geom_violin() +
  scale_fill_manual(values = c("Grey", "Red")) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 500)) +
  geom_signif(comparisons = list(c("29°C", "37°C")),
              test = "wilcox.test", na.rm = FALSE, map_signif_level = TRUE) +
  xlab(NULL) +
  ylab("Sampled biomass flux") +
  guides(fill = "none") +
  theme()
g

ggsave("results/figures/Fig2/biomass_plot_JB197.pdf", g, width = 2, height = 2)

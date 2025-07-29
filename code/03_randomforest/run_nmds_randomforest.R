library(tidyverse)
library(vegan)
library(randomForest)

# read data
jb197.29 <- read_delim("data/processed/RIPTiDe/jb197_29_20240227_130901/flux_samples.tsv") %>% select(-1)
jb197.37 <- read_delim("data/processed/RIPTiDe/jb197_37_20240227_130922/flux_samples.tsv") %>% select(-1)

# Whole distribution analysis
overlap <- intersect(colnames(jb197.29), colnames(jb197.37))
only29 <- setdiff(colnames(jb197.29), colnames(jb197.37))
only37 <- setdiff(colnames(jb197.37), colnames(jb197.29))
only29 <- jb197.29[, only29]
jb197.29 <- jb197.29[, overlap]
only37 <- jb197.37[, only37]
jb197.37 <- jb197.37[, overlap]

# Format row names
names_29 <- paste('temp29_', 1:nrow(jb197.29), sep='')
rownames(jb197.29) <- names_29
names_37 <- paste('temp37_', 1:nrow(jb197.37), sep='')
rownames(jb197.37) <- names_37

# Create metadata
metadata29 <- cbind(names_29, rep('29°C', length(names_29)))
metadata37 <- cbind(names_37, rep('37°C', length(names_37)))
metadata <- rbind(metadata29, metadata37)
colnames(metadata) <- c('label', 'group')
metadata <- as.data.frame(metadata)
write.csv(metadata, "data/processed/RandomForest/JB197_metadata.csv")

# Merge data and prep for unsupervised learning
all_samples <- rbind(jb197.29, jb197.37)
all_samples <- all_samples + abs(min(all_samples))
all_samples[1:10, 1:10]
write.csv(all_samples, "data/processed/RandomForest/JB197_all_samples.csv")

flux_dist <- vegdist(all_samples, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
flux_groups <- as.factor(c(rep('29°C',nrow(jb197.29)), rep('37°C',nrow(jb197.37))))
meandist(flux_dist, grouping=flux_groups)

# Unsupervised learning
set.seed(412)
fn <- metaMDS(flux_dist, k=2, trymax=50)
flux_nmds <- as.data.frame(fn$points)
write.csv(flux_nmds, "data/processed/RandomForest/JB197_nmds.csv")

# read data
flux_nmds <- as.data.frame(read.csv("data/processed/RandomForest/JB197_nmds.csv", row.names = 1))

# NMDS
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

# Subset axes
solo_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% names_29)
coculture_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% names_37)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
pval <- adonis(flux_dist ~ group, data=test, perm=999, method='bray')
pval <- round(pval$aov.tab[[6]][1], 4)
if (pval > 0.05) {pval <- 'n.s.'} else {pval <- as.character(pval)}

p <- paste("P = ", pval, sep = "")

head(metadata)
flux_nmds$Temp <- metadata$group
head(flux_nmds)

g2 <- ggplot() +
  geom_point(data=flux_nmds, aes(x=MDS1, y=MDS2, colour=Temp), alpha = 0.2, size = 2) +
  scale_colour_manual(values=c("Black", "Red")) + #色を指定
  theme_classic() +
  xlab("NMDS Axis 1") +
  ylab("NMDS Axis 2") +
  scale_x_continuous(limits = c(-0.015, 0.015)) +
  scale_y_continuous(limits = c(-0.015, 0.015)) +
  theme(
    legend.position = c(.15, .90),
    legend.background = element_blank(),
    legend.title = element_blank(),
  ) +
  annotate("text", x = -0.001, y = 0.0127, label = p, size = 2) +
  annotate("segment", x = -0.004, xend = -0.004, y = 0.0143, yend = 0.0111, linewidth = 0.2) +
  annotate("segment", x = -0.004, xend = -0.005, y = 0.0143, yend = 0.0143, linewidth = 0.2) +
  annotate("segment", x = -0.004, xend = -0.005, y = 0.0111, yend = 0.0111, linewidth = 0.2) +
  labs(colour="Temp") #凡例の名前をきちんとする
g2

ggsave("results/figures/Fig2/NMDS_plot_JB197.pdf", g2, width = 3, height = 3)


# RandomForest
df.rf <- cbind(metadata$group, all_samples)
head(df.rf)
colnames(df.rf)[1] = "group"
df.rf$group <- factor(df.rf$group)
rxns <- colnames(all_samples)[1]
for (i in 2:ncol(all_samples)) {
  rxns <- paste(rxns, colnames(all_samples)[i], sep='+')
}
fm <- as.formula(paste('group ~ ', rxns, sep=''))
set.seed(12)
temp_vs_rxn <- randomForest(formula = fm, data = df.rf, importance = TRUE, ntree = 500)
varImpPlot(temp_vs_rxn)

# 3. 特徴重要度の抽出
importance_values <- importance(temp_vs_rxn)
importance_df <- data.frame(
  Feature = rownames(importance_values),
  MeanDecreaseAccuracy = importance_values[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = importance_values[, "MeanDecreaseGini"]
)

varImpPlot(temp_vs_rxn, sort = TRUE)

# 重要度でソート
importance_df <- importance_df[order(-importance_df$MeanDecreaseAccuracy), ]

# 5. ファイルに保存
write.table(importance_df, file="data/processed/RandomForest/RandomForest_JB197_Reactions.tsv", sep="\t", quote=FALSE, row.names=FALSE)

imp <- head(importance_df, 10)


# ggplot2でプロット
mda_plot <- ggplot(imp, aes(y = reorder(Feature, MeanDecreaseAccuracy), x = MeanDecreaseAccuracy)) +
  geom_point(shape = 21, fill = "#999999", colour = "#000000", size = 3) +
  #geom_vline(xintercept = 4, linetype = "dashed", color = "black", linewidth = 0.3) +
  theme_classic() +
  xlab("Mean decreased accuracy (%)") +
  ylab(NULL) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey60', linetype = 'dotted'),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 9)
  ) +
  xlim(c(0, 6))

mda_plot

# 8. 結果をPDFに保存する場合
ggsave('results/figures/Fig2/RandomForest_JB197_reactions.pdf',mda_plot, width=3, height=3)

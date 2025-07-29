library(tidyverse)
library(ggsignif)
library(patchwork)
library(rstatix)

# read data
importance_df <- read_delim("data/processed/RandomForest/RandomForest_JB197_Reactions.tsv", delim = "\t")
df_plot <- head(importance_df, 10)

# read plot data
jb197.29 <- read_delim("data/processed/RIPTiDe/jb197_29_20240227_130901/flux_samples.tsv") %>% select(-1)
jb197.37 <- read_delim("data/processed/RIPTiDe/jb197_37_20240227_130922/flux_samples.tsv") %>% select(-1)
jb197.29$type <- "29°C"
jb197.37$type <- "37°C"

# Whole distribution analysis
overlap <- intersect(colnames(jb197.29), colnames(jb197.37))
only29 <- setdiff(colnames(jb197.29), colnames(jb197.37))
only37 <- setdiff(colnames(jb197.37), colnames(jb197.29))
only29 <- jb197.29[, only29]
jb197.29c <- jb197.29[, overlap]
only37 <- jb197.37[, only37]
jb197.37c <- jb197.37[, overlap]

# データフレームを結合
plot_samples <- rbind(jb197.29c, jb197.37c)

# 個別の箱ひげ図
# ggsignifを使った検定結果を表示する関数
plot_reaction_boxplot_sign <- function(data, reaction, ymin = -1000, ymax = 1000) {
  # チェック: 指定した反応がデータに存在するか確認
  if (!reaction %in% colnames(data)) {
    stop(paste("反応", reaction, "はデータ内に存在しません。"))
  }
  
  # チェック: 'type' 列が存在するか
  if (!"type" %in% colnames(data)) {
    stop("データに 'type' 列が存在しません。正しいデータを渡してください。")
  }
  
  # データフレームを明示的に作成
  plot_data <- data.frame(
    ReactionValue = as.numeric(data[[reaction]]),
    Condition = as.factor(data$type),
    stringsAsFactors = FALSE
  )
  
  # Dunn検定（多重比較）を実行
  dunn_test <- plot_data %>%
    dunn_test(ReactionValue ~ Condition, p.adjust.method = "bonferroni")
  
  # p値をアスタリスクに変換する関数
  pval_to_stars <- function(p) {
    if (p < 0.0001) return("****")
    else if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  
  # 検定結果をggsignif用に整形
  comparison_list <- list(
    c("29°C", "37°C")
  )
  star_annotations <- sapply(dunn_test$p.adj, pval_to_stars)  # p値をアスタリスクに変換
  
  # 箱ひげ図の作成
  p <- ggplot(plot_data, aes(x = Condition, y = ReactionValue, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +  # 箱ひげ図
    #geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +  # データ点の散布
    geom_signif(
      comparisons = comparison_list,
      annotations = star_annotations,  # アスタリスクを表示
      tip_length = 0.01,
      y_position = ymax * 0.8 - 10,  # p値の表示位置を調整
      size = 0.8
    ) +
    labs(
      title = reaction,
      x = NULL,
      y = "Sampled Flux"
    ) +
    theme_classic() +
    scale_fill_manual(values = c("29°C" = "white", "37°C" = "red")) +
    ylim(c(ymin, ymax)) +
    geom_hline(yintercept = 0,             # y = 0 の位置に水平線
               linetype = "dashed",        # 点線の指定
               color = "black",              # 線の色
               linewidth = 0.2) +               # 線の太さ
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p)
}


rxn_plot_only37 <- function(i, amin, amax, title) {
  # 必要な列の名前を取得
  rxn <- colnames(only37)[i]
  df37 <- only37 %>% select(rxn)
  
  # データフレーム作成
  df <- data.frame(
    "29°C" <- c(-1000),
    "37°C" <- c(df37)
  )
  colnames(df) <- c("29°C", "37°C")
  df.l <- pivot_longer(df, cols = c("29°C", "37°C"), names_to = "group", values_to = 'flux')
  df.l$group <- factor(df.l$group, levels = c("29°C", "37°C"))
  
  # ggplot作成
  g <- ggplot(df.l, aes(x = group, y = flux)) +
    # "DMC" の箱ひげ図
    geom_boxplot(data = df.l %>% filter(group == "37°C"),
                 aes(fill = group)) +
    # "In vitro" の列にラベル "Inactivated" を表示
    geom_text(data = df.l %>% filter(group == "29°C"),
              aes(x = group, y = amin + (amax - amin) * 0.1, label = "Inactivated"),
              color = "black", size = 2.5, hjust = 0.5, family = "Helvetica") +
    # カスタム設定
    scale_fill_manual(values = c("red")) +
    theme_classic() +
    # geom_signif(comparisons = list(c("In vitro", "DMC")),
    #             test = "wilcox.test", na.rm = FALSE, map_signif_level = TRUE) +
    ggtitle(title) +
    ylab("Sampled flux") +
    scale_y_continuous(limits = c(amin, amax)) +
    scale_x_discrete(limits = c("29°C", "37°C")) + # 表示順を指定
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10)
    ) +
    geom_hline(yintercept = 0,             # y = 0 の位置に水平線
               linetype = "dashed",        # 点線の指定
               color = "black",              # 線の色
               linewidth = 0.2) +               # 線の太さ
    guides(fill = "none")
  
  return(g)
}


indv <- plot_reaction_boxplot_sign(plot_samples, df_plot$Feature[3], -150, 30) +
  plot_reaction_boxplot_sign(plot_samples, df_plot$Feature[8], 0, 160) +
  plot_reaction_boxplot_sign(plot_samples, df_plot$Feature[9], 0, 160) +
  rxn_plot_only37(1, 0, 160, colnames(only37)[1]) +
  rxn_plot_only37(2, 0, 160, colnames(only37)[2]) +
  rxn_plot_only37(3, 0, 60, colnames(only37)[3]) +
  rxn_plot_only37(4, 0, 160, colnames(only37)[4]) +
  plot_layout(ncol = 7)

indv

# You can save this as PDF
# ggsave('rf_extraction_Temp_plots_focused250214.pdf',indv, width=14, height=2)


# 個別の結果をPDFに保存
ggsave(paste0("results/figures/Fig2/", df_plot$Feature[3], ".pdf"), plot_reaction_boxplot_sign(plot_samples, df_plot$Feature[3], -150, 30), width=2, height=2)
ggsave(paste0("results/figures/Fig2/", df_plot$Feature[8], ".pdf"), plot_reaction_boxplot_sign(plot_samples, df_plot$Feature[8], 0, 160), width=2, height=2)
ggsave(paste0("results/figures/Fig2/", df_plot$Feature[9], ".pdf"), plot_reaction_boxplot_sign(plot_samples, df_plot$Feature[9], 0, 160), width=2, height=2)
ggsave(paste0("results/figures/Fig2/", colnames(only37)[1], ".pdf"), rxn_plot_only37(1, 0, 160, colnames(only37)[1]), width=2, height=2)
ggsave(paste0("results/figures/Fig2/", colnames(only37)[2], ".pdf"), rxn_plot_only37(2, 0, 160, colnames(only37)[2]), width=2, height=2)
ggsave(paste0("results/figures/Fig2/", colnames(only37)[3], ".pdf"), rxn_plot_only37(3, 0, 60, colnames(only37)[3]), width=2, height=2)
ggsave(paste0("results/figures/Fig2/", colnames(only37)[4], ".pdf"), rxn_plot_only37(4, 0, 160, colnames(only37)[4]), width=2, height=2)

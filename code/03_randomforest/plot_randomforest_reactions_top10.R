library(tidyverse)
library(ggsignif)
library(patchwork)
library(rstatix)

importance_df <- read_delim("data/processed/RandomForest/RandomForest_JB197_Reactions.tsv", delim = "\t")
df_plot <- head(importance_df, 10)

# ggsignifを使った検定結果を表示する関数
plot_reaction_boxplot <- function(data, reaction, ymin = -1000, ymax = 1000) {
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
    # geom_signif(
    #   comparisons = comparison_list,
    #   annotations = star_annotations,  # アスタリスクを表示
    #   tip_length = 0.01,
    #   y_position = c(ymax * 0.8, ymax * 0.72, ymax * 0.85),  # p値の表示位置を調整
    #   size = 0.8
    # ) +
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

plot_reaction_boxplot(plot_samples, df_plot$Feature[1])


# すべての反応に対するプロットを作成し、patchwork で結合
plots <- list()  # プロットを格納するリスト

for (reaction in df_plot$Feature) {
  if (reaction %in% colnames(plot_samples)) {
    # 反応データの最小・最大値を取得
    reaction_values <- as.numeric(plot_samples[[reaction]])
    ymin <- min(reaction_values, na.rm = TRUE) * 1.2
    ymax <- max(reaction_values, na.rm = TRUE) * 1.2
    
    # プロットを作成
    p <- plot_reaction_boxplot(plot_samples, reaction, ymin, ymax)
    plots[[reaction]] <- p
  } else {
    message(paste("警告: 反応", reaction, "はデータ内に存在しません。"))
  }
}

# すべてのプロットをpatchworkで結合して表示
if (length(plots) > 0) {
  combined_plot <- patchwork::wrap_plots(plots, ncol = 5)  # 2列で表示
  print(combined_plot)
} else {
  message("表示するプロットがありません。")
}

# 結果をPDFに保存
ggsave('results/figures/FigS1/RandomForest_JB197_rections_plots.pdf',combined_plot, width=10, height=4)

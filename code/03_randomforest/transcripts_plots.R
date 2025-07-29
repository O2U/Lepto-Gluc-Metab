library(tidyverse)
library(ggsignif)
library(patchwork)


# GPRに基づく発現量変化
tr37 <- read_delim("data/processed/RIPTiDe/jb197_37_20240227_130922/transcriptome.tsv", delim = '\t', col_names = NA)
tr29 <- read_delim("data/processed/RIPTiDe/jb197_29_20240227_130901/transcriptome.tsv", delim = '\t', col_names = NA)

tr <- tr37 %>% left_join(tr29, by = "X1")
head(tr)
colnames(tr) <- c("ID", "37A", "37B", "37C", "37D", "29A", "29B", "29BC", "29D")


plot_gene_expression <- function(tr, target_id, title = target_id) {
  # 対象IDでフィルタ
  gene_data <- tr %>% filter(ID == target_id)
  
  if (nrow(gene_data) == 0) {
    stop("指定したIDが見つかりません")
  }
  
  # データをロング形式に変換して、グループ列を追加
  df_long <- gene_data %>%
    pivot_longer(cols = -ID, names_to = "Sample", values_to = "Expression") %>%
    mutate(Group = ifelse(str_detect(Sample, "^37"), "37°C", "29°C"))
  df_long$Sample <- factor(df_long$Sample, levels = c("29°C", "37°C"))
  
  # 各グループのデータ
  group_29 <- df_long %>% filter(Group == "29°C") %>% pull(Expression)
  group_37 <- df_long %>% filter(Group == "37°C") %>% pull(Expression)
  
  # Shapiro-Wilk 正規性検定
  sw_29 <- shapiro.test(group_29)
  sw_37 <- shapiro.test(group_37)
  
  # # 両方とも正規分布なら t検定、そうでなければWilcoxon検定
  # if (sw_29$p.value > 0.05 && sw_37$p.value > 0.05) {
  #   test_used <- "t.test"
  #   test_result <- t.test(Expression ~ Group, data = df_long)
  # } else {
  #   test_used <- "wilcox.test"
  #   test_result <- wilcox.test(Expression ~ Group, data = df_long)
  # }
  
  test_used <- "wilcox.test"
  test_result <- wilcox.test(Expression ~ Group, data = df_long)
  
  # 検定結果のp値をもとにアノテーションを決定
  p_value <- test_result$p.value
  annotation <- ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(p_value < 0.05, "*", "ns")))
  
  # グループごとの平均値を計算
  means <- df_long %>%
    group_by(Group) %>%
    summarise(Mean = mean(Expression), .groups = "drop")
  
  # Y軸の上限値
  y_max <- max(df_long$Expression) * 1.2
  
  # ドットプロットを作成（binaxis = "y"により縦ドット）
  ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2) +
    geom_errorbar(data = means, aes(x = Group, ymin = Mean, ymax = Mean), 
                  width = 0.5, color = "black", inherit.aes = FALSE) +
    theme_classic() +
    ylim(c(0, y_max)) +
    scale_fill_manual(values = c("grey80", "red")) +
    geom_signif(comparisons = list(c("29°C", "37°C")),
                annotations = annotation,
                y_position = max(df_long$Expression) * 1.05,
                tip_length = 0.01) +
    labs(title = paste(title),
         subtitle = NULL, #paste(test_used),
         x = NULL, y = "Gene Expression") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}


plt <-
  
  # Glucose
  plot_gene_expression(tr, "lbj:LBJ_2173") +
  
  # G6P
  plot_gene_expression(tr, "lbj:LBJ_1861") +
  plot_gene_expression(tr, "stur:STURON_00823") +
  plot_gene_expression(tr, "bcw:Q7M_741") +
  plot_gene_expression(tr, "lbj:LBJ_0600") +
  
  # F6P
  plot_gene_expression(tr, "lbj:LBJ_1370") +
  plot_gene_expression(tr, "lbj:LBJ_1326") +
  plot_gene_expression(tr, "lbj:LBJ_4094") +
  plot_layout(nrow = 1)

plt

# You can save this as PDF
# ggsave("JB197_glycolysis_transcripts_250422.pdf", plt, height = 3, width = 16)


plot_g <- plot_gene_expression(tr, "lbj:LBJ_2173")
plot_g
ggsave("JB197_glucose_transcripts.pdf", plot_g, height = 2, width = 2)

plot_g6p1 <- plot_gene_expression(tr, "bcw:Q7M_741")
plot_g6p2 <- plot_gene_expression(tr, "lbj:LBJ_0600")
plot_g6p1 + plot_g6p2

ggsave("JB197_g6p1_transcripts_250422.pdf", plot_g6p1, height = 2, width = 2)
ggsave("JB197_g6p2_transcripts_250422.pdf", plot_g6p2, height = 2, width = 2)

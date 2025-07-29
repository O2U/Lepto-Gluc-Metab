library(tidyverse)
library(ggpubr)

importance_df <- read_delim("data/processed/RandomForest/RandomForest_JB197_Reactions.tsv", delim = "\t")
imp <- head(importance_df, 10)
imp$Feature

rxns <- read_delim("data/raw/ModelSEED/reactions.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

rxn_imp <- rxns %>%
  filter(id %in% sub("_c$", "", imp$Feature)) %>% 
  mutate(id = factor(id, levels = sub("_c$", "", imp$Feature))) %>%
  arrange(id)
rxn_imp

rxn_imp_table <- rxn_imp %>% select(id, definition)
rxn_imp_table
colnames(rxn_imp_table) <- c("ModelSEED ID", "Definition")

# 追加する新しい行
new_row <- data.frame(`ModelSEED ID` = "cpd00027", Definition = "(1) D-Glucose[0] <=>")
colnames(new_row) <- c("ModelSEED ID", "Definition")

# 追加したい位置（たとえば2行目の後に）
position <- 2

# 行を追加
rxn_imp_table <- rbind(
  rxn_imp_table[1:position, ],
  new_row,
  rxn_imp_table[(position + 1):nrow(rxn_imp_table), ]
)

id_table <- ggtexttable(rxn_imp_table, rows = NULL)
id_table

ggsave("results/figures/Fig2/ModelSEED_ID_table.pdf", id_table, width = 9, height = 3)

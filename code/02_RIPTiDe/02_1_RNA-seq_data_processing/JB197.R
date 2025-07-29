library(tidyverse)

exp1 <- read_csv("JB197.csv")
head(exp1)

gb1 <- read_csv("../IDchanger/locus_tag_uniq_JB197.csv")
head(gb1)

brc1 <- read_csv("../IDchanger/BVBRC_genome_feature_JB197.csv")
head(brc1)

ids1 <- brc1 %>% left_join(gb1, by = c(`RefSeq Locus Tag` = "old_locus_tag"))

expt1 <- ids1 %>% left_join(exp1, by = c("locus_tag" = "GeneID"))
exp.table1 <- expt1 %>%
  select(`BRC ID`,
         JB197_29A, JB197_29B, JB197_29C, JB197_29D,
         JB197_37A, JB197_37B, JB197_37C, JB197_37D)
exp.table1.29 <- expt1 %>%
  select(`BRC ID`,
         JB197_29A, JB197_29B, JB197_29C, JB197_29D)
exp.table1.37 <- expt1 %>%
  select(`BRC ID`,
         JB197_37A, JB197_37B, JB197_37C, JB197_37D)

write_tsv(exp.table1.29, "jb197_29.tsv", na = "0", col_names = TRUE)
write_tsv(exp.table1.37, "jb197_37.tsv", na = "0", col_names = TRUE)

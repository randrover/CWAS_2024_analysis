## Preprocessing data of Velmeshev et al., Science (2023)

rm(list = ls())
library(tidyverse)
# Lineage and branch-specific genes.
df = readxl::read_excel('~/Dropbox/CWAS_paper_WD/Data/raw_annotations/Velmeshev2023_TableS2.xlsx')
#df %>%
#  group_by(lineage) %>%
#  dplyr::count() %>% View()
table(df$lineage)

mge = c('IN1', 'PV', 'PV_MP', 'SST', 'SST_RELN') # MGE & MGE-derived
cge = c('IN2', 'VIP', 'NOS', 'RELN', 'CALB2', 'CCK', 'SV2C') # CGE & CGE-derived
exn = c('L2_3', 'L4', 'L5', 'L6', 'L5_6_IT', 'SP') # Excitatory neurons (use each)
others = c('OL', 'MG', 'PER', 'END')
ast = c('AST', 'AST_FB', 'AST_PP')

# subtype
sub_df = df %>%
  filter(lineage %in% c(mge, cge, exn, others, ast))

# mge, cge, ast
df1 = data.frame(stringsAsFactors = F)
ls1 = list(MGE.dev = mge, CGE.dev = cge, AST = ast)
for(k in names(ls1)){
  tmp0 = ls1[[k]]
  tmp1 = sub_df %>%
    filter(lineage %in% tmp0) %>%
    pull(gene) %>%
    unique() %>%
    sort()
  tmp1 = as.data.frame(tmp1)
  colnames(tmp1) = 'Gene'
  tmp1 = tmp1 %>%
    mutate(Lineage = k)
  df1 = rbind(df1, tmp1)
}

# exn, others
df2 = data.frame(stringsAsFactors = F)
sub_df_others = sub_df %>%
  filter(lineage %in% c(exn, others))
for(k in unique(sub_df_others$lineage)){
  print(k)
  tmp = sub_df_others %>%
    filter(lineage==k) %>%
    pull(gene) %>%
    unique() %>%
    sort()
  tmp = as.data.frame(tmp)
  colnames(tmp) = 'Gene'
  tmp = tmp %>%
    mutate(Lineage = k)
  df2 = rbind(df2, tmp)
}

fin_df = rbind(df1, df2)

write.table(fin_df,
            '~/Dropbox/CWAS_paper_WD/Data/raw_annotations/Velmeshev2023_markers.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')


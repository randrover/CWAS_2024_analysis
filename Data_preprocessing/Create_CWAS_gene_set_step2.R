## Create a CWAS gene set: Step 2

rm(list = ls())
library(tidyverse)

gm0 = read.delim('~/Dropbox/gene-mat-recipe_v2/my_gene_matrix_231123.txt')
head(gm0)
gm0 = gm0 %>%
  dplyr::select(gene_id, gene_name, ProteinCoding,
                lincRNA = Long.ncRNA)

gm0$tmp_id = do.call(rbind.data.frame, strsplit(x = gm0$gene_id, split = '.', fixed = T))[[1]]

hgnc = read.delim('~/Dropbox/gene-mat-recipe_v2/data/hgnc_complete_set.txt')

# Werling 2018
ss = readxl::read_excel('~/Dropbox/Resources/Papers_suppl/Werling_2018/suppTables/SuppTable06_GeneLists.xlsx')
table(ss$CHD8_targets_Cotney2015_Sugathan2014)
table(ss$ASD_coexpression_networks_Willsey2013)
table(ss$BrainExpressed_Kang2011)

# An 2018
tt = readxl::read_excel('~/Dropbox/Resources/Papers_suppl/An_2018/aat6576_Table-S8.xlsx', sheet = 2)
tt$tmp_id = do.call(rbind.data.frame, strsplit(x = tt$EnsemblGeneId, split = '.', fixed = T))[[1]]
table(tt$`CHD8 targets`)
table(tt$`Midfetal co-expression genes`)
table(tt$`Brain expressed genes`)

g1 = ss %>% filter(BrainExpressed_Kang2011==1) %>% pull(Genes_background)
g2 = tt %>% filter(`Brain expressed genes`==1) %>% pull(Genes)



# DDD risk genes (Kaplanis 2020)
d1 = readxl::read_excel('~/Dropbox/Resources/Papers_suppl/Kaplanis_2020/41586_2020_2832_MOESM4_ESM.xlsx')
table(d1$significant)
d1 = d1 %>% filter(significant==TRUE)
table(d1$symbol %in% gm0$gene_name) # 3 non-match
d1$symbol[!(d1$symbol %in% gm0$gene_name)] # COL4A3BP, H3F3A, HIST1H1E
## gene change
#COL4A3BP (=CERT1)
d1$symbol = ifelse(d1$symbol=='COL4A3BP', 'CERT1', d1$symbol)
#H3F3A (=H3-3A)
d1$symbol = ifelse(d1$symbol=='H3F3A', 'H3-3A', d1$symbol)
#HIST1H1E (=H1-4)
d1$symbol = ifelse(d1$symbol=='HIST1H1E', 'H1-4', d1$symbol)

table(d1$symbol %in% gm0$gene_name) # all match

gm0$DDD285 = ifelse(gm0$gene_name %in% d1$symbol,
                    1, 0)
table(gm0$DDD285)


# CHD8 target genes (Cotney 2015, Sugathan 2014)
chd8 = tt %>%
  filter(`CHD8 targets`==1)
table(chd8$tmp_id %in% gm0$tmp_id) # all match
gm0$CHD8Common = ifelse(gm0$tmp_id %in% chd8$tmp_id, 1, 0)
table(gm0$CHD8Common)



# FMRP target genes (Darnell 2011)
FMRPs = read.delim('~/Dropbox/Resources/Papers_suppl/Darnell_2011/FMRP_gencode.v44_matched_20231123.txt')
table(FMRPs$Gene.Symbol %in% gm0$gene_name) # five non-match
FMRPs[!(FMRPs$Gene.Symbol %in% gm0$gene_name),] # 2 not matching as expected. (SCD2, CBX6-NPTXR)

gm0$FMRPDarnell = ifelse(gm0$gene_name %in% FMRPs$Gene.Symbol, 1, 0)
table(gm0$FMRPDarnell) # 838



# LOEUF<0.37
if(T){
  lof_metric = read.delim('~/Dropbox/Resources/gnomAD_resources/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
  table(lof_metric$oe_lof_upper<0.37)
  lof_metric %>% filter(oe_lof_upper<0.37) %>% pull(gene) %>% unique() %>% length()
  loeuf_g = lof_metric %>% filter(oe_lof_upper<0.37)
  
  table(loeuf_g$gene_id %in% gm0$gene_id)
  table(loeuf_g$gene_id %in% gm0$tmp_id) # 28 non-match
  # Non-matching gene ids
  #loeuf_g[!(loeuf_g$gene_id %in% gm0$tmp_id),] %>% View()
  non_match1 = loeuf_g[!(loeuf_g$gene_id %in% gm0$tmp_id),] # Keep
  
  # find duplicates
  which(duplicated(loeuf_g$gene_id))
  which(duplicated(loeuf_g$gene))
  loeuf_g[duplicated(loeuf_g$gene),] # DCAF8, IDS, SLC5A3, MDGA2
  #loeuf_g %>% filter(gene %in% loeuf_g[duplicated(loeuf_g$gene),]$gene) %>% arrange(gene) %>% View()
  # One gene has two gene ids with the same LOEUF.
  # Three genes have two gene ids each with the different LOEUF.
  #loeuf_g %>% filter(gene %in% loeuf_g[duplicated(loeuf_g$gene),]$gene) %>% arrange(gene) %>%
  #  dplyr::select(gene, gene_id, oe_lof_upper, pLI) %>%
  #  mutate(is_gencode = gene_id %in% gm0$tmp_id) %>% View()
  # Only use gene ids that match with the LOEUF gene ids. Some gene ids match with gene ids in GENCODE but gene name not matching.
  gm0 = gm0 %>%
    mutate(LOEUF37 = ifelse(tmp_id %in% loeuf_g$gene_id, 1, 0))
  table(gm0$LOEUF37) # 3166
  
  table(non_match1$gene %in% gm0$gene_name) # 21 match
  non_match1 = non_match1 %>% dplyr::select(gene, gene_id, oe_lof_upper, pLI)
  #gm0 %>% filter(gene_name %in% non_match1$gene) %>% View()
  table(non_match1$gene_id %in% gm0$tmp_id)
  # Set 1 if same gene name
  gm0$LOEUF37 = ifelse(gm0$gene_name %in% non_match1$gene,
                       1,
                       gm0$LOEUF37)
  table(gm0$LOEUF37) # 3185
  
  #gm0 %>% filter(LOEUF37==1) %>% pull(gene_name) %>% unique() %>% length()
  #gm0 %>% filter(LOEUF37==1) %>% pull(gene_id) %>% unique() %>% length()
  
  non_match2 = loeuf_g %>%
    filter(!(gene_id %in% gm0$tmp_id) & !(gene %in% gm0$gene_name)) %>%
    dplyr::select(gene, gene_id, oe_lof_upper, pLI)
  # Check with hgnc
  hgnc = read.delim('~/Dropbox/gene-mat-recipe_v2/data/hgnc_complete_set.txt')
  table(non_match2$gene %in% hgnc$prev_symbol)
  hgnc_with_loeuf = hgnc %>% filter(prev_symbol %in% non_match2$gene)
  gm0$LOEUF37 = ifelse(gm0$gene_name %in% hgnc_with_loeuf$symbol,
                       1,
                       gm0$LOEUF37)
  table(gm0$LOEUF37)
}




gm = gm0

## ASD 72 -> 185 genes (Fu et al., 2022)
asd_g = read_csv('~/Dropbox/CWAS_paper_WD/Data/raw_annotations/Fu_ASD_185_risk_genes.csv')
## gene change
asd_g$gene = ifelse(asd_g$gene=='H3F3A', 'H3-3A', asd_g$gene)
asd_g$gene = ifelse(asd_g$gene=='SUV420H1', 'KMT5B', asd_g$gene)
asd_g$gene = ifelse(asd_g$gene=='SRPR', 'SRPRA', asd_g$gene)
asd_g$gene = ifelse(asd_g$gene=='C16orf72', 'HAPSTR1', asd_g$gene)
gm$ASD185 = ifelse(gm$gene_name %in% asd_g$gene,
                   1,
                   0)
table(gm$ASD185)




## Cell type gene list (Velmeshev et al., 2023)
ct_g = read.delim('~/Dropbox/CWAS_paper_WD/Data/raw_annotations/Velmeshev2023_markers.txt')
ct_g$Lineage = gsub(x = ct_g$Lineage, pattern = '_', replacement = '.', fixed = T)

alt_g = readxl::read_excel('~/Dropbox/CWAS_paper_WD/Data/raw_annotations/Velmeshev_nonmatch_gene.xlsx')
table(is.na(alt_g$`New Gene ID`))
alt_g = alt_g %>%
  filter(!is.na(`New Gene ID`))
for(m in 1:nrow(alt_g)){
  print(m)
  ct_g$Gene = ifelse(ct_g$Gene==alt_g[m,]$`Gene ID`, alt_g[m,]$`New Gene ID`, ct_g$Gene)
}

table(unique(ct_g$Gene) %in% gm$gene_name)
ct_g %>%
  filter(!(Gene %in% gm$gene_name)) %>%
  dplyr::select(Gene) %>%
  unique() %>% nrow()

for(k in unique(ct_g$Lineage)){
  g0 = ct_g %>%
    filter(Lineage==k) %>%
    pull(Gene)
  gm[, k] <- ifelse(gm$gene_name %in% g0, 1, 0)
}


write.table(gm,
            '~/Dropbox/CWAS_paper_WD/Data/cwas_annotation.v5/gene_matrix_20231221.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')


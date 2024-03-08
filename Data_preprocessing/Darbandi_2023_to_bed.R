## Download the supplementary table for Darbandi et al., bioRxiv (2023) and convert it to bed files
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248876

library(tidyverse)

d0 = readxl::read_excel('~/Dropbox/Resources/Human_5TF_20230829/TableS3_5TRa_Peaks_May23.xlsx', sheet = 'Human_5TRa_Wide')
table(d0$prox_distal)
# promoter
t1 = d0 %>%
  filter(prox_distal == 'proximal') %>%
  dplyr::select(`#chr_hg38`, start_hg38, end_hg38)
# enhancer
t2 = d0 %>%
  filter(prox_distal == 'distal') %>%
  dplyr::select(`#chr_hg38`, start_hg38, end_hg38)

table(d0$`#chr_hg38`)

write.table(t1,
            '~/Dropbox/Resources/Human_5TF_20230829/ASD5TF_proximal.bed',
            quote = F, row.names = F, col.names = F, sep = '\t')
write.table(t2,
            '~/Dropbox/Resources/Human_5TF_20230829/ASD5TF_distal.bed',
            quote = F, row.names = F, col.names = F, sep = '\t')

# OCR
d1 = readxl::read_excel('~/Dropbox/Resources/Human_5TF_20230829/TableS1_Peaks_Cortex_May23.xlsx', sheet = 'ATAC_GW18_GW19_hg38')
d2 = d1 %>%
  dplyr::select(chromosome, start_coord, end_coord)
write.table(d2,
            '~/Dropbox/Resources/Human_5TF_20230829/ATAC_GW18_GW19_hg38.bed',
            quote = F, row.names = F, col.names = F, sep = '\t')

# TFs
TFs = c('ARID1B_GW23_hg38', 'BCL11A_GW23_hg38', 'FOXP1_GW23_hg38', 'TBR1_GW23_hg38', 'TCF7L2_GW23_hg38')
for(k in TFs){
  print(k)
  tmp = readxl::read_excel('~/Dropbox/Resources/Human_5TF_20230829/TableS1_Peaks_Cortex_May23.xlsx', sheet = k)
  tmp2 = tmp %>%
    dplyr::select(chromosome, start_coord, end_coord)
  write.table(tmp2,
              paste('~/Dropbox/Resources/Human_5TF_20230829/unfiltered.', k, '.bed', sep = ''),
              quote = F, row.names = F, col.names = F, sep = '\t')
}
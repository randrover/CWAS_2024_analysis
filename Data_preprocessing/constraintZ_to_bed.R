## Download the supplementary data of Chen et al., Nature (2023) and convert to bed file

library(tidyverse)

d0 = data.table::fread('~/Dropbox/Resources/Papers_suppl/Chen_constraint_2022/constraint_z_genome_1kb.qc.download.txt.gz')
class(d0$z)
d1 = d0 %>% filter(z >= 4)

d2 = d1 %>% dplyr::select(chrom, start, end, z)
write.table(d2,
            '~/Dropbox/Resources/Papers_suppl/Chen_constraint_2022/constraint_z_genome_1kb.qc.bed',
            quote = F, row.names = F, col.names = F, sep = '\t')
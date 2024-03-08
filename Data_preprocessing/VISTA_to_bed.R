## Visel et al., Nucleic Acids Research (2007)
## Download VISTA database from https://enhancer.lbl.gov/frnt_page_n.shtml

library(tidyverse)

d0 = read.delim('~/Dropbox/Resources/VISTA/Human.pos.brain.VISTA', header = F)
d0$POS = do.call(rbind.data.frame, strsplit(x = d0$V1, split = '|', fixed = T))[[2]]
d0$CHR = do.call(rbind.data.frame, strsplit(x = d0$POS, split = ':', fixed = T))[[1]]
d0$POS2 = do.call(rbind.data.frame, strsplit(x = d0$POS, split = ':', fixed = T))[[2]]
d0$START = do.call(rbind.data.frame, strsplit(x = d0$POS2, split = '-', fixed = T))[[1]]
d0$END = do.call(rbind.data.frame, strsplit(x = d0$POS2, split = '-', fixed = T))[[2]]

write.table(d0 %>%
              dplyr::select(CHR, START, END),
            '~/Dropbox/Resources/VISTA/Human.pos.brain.VISTA.bed',
            quote = F, row.names = F, col.names = F, sep = '\t')
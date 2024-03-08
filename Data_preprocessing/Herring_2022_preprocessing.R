## Preprocessing data of Herring et al., Cell (2022)

library(GenomicRanges)
library(data.table)
library(dplyr)
library(stringr)

atac_dir <- '~/Dropbox/CWAS_paper_WD/Data/raw_annotations'
peaks <- readRDS(paste0(atac_dir, "/Herring2022_cell_type_atac_peaks_filtered_anno_gr.rds"))

cell_types <- names(peaks)
stages <- c("Fetal")
res_dir <- paste0(atac_dir, "/Peaks")

## Peaks
for (c in cell_types){

    print(c)

    c_peaks <- as.data.frame(peaks[[c]])
    write.table(c_peaks[c("seqnames","start","end")], file = paste0(res_dir, "/", c, "_peaks_final_merge_filtered.bed"), sep = "\t", quote = F, row.names = F, col.names = F)

    for (s in stages){
        
        print(s)
        outdir <- paste0(res_dir, "/", s)

        if (!dir.exists(outdir)){
            dir.create(outdir)
        }
        
        cs_peaks <- c_peaks[unlist(as.vector(c_peaks[s])), c("seqnames","start","end")]

        print(length(rownames(cs_peaks)))
        s = str_replace(s, "\\.", "_")
        if (length(rownames(cs_peaks))>0){
            write.table(cs_peaks, file = paste0(outdir, "/", s, "_", c, "_peaks_final_filtered.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
        }
    }
}

cre_dir = paste0(atac_dir, '/CREs2/')
if (!dir.exists(cre_dir)){
  dir.create(cre_dir)
}

## CREs
for (s in stages){
    print(s)

    outdir <- paste0(cre_dir, s)

    if (!dir.exists(outdir)){
        dir.create(outdir)
    }

    for (c in cell_types){

        print(c)

        c_peaks <- as.data.frame(peaks[[c]])
        cs_peaks <- c_peaks[unlist(as.vector(c_peaks[s])), c("seqnames","start","end","CRE")]
        cs_peaks <- cs_peaks[!is.na(cs_peaks$CRE),]

        print(length(rownames(cs_peaks)))
        
        s = str_replace(s, " ", "_")
        if (length(rownames(cs_peaks))>0){
            write.table(cs_peaks, file = paste0(outdir, "/", c, "_sig_gene-peak_pairs.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
        }
    }
}

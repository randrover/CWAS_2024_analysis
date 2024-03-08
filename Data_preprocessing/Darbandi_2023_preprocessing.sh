#!/bin/sh

## Preprocessing data of Darbandi et al., bioRxiv (2023)
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248876

Rscript Darbandi_2023_to_bed.R

## Common promoter & enhancer of five ASD TF (TableS3_5TRa_Peaks_May23.xlsx)
# merge bed files
cat ASD5TF_proximal.bed ASD5TF_distal.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | mergeBed -i stdin > ASD5TF_common.bed
# sort bed files
cat ASD5TF_common.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.ASD5TF_common.bed.gz
tabix -p bed sorted.ASD5TF_common.bed.gz

## Binding sites of five ASD TF (TableS1_Peaks_Cortex_May23.xlsx)
# Presort bed files
cat ATAC_GW18_GW19_hg38.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > ATAC_GW18_GW19_hg38.bed.gz
tabix -p bed ATAC_GW18_GW19_hg38.bed.gz

cat unfiltered.ARID1B_GW23_hg38.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > unfiltered.ARID1B_GW23_hg38.bed.gz
tabix -p bed unfiltered.ARID1B_GW23_hg38.bed.gz

cat unfiltered.BCL11A_GW23_hg38.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > unfiltered.BCL11A_GW23_hg38.bed.gz
tabix -p bed unfiltered.BCL11A_GW23_hg38.bed.gz

cat unfiltered.FOXP1_GW23_hg38.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > unfiltered.FOXP1_GW23_hg38.bed.gz
tabix -p bed unfiltered.FOXP1_GW23_hg38.bed.gz

cat unfiltered.TBR1_GW23_hg38.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > unfiltered.TBR1_GW23_hg38.bed.gz
tabix -p bed unfiltered.TBR1_GW23_hg38.bed.gz

cat unfiltered.TCF7L2_GW23_hg38.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > unfiltered.TCF7L2_GW23_hg38.bed.gz
tabix -p bed unfiltered.TCF7L2_GW23_hg38.bed.gz

# Find intersected regions using bedtools intersect
bedtools intersect -a unfiltered.ARID1B_GW23_hg38.bed.gz -b ATAC_GW18_GW19_hg38.bed.gz | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > filtered.ARID1B_GW23_hg38.bed.gz
bedtools intersect -a unfiltered.TBR1_GW23_hg38.bed.gz -b ATAC_GW18_GW19_hg38.bed.gz | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > filtered.TBR1_GW23_hg38.bed.gz
bedtools intersect -a unfiltered.BCL11A_GW23_hg38.bed.gz -b ATAC_GW18_GW19_hg38.bed.gz | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > filtered.BCL11A_GW23_hg38.bed.gz
bedtools intersect -a unfiltered.TCF7L2_GW23_hg38.bed.gz -b ATAC_GW18_GW19_hg38.bed.gz | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > filtered.TCF7L2_GW23_hg38.bed.gz
bedtools intersect -a unfiltered.FOXP1_GW23_hg38.bed.gz -b ATAC_GW18_GW19_hg38.bed.gz | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > filtered.FOXP1_GW23_hg38.bed.gz

tabix -p bed filtered.ARID1B_GW23_hg38.bed.gz
tabix -p bed filtered.BCL11A_GW23_hg38.bed.gz
tabix -p bed filtered.FOXP1_GW23_hg38.bed.gz
tabix -p bed filtered.TBR1_GW23_hg38.bed.gz
tabix -p bed filtered.TCF7L2_GW23_hg38.bed.gz
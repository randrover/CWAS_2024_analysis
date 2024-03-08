#!/bin/sh

## Download JARVIS data (version 1.1) from Vitsios et al., Nature Communications (2021)
## Preprocessing

gzcat *hg38.tsv.gz | awk 'BEGIN{OFS="\t"} {print "chr"$1, $2-1, $2, $3}' | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.jarvis.all_chr.bed.gz
tabix -p bed sorted.jarvis.all_chr.bed.gz

gzcat sorted.jarvis.all_chr.bed.gz | awk 'BEGIN{OFS="\t"} ($4>=0.99)' | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.jarvis.all_chr.099_.bed.gz
tabix -p bed sorted.jarvis.all_chr.099_.bed.gz
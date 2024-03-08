#!/bin/sh

## Conservation score from the following studies
## 1. Kuderna et al., Nature (2023)
## 2. Siepel et al., Genome Research (2005)

# BigWig → Wig (https://www.biostars.org/p/71692/)
cd lib
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
# Executable
chmod +x bigWigToWig

# bigWigToWig.pl tmp.bw > tmp.wig
bigWigToWig.pl vertebrate.phastCons46way.hg19ToHg38.bw > vertebrate.phastCons46way.hg19ToHg38.wig
bigWigToWig.pl vertebrate.phyloP46way.hg19ToHg38.bw > vertebrate.phyloP46way.hg19ToHg38.wig

# Wig → BED
wig2bed --do-not-sort < vertebrate.phastCons46way.hg19ToHg38.wig > vertebrate.phastCons46way.hg19ToHg38.bed
wig2bed --do-not-sort < vertebrate.phyloP46way.hg19ToHg38.wig > vertebrate.phyloP46way.hg19ToHg38.bed

# Filter by cutoff
cat vertebrate.phastCons46way.hg19ToHg38.bed | awk -F"\t" 'BEGIN{OFS="\t";} {if ($5>=0.2) print $1, $2, $3, $4, $5}' | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > vertebrate.phastCons46way.hg19ToHg38.over02.bed.gz
tabix -p bed vertebrate.phastCons46way.hg19ToHg38.over02.bed.gz

cat vertebrate.phyloP46way.hg19ToHg38.bed | awk -F"\t" 'BEGIN{OFS="\t";} {if ($5>=2) print $1, $2, $3, $4, $5}' | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > vertebrate.phyloP46way.hg19ToHg38.over2.bed.gz
tabix -p bed vertebrate.phyloP46way.hg19ToHg38.over2.bed.gz
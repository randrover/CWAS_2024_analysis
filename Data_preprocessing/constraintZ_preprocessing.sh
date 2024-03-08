#!/bin/sh

## Preprocessing data of Chen et al., Nature (2023)

Rscript constraintZ_to_bed.R

cat constraint_z_genome_1kb.qc.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.constraint_z_genome_1kb.qc.bed.gz
tabix -p bed sorted.constraint_z_genome_1kb.qc.bed.gz
#!/bin/sh

## Preprocessing data of Visel et al., Nucleic Acids Research (2007)

# Go to link (https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page=1;show=1;form=search;action=search;page_size=100;search.form=no;search.result=yes;search.sequence=1)
# Save data to text file

# Filter only human enhancers.
cat VISTA| grep Human > Human.VISTA

# Filter only positive enhancers.
cat Human.VISTA | grep positive > Human.pos.VISTA

# Filter only brain regions.
cat Human.pos.VISTA | grep -E 'brain|neural' > Human.pos.brain.VISTA

# Save as bed file.
Rscript VISTA_to_bed.R

# Sort file.
cat Human.pos.brain.VISTA.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.Human.pos.brain.VISTA.bed.gz
tabix -p bed sorted.Human.pos.brain.VISTA.bed.gz

# liftover to hg38.
# Download liftover tool
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/liftOver
chmod +x liftOver # Add execute mode
# or use conda liftover

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
resource_path=/Users/yujinkim/Dropbox/Resources/
chain=$resource_path"hg19ToHg38.over.chain.gz"

gunzip -k sorted.Human.pos.brain.VISTA.bed.gz

~/lib/liftOver sorted.Human.pos.brain.VISTA.bed $chain sorted.Human.pos.brain.VISTA.hg38.lifted.bed sorted.Human.pos.brain.VISTA.hg38.unlifted.bed

# sort file again.
cat sorted.Human.pos.brain.VISTA.hg38.lifted.bed | sort -k1,1V -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.Human.pos.brain.VISTA.hg38.lifted.bed.gz
tabix -p bed sorted.Human.pos.brain.VISTA.hg38.lifted.bed.gz

rm sorted.Human.pos.brain.VISTA.hg38.lifted.bed sorted.Human.pos.brain.VISTA.bed
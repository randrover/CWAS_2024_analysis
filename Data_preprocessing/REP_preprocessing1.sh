#!/bin/bash

## Preprocessing data of Roadmap Epigenomics Project
## Among these, we only used fetal stage samples

chain=/data1/Resources/liftOver_chain/hg19ToHg38.over.chain.gz

echo $(date)

for m in $(seq 1 1 9)
do
    echo "liftover E00${m}"
    liftOver regions_prom_E00${m}.bed.gz $chain regions_prom_E00${m}.hg38.lifted.bed regions_prom_E00${m}.hg38.unlifted.bed

    echo "Sort E00${m}"

    cat regions_prom_E00${m}.hg38.lifted.bed | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.regions_prom_E00${m}.hg38.lifted.bed.gz
    tabix -p bed sorted.regions_prom_E00${m}.hg38.lifted.bed.gz
    
    # Check the exit status of wget
    if [ $? -ne 0 ]; then
        echo "Failed to liftover: E00${m}"
        continue
    fi
done

for m in $(seq 10 1 99)
do
    echo "liftover E0${m}"
    liftOver regions_prom_E0${m}.bed.gz $chain regions_prom_E0${m}.hg38.lifted.bed regions_prom_E0${m}.hg38.unlifted.bed

    echo "Sort E0${m}"

    cat regions_prom_E0${m}.hg38.lifted.bed | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.regions_prom_E0${m}.hg38.lifted.bed.gz
    tabix -p bed sorted.regions_prom_E0${m}.hg38.lifted.bed.gz

    # Check the exit status of wget
    if [ $? -ne 0 ]; then
        echo "Failed to liftover: E0${m}"
        continue
    fi
done

for m in $(seq 100 1 129)
do
    echo "liftover E${m}"
    liftOver regions_prom_E${m}.bed.gz $chain regions_prom_E${m}.hg38.lifted.bed regions_prom_E${m}.hg38.unlifted.bed

    echo "Sort E${m}"

    cat regions_prom_E${m}.hg38.lifted.bed | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.regions_prom_E${m}.hg38.lifted.bed.gz
    tabix -p bed sorted.regions_prom_E${m}.hg38.lifted.bed.gz

    # Check the exit status of wget
    if [ $? -ne 0 ]; then
        echo "Failed to liftover: E${m}"
        continue
    fi
done


for m in $(seq 1 1 9)
do
    echo "liftover E00${m}"
    liftOver regions_enh_E00${m}.bed.gz $chain regions_enh_E00${m}.hg38.lifted.bed regions_enh_E00${m}.hg38.unlifted.bed

    echo "Sort E00${m}"

    cat regions_enh_E00${m}.hg38.lifted.bed | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.regions_enh_E00${m}.hg38.lifted.bed.gz
    tabix -p bed sorted.regions_enh_E00${m}.hg38.lifted.bed.gz

    # Check the exit status of wget
    if [ $? -ne 0 ]; then
        echo "Failed to liftover: E00${m}"
        continue
    fi
done

for m in $(seq 10 1 99)
do
    echo "liftover E0${m}"
    liftOver regions_enh_E0${m}.bed.gz $chain regions_enh_E0${m}.hg38.lifted.bed regions_enh_E0${m}.hg38.unlifted.bed

    echo "Sort E0${m}"

    cat regions_enh_E0${m}.hg38.lifted.bed | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.regions_enh_E0${m}.hg38.lifted.bed.gz
    tabix -p bed sorted.regions_enh_E0${m}.hg38.lifted.bed.gz

    # Check the exit status of wget
    if [ $? -ne 0 ]; then
        echo "Failed to liftover: E0${m}"
        continue
    fi
done

for m in $(seq 100 1 129)
do
    echo "liftover E${m}"
    liftOver regions_enh_E${m}.bed.gz $chain regions_enh_E${m}.hg38.lifted.bed regions_enh_E${m}.hg38.unlifted.bed

    echo "Sort E${m}"

    cat regions_enh_E${m}.hg38.lifted.bed | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > sorted.regions_enh_E${m}.hg38.lifted.bed.gz
    tabix -p bed sorted.regions_enh_E${m}.hg38.lifted.bed.gz

    # Check the exit status of wget
    if [ $? -ne 0 ]; then
        echo "Failed to liftover: E${m}"
        continue
    fi
done

echo $(date)
echo "Done"
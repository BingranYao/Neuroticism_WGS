#!/bin/bash

pheno=$1
chr=$2
type=$3

# Quantitative Trait
regenie \
    --step 2 \
    --chr ${chr} \
    --keep /home/dnanexus/data/Neu_id.txt \
    --bed /home/dnanexus/Q0_unre_Caucasian_EUR_c${chr} \
    --phenoFile /home/dnanexus/data/Neuroticism_qt.txt \
    --phenoCol ${pheno} \
    --covarFile /home/dnanexus/data/cov_pc20.txt \
    --catCovarList sex \
    --pred /home/dnanexus/step1/qt_pred.list \
    --anno-file /home/dnanexus/Anno_New/chr${chr}/Main/${type}_chr${chr}.txt \
    --set-list  /home/dnanexus/Anno_New/chr${chr}/Main/chr${chr}_${type}.setlist \
    --mask-def /home/dnanexus/Anno_New/Mask/Mask_${type}.txt \
    --aaf-bins 0.01 \
    --vc-tests skato,acato-full \
    --vc-maxAAF 0.01 \
    --minMAC 1 \
    --qt \
    --rgc-gene-p \
    --joint minp,acat,sbat \
    --bsize 200 \
    --threads 32 \
    --write-mask \
    --check-burden-files \
    --write-samples \
    --write-mask-snplist \
    --print-pheno \
    --condition-list /home/dnanexus/Cond_clump/${pheno}/snplist/chr${chr}_snplist.txt
    --out /home/dnanexus/step2/qt/chr${chr}_${type}

# Binary Trait
regenie \
    --step 2 \
    --chr ${chr} \
    --keep /home/dnanexus/data/Neu_id.txt \
    --bed /home/dnanexus/Q0_unre_Caucasian_EUR_c${chr} \
    --phenoFile /home/dnanexus/data/Neuroticism_bt.txt \
    --phenoCol ${pheno} \
    --covarFile /home/dnanexus/data/cov_pc20.txt \
    --catCovarList sex \
    --pred /home/dnanexus/step1/bt_pred.list \
    --anno-file /home/dnanexus/Anno_New/chr${chr}/Main/${type}_chr${chr}.txt \
    --set-list  /home/dnanexus/Anno_New/chr${chr}/Main/chr${chr}_${type}.setlist \
    --mask-def /home/dnanexus/Anno_New/Mask/Mask_${type}.txt \
    --aaf-bins 0.01 \
    --vc-tests skato,acato-full \
    --vc-maxAAF 0.01 \
    --minMAC 1 \
    --bt \
    --firth --approx \
    --firth-se \
    --pThresh 0.999 \
    --rgc-gene-p \
    --joint minp,acat \
    --bsize 200 \
    --threads 32 \
    --write-mask \
    --check-burden-files \
    --write-samples \
    --write-mask-snplist \
    --print-pheno \
    --condition-list /home/dnanexus/Cond_clump/${pheno}/snplist/chr${chr}_snplist.txt
    --out /home/dnanexus/step2/bt/chr${chr}_${type}

  
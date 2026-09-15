#!/bin/bash

# Quantitative Trait
regenie \
    --step 2 \
    --chr ${chr} \
    --keep /home/dnanexus/data/Neu_id.txt \
    --bed /home/dnanexus/Q0_unre_Caucasian_EUR_c${chr} \
    --phenoFile /home/dnanexus/data/qt.txt \
    --phenoCol ${pheno} \
    --covarFile /home/dnanexus/data/cov_pc20.txt \
    --catCovarList sex \
    --pred /home/dnanexus/step1/qt_pred.list \
    --anno-file /home/dnanexus/lovofiles/${base_mask}_chr${chr}.txt \
    --set-list /home/dnanexus/lovofiles/chr${chr}_${base_mask}.setlist \
    --mask-def /home/dnanexus/Anno_New/Mask/Mask_${base_mask}.txt \
    --aaf-bins 0.01 \
    --vc-tests skato,acato-full \
    --vc-maxAAF 0.01 \
    --minMAC 1 \
    --qt \
    --bsize 200 \
    --threads 8 \
    --mask-lovo ${lovo_spec} \
    --write-samples \
    --print-pheno \
    --check-burden-files \
    --out /home/dnanexus/step2/qt/chr${chr}_${base_mask}

# Binary Trait
regenie \
    --step 2 \
    --chr ${chr} \
    --keep /home/dnanexus/data/Neu_id.txt \
    --bed /home/dnanexus/Q0_unre_Caucasian_EUR_c${chr} \
    --phenoFile /home/dnanexus/data/bt.txt \
    --phenoCol ${pheno} \
    --covarFile /home/dnanexus/data/cov_pc20.txt \
    --catCovarList sex \
    --pred /home/dnanexus/step1/bt_pred.list \
    --anno-file /home/dnanexus/lovofiles/${base_mask}_chr${chr}.txt \
    --set-list /home/dnanexus/lovofiles/chr${chr}_${base_mask}.setlist \
    --mask-def /home/dnanexus/Anno_New/Mask/Mask_${base_mask}.txt \
    --aaf-bins 0.01 \
    --vc-tests skato,acato-full \
    --vc-maxAAF 0.01 \
    --minMAC 1 \
    --bt \
    --firth --approx \
    --firth-se \
    --pThresh 0.999 \
    --bsize 200 \
    --threads 8 \
    --mask-lovo ${lovo_spec} \
    --write-samples \
    --print-pheno \
    --check-burden-files \
    --out /home/dnanexus/step2/bt/chr${chr}_${base_mask}

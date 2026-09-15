#!/bin/bash

chr=$1
# Quantitative Trait
regenie \
--step 2 \
--chr ${chr} \
--bed /home/dnanexus/Q0_unre_Caucasian_EUR_c${chr} \
--keep /home/dnanexus/data/Neu_id.txt \
--phenoFile /home/dnanexus/data/Neuroticism_qt.txt \
--covarFile /home/dnanexus/data/cov_pc20.txt \
--catCovarList sex \
--pred /home/dnanexus/step1/qt_pred.list \
--qt \
--minMAC 20 \
--bsize 1000 \
--threads 16 \
--write-samples \
--print-pheno \
--out /home/dnanexus/step2/qt/step2_chr${chr}

# Binary Trait
regenie \
--step 2 \
--chr ${chr} \
--bed /home/dnanexus/Q0_unre_Caucasian_EUR_c${chr} \
--keep /home/dnanexus/data/Neu_id.txt \
--phenoFile /home/dnanexus/data/Neuroticism_bt.txt \
--covarFile /home/dnanexus/data/cov_pc20.txt \
--catCovarList sex \
--pred /home/dnanexus/step1/bt_pred.list \
--bt \
--firth \
--approx \
--firth-se \
--pThresh 0.999 \
--minMAC 20 \
--bsize 1000 \
--threads 16 \
--write-samples \
--print-pheno \
--out /home/dnanexus/step2/bt/step2_chr${chr}


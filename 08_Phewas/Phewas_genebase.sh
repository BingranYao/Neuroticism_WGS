#!/bin/bash
p=parameter1
chr=parameter2
pheno_inc=parameter3
mask=parameter4
ensg=parameter5

regenie \
  --step 2 \
  --chr ${chr} \
  --bed /Q0_unre_Caucasian_c${chr} \
  --extract-setlist ${ensg} \
  --phenoFile /data/CaucasianSample.FullDisease_p${p}.txt \
  --phenoColList ${pheno_inc} \
  --covarFile /data/wgs_cov.txt \
  --catCovarList sex_0fe,batch_32053 --maxCatLevels 30 \
  --anno-file /Anno_New/chr${chr}/Main/${mask}_chr${chr}.txt \
  --set-list /Anno_New/chr${chr}/Main/chr${chr}_${mask}.setlist \
  --mask-def /Anno_New/Mask/Mask_${mask}.txt \
  --aaf-bins 0.01 \
  --vc-tests skato,acato-full \
  --vc-maxAAF 0.01 \
  --joint minp,acat \
  --rgc-gene-p \
  --check-burden-files \
  --write-mask \
  --pred /step1/all${p}/step1_all${p}_l1_pred.list \
  --bt \
  --bsize 200 \
  --threads 31 \
  --firth --firth-se --approx \
  --pThresh 0.001 \
  --write-samples \
  --write-mask-snplist \
  --print-pheno \
  --out d/result/all${p}_${gene_symbol}/${mask}/${gene_symbol}_all${p}_chr${chr}_${mask}
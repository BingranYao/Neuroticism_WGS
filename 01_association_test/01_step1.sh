#!/bin/bash

# Quantitative Trait
regenie \
--step 1 \
--bed /home/dnanexus/ukb_cal_allChrs_hg38 \
--extract /home/dnanexus/data/qc_pass_EUR.snplist \
--keep /home/dnanexus/data/Neu_id.txt \
--phenoFile /home/dnanexus/data/Neuroticism_qt.txt \
--covarFile /home/dnanexus/data/cov_pc20.txt \
--catCovarList sex \
--qt \
--bsize 1000 \
--threads 100 \
--lowmem \
--lowmem-prefix /home/dnanexus/step1/qt/tmp \
--out /home/dnanexus/step1/qt

# Binary Trait
#!/bin/bash
regenie \
--step 1 \
--bed /home/dnanexus/ukb_cal_allChrs_hg38 \
--extract /home/dnanexus/data/qc_pass_EUR.snplist \
--keep /home/dnanexus/data/Neu_id.txt \
--phenoFile /home/dnanexus/data/Neuroticism_bt.txt \
--covarFile /home/dnanexus/data/cov_pc20.txt \
--catCovarList sex \
--bt \
--bsize 1000 \
--threads 100 \
--lowmem \
--lowmem-prefix /home/dnanexus/step1/bt/tmp \
--out /home/dnanexus/step1/bt

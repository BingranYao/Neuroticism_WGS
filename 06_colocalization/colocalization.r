#!/usr/bin/env Rscript
options(warn=1, showWarnCalls=TRUE)

library(data.table)
library(dplyr)
library(coloc)

args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1]
p1 <- 1e-4; p2 <- 1e-4; p12 <- 1e-5

for(disease_file in disease_files){
  disease_name <- sub("^finngen_R12_", "", disease_file)
  Nsum <- R12_manifest[phenocode==disease_name]$num_sum

  outcome_gwas <- fread(file.path(disease_dir, disease_file))
  colnames(outcome_gwas)[1] <- "chr"
  outcome_gwas[, ID2 := paste(chr,pos,ref,alt,sep=":")]

  coloc_results <- list()

  for(snv in leadsnp$ID2){
    lead_info <- leadsnp[ID2==snv]
    lead_chr <- lead_info$chr
    start <- lead_info$start
    end <- lead_info$end

    pheno_sub <- pheno_gwas %>% filter(chr==lead_chr, pos>=start, pos<=end)
    outcome_sub <- outcome_gwas %>% filter(chr==lead_chr, pos>=start, pos<=end)
    common_snps <- intersect(pheno_sub$ID2, outcome_sub$ID2)
    pheno_snps$P <- 10^(-pheno_snps$LOG10P)

    coloc_res <- tryCatch({
      coloc.abf(
        dataset1=list(
          snp=pheno_snps$ID2,
          beta=pheno_snps$BETA,
          varbeta=pheno_snps$SE^2,
          p=pheno_snps$P,
          MAF=pheno_snps$A1FREQ,
          N=338786,
          type="cc"
        ),
        dataset2=list(
          snp=outcome_snps$ID2,
          beta=outcome_snps$beta,
          varbeta=outcome_snps$sebeta^2,
          MAF=outcome_snps$af_alt,
          N=Nsum,
          type="cc"
        ),
        p1=p1, p2=p2, p12=p12
      )
    }, error=function(e){
      message("WARNING: coloc error for ", snv, ": ", e$message)
      return(NULL)
    })

    if(!is.null(coloc_res)) coloc_results[[snv]] <- coloc_res
  }
}
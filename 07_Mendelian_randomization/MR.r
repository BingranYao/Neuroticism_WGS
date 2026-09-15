suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(readr))

R12_manifest <- fread(manifest_path)
target_pairs <- read.csv(pairs_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
target_pairs$exposure <-"Neuroticism_score_rint"
target_pairs <- target_pairs[!is.na(target_pairs$exposure) & !is.na(target_pairs$outcome), ]
target_pairs <- unique(target_pairs)

unique_outcomes <- unique(target_pairs$outcome)
total_unique_outcomes <- length(unique_outcomes)

chunk_indices <- split(seq_along(unique_outcomes), ceiling(seq_along(unique_outcomes) / (length(unique_outcomes) / total_chunks)))
my_outcome_indices <- chunk_indices[[chunk_id]]
my_outcomes <- unique_outcomes[my_outcome_indices]

for (j in seq_along(my_outcomes)) {
  d <- my_outcomes[j]
  
  # Find all exposures needed for this specific outcome
  exposures_for_d <- target_pairs$exposure[target_pairs$outcome == d]
  
  tryCatch({
    output_folder <- file.path(output_folder_base, d)
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Check if all exposure results for this outcome already exist
    all_exist <- TRUE
    for (m in exposures_for_d) {
      if (!dir.exists(file.path(output_folder, m))) {
        all_exist <- FALSE
        break
      }
    }
    if (all_exist) {
      cat(sprintf("  All exposures for outcome %s already completed. Skipping GWAS load.\n", d))
      next
    }
    
    # Load Outcome Data
    cat("  Loading FinnGen GWAS outcome data...\n")
    N <- R12_manifest[R12_manifest$phenocode == as.character(d),]$num_sum
    if(length(N) == 0) N <- NA # Handle cases where sample size is missing
    
    file_path_d = paste0(finngen_dir, "finngen_R12_", d)
    if(!file.exists(file_path_d)) {
      cat(sprintf("  ERROR: Outcome file not found: %s\n", file_path_d))
      next
    }
    
    outcome_raw <- fread(file_path_d)
    setnames(outcome_raw, "#chrom", "chr")
    outcome_raw[, SNP := paste(chr, pos, ref, alt, sep = ":")]
    outcome_raw <- as.data.frame(outcome_raw)
    
    outcome_data <- format_data(
      outcome_raw,
      type = "outcome",
      snp_col = "SNP",
      effect_allele_col = "alt", 
      other_allele_col = "ref",
      pval_col = "pval",
      chr_col = "chr",
      beta_col = "beta", 
      se_col = "sebeta",
      eaf_col = "af_alt"
    )
    outcome_data$outcome <- d
    
    # Process each exposure for this outcome
    for (i in seq_along(exposures_for_d)) {
      m <- exposures_for_d[i]
      cat(sprintf("  -> Exposure %d/%d: %s\n", i, length(exposures_for_d), m))
      
      full_output_f = file.path(output_folder, m)
      leadsnp_file = file.path(leadSNP_folder, paste0(m, "_exposure.txt"))
      
      dir.create(full_output_f, recursive = TRUE)
      
      exposure_data <- read_exposure_data(
        filename = leadsnp_file,
        sep = "\t",
        snp_col = "ID2",        
        beta_col = "BETA", 
        se_col = "SE", 
        effect_allele_col = "ALLELE1", 
        other_allele_col = "ALLELE0", 
        pval_col = "P",
        eaf_col = "A1FREQ"
      ) 
      exposure_data$exposure <- m
      write.csv(exposure_data, file = file.path(full_output_f, "exposure_data.csv"), row.names = FALSE)
      
      dat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data) 
      
      if (nrow(dat) == 0) {
        cat("     0 harmonised SNPs. Skipping MR.\n")
        next
      } else if (nrow(dat) == 1) {
        # Wald ratio for single SNP
        results <- mr(dat, method_list=c("mr_wald_ratio"))
        if (nrow(results) > 0) {
          results.withOR <- generate_odds_ratios(results)
          write.csv(results.withOR, file=file.path(full_output_f, 'wald_results.csv'), row.names=FALSE)
        }
      } else if (nrow(dat) > 1) {
        # Standard MR for multiple SNPs
        results <- mr(dat)
        results.withOR <- generate_odds_ratios(results)
        
        dat$samplesize.outcome <- N  
        dat$samplesize.exposure <- 338786
        
        dat.mr_heterogeneity <- mr_heterogeneity(dat)   
        dat.mr_pleiotropy_test <- mr_pleiotropy_test(dat) 
        leaveoneout <- mr_leaveoneout(dat)
        leaveoneout.withOR <- generate_odds_ratios(leaveoneout)
        single_snp_analysis <- mr_singlesnp(dat) 
        results.single_snp_analysis.withOR <- generate_odds_ratios(single_snp_analysis)
        
        write.csv(dat.mr_heterogeneity, file=file.path(full_output_f,"dat.mr_heterogeneity.csv"), row.names = FALSE)
        write.csv(dat.mr_pleiotropy_test, file=file.path(full_output_f,"dat.mr_pleiotropy_test.csv"), row.names = FALSE)
        write.csv(leaveoneout.withOR, file=file.path(full_output_f,"leaveoneout.withOR.csv"), row.names = FALSE)
        write.csv(results.single_snp_analysis.withOR, file=file.path(full_output_f,"single_snp_analysis.withOR.csv"), row.names = FALSE)
        
        # Format comprehensive results
        if(nrow(dat.mr_pleiotropy_test) > 0) {
          results.withOR$egger_intercept[results.withOR$method == 'MR Egger'] = dat.mr_pleiotropy_test$egger_intercept
          results.withOR$egger_intercept.pval[results.withOR$method == 'MR Egger'] = dat.mr_pleiotropy_test$pval
        }
        if(nrow(dat.mr_heterogeneity) > 0) {
          results.withOR$MR_Egger.Q[results.withOR$method == 'MR Egger'] = dat.mr_heterogeneity$Q[dat.mr_heterogeneity$method == 'MR Egger']
          results.withOR$MR_Egger.Q_df[results.withOR$method == 'MR Egger'] = dat.mr_heterogeneity$Q_df[dat.mr_heterogeneity$method == 'MR Egger']
          results.withOR$MR_Egger.Q_pval[results.withOR$method == 'MR Egger'] = dat.mr_heterogeneity$Q_pval[dat.mr_heterogeneity$method == 'MR Egger']
          results.withOR$Inverse_variance_weighted.Q[results.withOR$method == 'Inverse variance weighted'] = dat.mr_heterogeneity$Q[dat.mr_heterogeneity$method == 'Inverse variance weighted']
          results.withOR$Inverse_variance_weighted.Q_df[results.withOR$method == 'Inverse variance weighted'] = dat.mr_heterogeneity$Q_df[dat.mr_heterogeneity$method == 'Inverse variance weighted']
          results.withOR$Inverse_variance_weighted.Q_pval[results.withOR$method == 'Inverse variance weighted'] = dat.mr_heterogeneity$Q_pval[dat.mr_heterogeneity$method == 'Inverse variance weighted']
        }
        
        write.csv(results.withOR, file=file.path(full_output_f,"results.csv"), row.names = FALSE)
        write.csv(dat, file=file.path(full_output_f,"dat.csv"), row.names = FALSE)
      }
      cat("MR completed successfully.\n")
    }
    
    # 6. Organize and combine MR results for this outcome
    cat("Combining results for outcome...\n")
    setwd(output_folder)
    folders <- list.dirs(full.names = FALSE, recursive = FALSE)
    csv_names <- c("results.csv", "dat.mr_heterogeneity.csv", "dat.mr_pleiotropy_test.csv")
    wald_csv_name <- "wald_results.csv"
    all_combined <- list()
    
    for (folder in folders) {
      folder_path <- file.path(output_folder, folder)
      file_exist <- all(file.exists(file.path(folder_path, csv_names)))
      wald_exist <- file.exists(file.path(folder_path, wald_csv_name))
      
      combined <- data.frame(nsnp=NA, Egger_p=NA, Weightedmedian_p=NA, IVW_p=NA, Simplemode_p=NA, Weightedmode_p=NA, Egger_b=NA, Weightedmedian_b=NA, IVW_b=NA, Simplemode_b=NA, Weightedmode_b=NA, hetero_Egger_p=NA, hetero_IVW_p=NA, pleio_p=NA, wald_b=NA, wald_p=NA)
      
      if (file_exist) {
        df_hetero <- tryCatch({ read.csv(file.path(folder_path, "dat.mr_heterogeneity.csv")) }, error = function(e) NULL)
        if (!is.null(df_hetero) && nrow(df_hetero) > 0 && "Q_pval" %in% names(df_hetero)) {
          df_hetero_temp <- as.data.frame(t(data.frame(df_hetero$Q_pval)))
          if (ncol(df_hetero_temp) == 2) {
            colnames(df_hetero_temp) <- c("hetero_Egger_p", "hetero_IVW_p")
            combined[c("hetero_Egger_p", "hetero_IVW_p")] <- df_hetero_temp
          } else if (ncol(df_hetero_temp) == 1) {
            colnames(df_hetero_temp) <- "hetero_IVW_p"
            combined["hetero_IVW_p"] <- df_hetero_temp
          }
        }
        
        df_pleio <- tryCatch({ read.csv(file.path(folder_path, "dat.mr_pleiotropy_test.csv")) }, error = function(e) NULL)
        if (!is.null(df_pleio) && nrow(df_pleio) > 0 && "pval" %in% names(df_pleio)) {
          combined["pleio_p"] <- df_pleio$pval[1]
        }
        
        df_results <- read.csv(file.path(folder_path, "results.csv"))
        if (!is.null(df_results) && all(c("nsnp", "pval", "b") %in% names(df_results))) {
          if (nrow(df_results) > 1) {
            combined[c("nsnp", "Egger_p","Weightedmedian_p", "IVW_p", "Simplemode_p","Weightedmode_p","Egger_b", "Weightedmedian_b", "IVW_b", "Simplemode_b","Weightedmode_b")] <- 
              c(df_results[1, "nsnp"], df_results$pval[1], df_results$pval[2], df_results$pval[3], df_results$pval[4], df_results$pval[5], df_results$b[1], df_results$b[2], df_results$b[3], df_results$b[4], df_results$b[5])
          } else if (nrow(df_results) == 1) {
            combined[c("nsnp", "IVW_p", "IVW_b")] <- c(df_results[1, "nsnp"], df_results$pval[1], df_results$b[1])
          }
        }
      }
      
      if (wald_exist) {
        df_wald <- read.csv(file.path(folder_path, "wald_results.csv"))
        if (!is.null(df_wald) && nrow(df_wald) > 0 && all(c("b", "pval") %in% names(df_wald))) {
          combined[c("nsnp", "wald_b", "wald_p")] <- c(df_wald$nsnp[1], df_wald$b[1], df_wald$pval[1])
        }
      }
      all_combined[[folder]] <- combined
    }
    
    if(length(all_combined) > 0) {
      final_combined <- do.call(rbind, all_combined)
      final_combined$ID <- rownames(final_combined)  
      final_combined <- final_combined[c("ID", names(final_combined)[names(final_combined) != "ID"])]
      combined_csv_path <- file.path(output_folder, "MR_results_combined.csv")
      write.csv(final_combined, combined_csv_path, row.names = FALSE)
    }

  }, error = function(e) {
    cat(sprintf("  [ERROR] Failed processing disease %s: %s\n", d, e$message))
  })
}
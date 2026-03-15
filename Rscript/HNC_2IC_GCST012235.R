
# devtools::install_github("mrcieu/ieugwasr",force = TRUE)
# install.packages("plinkbinr-master.zip")
# devtools::install_github("explodecomputer/plinkbinr")

library(future)
library(furrr)
library(ieugwasr)
library(plinkbinr)
library(TwoSampleMR)
library(data.table)
library(MRPRESSO)
library(tidyverse)
library(R.utils)
library(readr)
library(stringr)
library(dplyr)

plan(multisession, workers = 6)

clump_p_values <- c(1e-6, 5e-6, 5e-7, 5e-8)  # 1e-6, 5e-6, 5e-7,
## these 2 line for test
# clump_p_value <- clump_p_values(4)
# file <- fileslist[1]
# Define a function that processes each exposure file
for (clump_p_value in clump_p_values){
  process_files <- function(file) {
    Index_IC <- fread("[ImmuneCells.csv]", data.table = FALSE)
    foldername <- Index_IC[Index_IC[, 5] == substring(file, 1, 12), 4]
    # fileslist <- fileslist[!grepl("\\.gz", fileslist)]
    # fileslist <- fileslist[!grepl("\\.csv", fileslist)]
    # fileslist <- fileslist[!grepl("\\.xslx", fileslist)]
    add_exposure <- "HNC/GCST012235.csv"
    add_outcome <- paste("ImmuneCells/",file,sep = "")
    pre_result_path <- paste('Results/HNC2IC_',as.character(clump_p_value),sep = "")
    dir.create(pre_result_path)
    result_path <- paste('Results/HNC2IC_',as.character(clump_p_value),'/HNC_TO_',foldername,sep = "")
    dir.create(result_path)
    setwd(result_path)
    
    # #如果存在则跳过！！！！测试！！！！
    # if (file.exists(paste(result_path,"/MR_dat.csv",sep = ""))) {
    #   message(paste("Skipping processing for", result_path, "as it already exists"))
    #   return(NULL)  # Exit the function prematurely
    # }
    
    #download EXPOSURE
    
    exposure <- read_exposure_data(filename= add_exposure,
                                   sep=",",#GCST:","
                                   snp_col = "hm_rsid",#GCST:"variant_id"
                                   beta_col = "hm_beta",#GCST:"beta"
                                   se_col = "standard_error",#GCST:"standard_error"
                                   effect_allele_col = "hm_effect_allele",#GCST:"effect_allele"
                                   other_allele_col = "hm_other_allele",#GCST:"other_allele"
                                   pval_col = "p_value",#GCST:"p_value"
                                   eaf_col="hm_effect_allele_frequency", # eaf不是必须的，前面整理进去也可,
                                   # gene_col = "Gene",
                                   samplesize_col = "n",
                                   clump=FALSE
    )
    # Clumping via IEU GWAS database
    exposure_clumped <- ld_clump(dplyr::tibble(rsid=exposure$SNP, pval=exposure$pval.exposure, id=exposure$id.exposure,),
                                 #get_plink_exe()
                                 plink_bin = "[plink_Windows.exe]",
                                 #location of EUR reference gene
                                 bfile = "SNPannotations/data_for_clump/EUR",
                                 clump_kb = 10000,
                                 clump_r2 = 0.01,
                                 clump_p = clump_p_value,
    )
    exposure_extracted <- merge(exposure_clumped, exposure, by.x = "rsid", by.y = "SNP") #extract from exposure
    fwrite(exposure_extracted, file="exposure_extracted.csv") #save
    exposure <- NULL
    exposure_dat <- read_exposure_data(filename= "exposure_extracted.csv",
                                       sep=",",
                                       snp_col = "rsid",
                                       beta_col = "beta.exposure",
                                       se_col = "se.exposure",
                                       effect_allele_col = "effect_allele.exposure",
                                       other_allele_col = "other_allele.exposure",
                                       pval_col = "pval.exposure",
                                       eaf_col='eaf.exposure' , # eaf is non-necessary
                                       samplesize_col = 'samplesize.exposure',
                                       # gene_col = "Gene",
                                       clump=FALSE
    )
    # head(exposure_dat)
    # exposure_dat$samplesize.exposure <- 289268
    exposure_dat$id.exposure <- 'GCST012235'
    exposure_dat$exposure <- 'OSCCandOPSCC'
    # default using EUR population from 1000 Genomes project
    # Read reference for suitable r2
    
    # read OUTCOME
    outcome <- fread(add_outcome,data.table = FALSE) 
    # head(outcome)
    outcome_extracted <- merge(exposure_dat, outcome, by.x = "SNP", by.y = "variant_id") #extract outcome
    outcome <- NULL
    fwrite(outcome_extracted, file="outcome_extracted.csv") #save
    outcome_dat<-read_outcome_data(snps = exposure_dat$SNP, filename= "outcome_extracted.csv",
                                   sep=",",#GCST:","
                                   snp_col = "SNP",#GCST:"variant_id"
                                   beta_col = "beta",#GCST:"beta"
                                   se_col = "standard_error",#GCST:"standard_error"
                                   effect_allele_col = "effect_allele",#GCST:"effect_allele"
                                   other_allele_col = "other_allele",#GCST:"other_allele"
                                   pval_col = "p_value",#GCST:"p_value"
                                   eaf_col="effect_allele_frequency",
                                   # gene_col = "Gene",
                                   samplesize_col = "n",
    ) #read SNP 
    outcome_dat$id.outcome <- substring(file,1,12)
    outcome_dat$outcome <- foldername
    
    dat <- harmonise_data(
      exposure_dat = exposure_dat, 
      outcome_dat = outcome_dat)
    dat <- dat[!is.na(dat$beta.exposure),]
    dat <- dat[!is.na(dat$beta.outcome),]
    dat <- steiger_filtering(dat)
    if (nrow(dat[dat$steiger_dir==FALSE,]) > 0){
      dat[dat$steiger_dir==FALSE,][,"mr_keep"] <- FALSE
    }
    
    r2 <- get_r_from_pn(exposure_dat$pval.exposure,exposure_dat$samplesize.exposure)
    r2_out <- get_r_from_pn(outcome_dat$pval.outcome,outcome_dat$samplesize.outcome)
    dat$r.outcome <- r2_out
    dat$r.exposure <- r2[match(dat$SNP,exposure_dat$SNP)]
    f_stat <- r2 * (exposure_dat$samplesize.exposure - 2) / (1 - r2)
    r2_f_matrix <- cbind(exposure_dat$SNP,r2,f_stat)
    MR_dat <- merge(dat,r2_f_matrix,by.x = "SNP",by.y = "V1")
    fwrite(MR_dat,"MR_dat.csv")
    #MR report (saved as html、word or pdf, default is html)
    mr_report(dat,output_path = result_path,output_type = "html")
    # head(dat)
    presso <- mr_presso(BetaOutcome = "beta.outcome", 
                        BetaExposure = "beta.exposure", 
                        SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE, 
                        DISTORTIONtest = TRUE, 
                        data = dat, 
                        NbDistribution = 5000,  
                        SignifThreshold = 0.05
    )
    presso_1 <- presso$`Main MR results`
    presso_2 <- presso$`MR-PRESSO results`$`Global Test`$RSSobs
    presso_3 <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
    presso_all <- c("RSSobs")
    presso_all <- append(presso_all,presso_2)
    presso_all <- append(presso_all,"Pvalue")
    presso_all <- append(presso_all,presso_3)
    presso_all <- append(presso_1,presso_all)
    write.table(presso_all,"MR_presso.txt")
    print(paste("processing finished:", file,"_", clump_p_value))
  }
  
  # Apply the function in parallel to fileslist
  fileslist <- list.files("ImmuneCells/",pattern = ".gz")
  future_map(fileslist, function(file) {
    print(paste("processing:", file,"_", clump_p_value))
    tryCatch(process_files(file), error = function(e) {
      message(paste("Error occurred for file:", file))
      message("Error message:", conditionMessage(e))
      return(NA)
    })
  })
}


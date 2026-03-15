
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
library(readr)
library(stringr)
library(dplyr)

plan(multisession, workers = 6)#, workers = 3

# Define a function that processes each exposure file
clump_p_values <- c(1e-6, 5e-6, 5e-7, 5e-8) #1e-5, 
clump_p_value <- clump_p_values[1]
exposure_file <- exposure_fileslist[1]
for (clump_p_value in clump_p_values){
  process_exposure <- function(exposure_file) {
    Index_IC <- fread("[ImmuneCells.csv]", data.table = FALSE)
    exposure_foldername <- Index_IC[Index_IC[, 5] == substring(exposure_file, 1, 12), 4]
    name_IC <- fread("[HNC2IC/recorded_table.csv]", data.table = FALSE)
    name_IC <- gsub("_TO_HNC|HNC_TO_","",name_IC[,1])
    if (!(exposure_foldername %in% name_IC)){
      return(NULL)
    }
    # exposure_fileslist <- exposure_fileslist[!grepl("\\.gz", exposure_fileslist)]
    # exposure_fileslist <- exposure_fileslist[!grepl("\\.csv", exposure_fileslist)]
    # exposure_fileslist <- exposure_fileslist[!grepl("\\.xslx", exposure_fileslist)]
    add_exposure <- paste("ImmuneCells/",exposure_file,sep = "")
    pre_result_path <- paste('Results/IC2IF_',as.character(clump_p_value),sep = "")
    if (!dir.exists(pre_result_path)){
      dir.create(pre_result_path)
    }
    setwd(pre_result_path)
    #download EXPOSURE
    a <- fread(add_exposure,data.table = FALSE)
    add_exposure <- paste("ImmuneCells/",substr(exposure_file,1,(nchar(exposure_file)-3)),".tsv",sep = "")
    fwrite(a,add_exposure)
    a <- NULL
    exposure <- read_exposure_data(filename= add_exposure,
                                   sep=",",#GCST:","
                                   snp_col = "variant_id",#GCST:"variant_id"
                                   beta_col = "beta",#GCST:"beta"
                                   se_col = "standard_error",#GCST:"standard_error"
                                   effect_allele_col = "effect_allele",#GCST:"effect_allele"
                                   other_allele_col = "other_allele",#GCST:"other_allele"
                                   pval_col = "p_value",#GCST:"p_value"
                                   eaf_col="effect_allele_frequency",
                                   # gene_col = "Gene",
                                   samplesize_col = "n",
                                   clump=FALSE
    )
    file.remove(add_exposure)
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
    # exposure_dat$samplesize.exposure <- 133251
    exposure_dat$id.exposure <- substring(exposure_file,1,12)
    exposure_dat$exposure <- exposure_foldername
    # default using EUR population from 1000 Genomes project
    # Read reference for suitable r2
    
    outcome_fileslist <- list.files("InflammatoryFactors/",pattern = ".tsv")
    Index_IF <- fread("InflammatoryFactors.csv", data.table = FALSE)
    name_IF <- fread("recorded_table.csv", data.table = FALSE)
    name_IF <- gsub("_TO_HNC|HNC_TO_","",name_IF[,1])
    outcome_file <- outcome_fileslist[1]
    for (outcome_file in outcome_fileslist){
      outcome_foldername <- Index_IF[Index_IF[, 1] == substring(outcome_file, 1, 12), 4]
      #CP2IF: c("C-X-C motif chemokine 10 levels", "C-X-C motif chemokine 9 levels", "Transforming growth factor-alpha levels", "Matrix metalloproteinase-10 levels", "TNF-beta levels", "Interleukin-10 levels")
      #IF2CP: c("Eotaxin levels", "Sulfotransferase 1A1 levels", "Protein S100-A12 levels", "Programmed cell death 1 ligand 1 levels", "Neurturin levels", "Tumor necrosis factor receptor superfamily member 9 levels", "T-cell surface glycoprotein CD6 isoform levels", "Interleukin-4 levels", "Interleukin-10 receptor subunit alpha levels", "Cystatin D levels") 
      if (!(outcome_foldername %in% name_IF)){
        next
      }
      add_outcome <- paste("InflammatoryFactors/",outcome_file,sep = "")
      result_path <- paste('Results/IC2IF_',as.character(clump_p_value),'/',exposure_foldername,'_TO_',outcome_foldername,sep = "")
      if (!file.exists(result_path)){
        dir.create(result_path)
      }
      setwd(result_path)
      # read OUTCOME
      outcome <- fread(add_outcome,data.table = FALSE) 
      head(outcome)
      outcome_extracted <- merge(exposure_dat, outcome, by.x = "SNP", by.y = "rsid") #extract outcome
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
      )
      # outcome_dat$samplesize.outcome <- 263668
      outcome_dat$id.outcome <- substring(outcome_file,1,12)
      outcome_dat$outcome <- outcome_foldername
      
      dat <- harmonise_data(
        exposure_dat = exposure_dat, 
        outcome_dat = outcome_dat)
      
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
      print(paste("processing finished:",  exposure_foldername, "_", outcome_foldername, "_",clump_p_value))
    }
  }
  # Apply the function in parallel to exposure_fileslist
  exposure_fileslist <- list.files("ImmuneCells/",pattern = ".gz")
  future_map(exposure_fileslist, function(exposure_file) {
    print(paste("processing:", exposure_file,"_", clump_p_value))
    tryCatch(process_exposure(exposure_file), error = function(e) {
      message(paste("Error occurred for exposure_file:", exposure_file))
      message("Error message:", conditionMessage(e))
      return(NA)
    })
  })
}

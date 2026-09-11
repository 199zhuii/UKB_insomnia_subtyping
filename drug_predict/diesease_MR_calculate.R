
# ======================== MR Pipeline (Simplified & Modularized) ========================
set.seed(2025)
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(fs)
  library(data.table)
  library(stringr)
  library(vroom)
  library(ieugwasr)
  library(pbapply)
})

# ---- Load data ----
ieugwasr::user()
timestamp <- format(Sys.time(), "%y%m%d_%H%M")
all_res <- paste0("pathway/drug_position/reverse_MR")
dir_create <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# ---- Utility for console output ----
cat_flush <- function(core_id, ...) {
  cat("[Core", core_id, "]", ..., "\n")
  flush.console()
}

# 定义两个新的目录
rsid_res_dir <- file.path(all_res, "rsid")
dir_create(rsid_res_dir)


prot_res_path <- file.path(all_res, "prot")
dir_create(prot_res_path)


source_folder <- "pathway/UKB_ppp/diesease_expose/rsid/"
prot_all  <- read.csv("pathway/drug_protein/all_sig.csv", header = TRUE)  #lm结果中显著的蛋白


protein_names <- prot_all[, 1]  # 调整列索引

all_files <- dir_ls(source_folder, regexp = "*_PROTEIN_pQTL_with_RSID.csv")


files_to_copy <- all_files[sapply(all_files, function(file_path) {
  protein_name <- str_extract(basename(file_path), "^[^_]+")  
  protein_name %in% protein_names
})]


file_copy(files_to_copy, rsid_res_dir)



############################## MR分析 #########################################
iv <- read.csv("pathway/diesease_expose/insomnia_iv.csv", header = TRUE)

rsid_files <- list.files(rsid_res_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
process_protein <- function(file_path, core_id = 1) {
  protein <- str_extract(basename(file_path), "^[^_]+")  
  res_path <- file.path(prot_res_path, protein)
 
  if(file.exists(res_path)){
    cat_flush(paste("SKIP:", protein, "(result exists)"))
    return(NULL)
  }
  
  dir_create(res_path)
  cat_flush(core_id, "Processing protein:", protein)
  
  outcome <- read_outcome_data(file.path(file_path),
                               sep = ",", snp_col = "RSID",
                               beta_col = "BETA", se_col = "SE",
                               effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
                               eaf_col = "A1FREQ", pval_col = "P")
  
  mrdata <- harmonise_data(iv, outcome)
  

  mr_methods <- c(
    "mr_ivw",                      # Inverse variance weighted (默认)
    "mr_egger_regression",         # MR-Egger
    "mr_weighted_median",          # Weighted median
    "mr_weighted_mode",            # Weighted mode
    "mr_simple_mode"               # Simple mode
  )
  

  result_all <- mr(mrdata, method_list = mr_methods)
  
  OR_res <- generate_odds_ratios(result_all)
  
  hetero <- mr_heterogeneity(mrdata)
  pleio  <- mr_pleiotropy_test(mrdata)   
  
  # --- 保存完整结果 ---
  saveRDS(list(
    result_raw = result_all,
    OR = OR_res,
    hetero = hetero,
    pleio = pleio,
    path = res_path
  ), file = file.path(res_path, "result.rds"))
  
  cat_flush(core_id, "Finished protein:", protein)
}
invisible(
  for (i in a:b) {
    rsid_file <- rsid_files[i]
    process_protein(rsid_file)
  }
)



# ========================== RESULT SAVING & SIGNIFICANCE CHECK ==========================
sig_folder <- file.path(all_res, "all/0-SIG")
nsensitive_folder <- file.path(all_res, "all/0-nsensitive")
dir_create(sig_folder); dir_create(nsensitive_folder)

rsd_res <- prot_res_path  

tmp_files <- list.files(rsd_res, pattern = ".rds$", recursive = TRUE, full.names = TRUE)
for (f in tmp_files) {
  res <- readRDS(f)
  with(res, {
    write.csv(result_raw, file.path(path, "MR_result.csv"), row.names = FALSE)
    write.csv(OR, file.path(path, "MR_OR_result.csv"), row.names = FALSE)
    write.csv(hetero, file.path(path, "MR_heter.csv"), row.names = FALSE)
    write.csv(pleio, file.path(path, "MR_pleio.csv"), row.names = FALSE)
  })
}


################  FDR校正P值  ##################

sub_folders <- list.dirs(rsd_res, full.names = TRUE, recursive = FALSE)

all_results <- do.call(rbind, pblapply(sub_folders, function(folder) {
  detection_name <- str_extract(basename(folder), "^[^_]+") 
  csv_file <- file.path(folder, paste0("MR_result.csv"))
  
  if (file.exists(csv_file)) {
    data <- read.csv(csv_file)
    ivw_rows <- data[data$method == "Inverse variance weighted", ]
    
    if (nrow(ivw_rows) > 0) {
      ivw_rows$Detection <- detection_name
      return(ivw_rows)
    }
  }
  return(NULL)
}))
all_results$FDR <- p.adjust(all_results$pval, method = "BH")
write.csv(all_results, file.path(all_res, "MR_result_all_FDR.csv"), row.names = FALSE)

tmp_files <- list.files(rsd_res, pattern = ".rds$", recursive = TRUE, full.names = TRUE)
for (f in tmp_files) {
  if (any(all_results$FDR < 0.05, na.rm = TRUE)) {
    hetero_issue <- any(hetero$Q_pval < 0.05, na.rm = TRUE)
    pleio_issue  <- any(pleio$pval < 0.05, na.rm = TRUE)
    dest <- if (hetero_issue | pleio_issue) nsensitive_folder else sig_folder
    dest_path <- file.path(dest, basename(path))
    if (dir.exists(dest_path)) unlink(dest_path, recursive = TRUE)
    file.rename(path, dest_path)
  }
}


# ---- Copy script for reproducibility ----
script_path <- sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])
if (length(script_path) > 0) {
  file.copy(script_path, file.path(all_res, paste0("MR_script_run_", timestamp, ".R")), overwrite = TRUE)
}
cat("=== MR pipeline completed ===\n")



# ========================== COLOC ANALYSIS FOR PROTEINS ==========================
library(coloc)
# ---------------------- Load outcome GWAS (FinnGen) ----------------------
finn_oc <- vroom("pathway/finn-b-F5_INSOMNIA/finngen_R10_F5_INSOMNIA", 
                 delim = "\t",  
                 show_col_types = FALSE,
                 progress = TRUE)
colnames(finn_oc)[1] <- "CHP"
ncase <- 4801
ncontrol <- 405229
samplesize <- ncase + ncontrol

outcome <- finn_oc %>% 
  dplyr::select(rsids, CHP, pos, alt, ref, af_alt, beta, sebeta, pval) %>%
  dplyr::rename(SNP = rsids,
                chrom = CHP,
                effect_allele = alt,
                other_allele = ref,
                eaf = af_alt,
                se = sebeta,
                P = pval) %>%
  unique() %>%
  mutate(varbeta = se^2,
         MAF = ifelse(eaf < 0.5, eaf, 1 - eaf),
         s = ncase / samplesize,
         z = beta / se) %>%
  na.omit()
outcome$samplesize <- samplesize

# ---------------------- Define proteins ----------------------
subfolders <- list.dirs(sig_folder, full.names = FALSE, recursive = FALSE)

prots <- str_extract(subfolders, "^[^_]+")  
# 去重和清理
prots <- unique(na.omit(prots))  


# ---------------------- Paths ----------------------
coloc_res_dir <- file.path(all_res, "all/coloc_all")
dir.create(coloc_res_dir, recursive = TRUE, showWarnings = FALSE)

rsid_files <- list.files(rsid_res_dir, pattern = "_RSID\\.csv$", full.names = TRUE)

filtered_rsid_files <- rsid_files[sapply(rsid_files, function(file_path) {
  file_name <- basename(file_path)
  protein_name <- gsub("_PROTEIN_pQTL_with_RSID\\.csv$", "", file_name)
  protein_name %in% prots
})]


# ---------------------- COLOC analysis function ----------------------
coloc_protein <- function(file_path, protein_name){
  # Create result folder
  res_file <- file.path(coloc_res_dir, paste0(protein_name, "_snp_pph4.csv"))
  if(file.exists(res_file)){
    cat_flush(paste("SKIP:", protein_name, "(result exists)"))
    return(NULL)
  }
  
  cat_flush(paste("Processing protein:", protein_name))
  
  # ---- Load exposure pQTL ----
  exp <- fread(file_path)
  exp$end <- exp$GENPOS
  colnames(exp) <- c('SNP','CHP','start','ID','other_allele','effect_allele',
                     'eaf','INFO','samplesize','TEST','beta','se','CHISQ','nlogP','P','end')
  exp$varbeta <- exp$se^2
  exp$MAF <- ifelse(exp$eaf < 0.5, exp$eaf, 1 - exp$eaf)
  exp$z <- exp$beta / exp$se
  
  # ---- Identify lead SNP and select +/-1Mb region ----
  lead <- exp %>% arrange(P) %>% slice(1)
  leadchr <- lead$CHP
  leadstart <- lead$start
  leadend <- lead$end
  
  QTLdata <- exp %>%
    filter(CHP == leadchr,
           start > leadstart - 1e6,
           end < leadend + 1e6) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    na.omit()
  
  # Replace P=0 with minimal P
  QTLdata$P[QTLdata$P == 0] <- NA
  min_p <- min(QTLdata$P, na.rm = TRUE)
  QTLdata$P[is.na(QTLdata$P)] <- min_p * 0.01
  
  # ---- Match SNPs with outcome ----
  shared_SNP <- intersect(QTLdata$SNP, outcome$SNP)
  QTLdata <- QTLdata %>% filter(SNP %in% shared_SNP) %>% arrange(SNP)
  GWASdata <- outcome %>% filter(SNP %in% shared_SNP) %>% arrange(SNP)
  
  GWASdata <- GWASdata %>%
    arrange(P) %>%  
    distinct(SNP, .keep_all = TRUE)  
  
  QTLdata <- QTLdata %>% 
    filter(SNP %in% shared_SNP) %>% 
    arrange(SNP) %>%
    distinct(SNP, .keep_all = TRUE)  
  
  # ---- Run coloc ----
  coloc_res <- coloc.abf(
    dataset1 = list(pvalues = GWASdata$P, snp = GWASdata$SNP,
                    type = "cc", s = GWASdata$s[1], N = GWASdata$samplesize[1]),
    dataset2 = list(pvalues = QTLdata$P, snp = QTLdata$SNP,
                    type = "quant", N = QTLdata$samplesize[1]),
    MAF = QTLdata$MAF
  )
  
  # Sort results by posterior probability for H4
  SNP_result <- coloc_res$results %>% arrange(desc(SNP.PP.H4))
  
  # summary data
  summary_result <- as.data.frame(t(coloc_res$summary))
  fwrite(summary_result, file = file.path(coloc_res_dir, paste0(protein_name, "_coloc_summary.csv")))
  # ---- Save result ----
  fwrite(SNP_result, file = res_file)
  cat_flush(paste("Saved COLOC results for", protein_name, "->", res_file))
}

# ---------------------- Run analysis for all proteins ----------------------
for(i in seq_along(filtered_rsid_files)){
  protein_name <- str_remove(basename(filtered_rsid_files[i]), "\\.csv$")
  protein_name <- str_extract(protein_name, "^[^_]+")
  coloc_protein(filtered_rsid_files[i], protein_name)
}

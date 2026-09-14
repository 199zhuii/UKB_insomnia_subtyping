# 加载必要的包
library(dplyr)
library(readr)
library(purrr)
library(stringr)


cis_file <- "pathway/UKB_ppp/all_UKB_cis_protein.csv"  
rsid_folder <- "pathway/UKB_ppp/protein_expose/rsid/"  
trans_folder <- "pathway/UKB_ppp/protein_expose/trans/"  

cis_info <- read_csv(cis_file)


filter_trans_pqtl <- function(protein_name, cis_chr, cis_start, cis_end) {
  b_file_name <- sprintf("%s_PROTEIN_pQTL_with_RSID.csv", protein_name)
  b_file_path <- file.path(rsid_folder, b_file_name)
  
  cat(sprintf("dealng protein: %s\n", protein_name))
  
  tryCatch({
    pqtl_data <- read_csv(b_file_path)
    
    trans_pqtl <- pqtl_data %>%
      filter(
        !(CHROM == cis_chr & GENPOS >= cis_start & GENPOS <= cis_end)
      )
    
    n_total <- nrow(pqtl_data)
    n_trans <- nrow(trans_pqtl)
    
    cat(sprintf("  原始SNP数: %d, trans-pQTL数: %d (%.1f%%)\n", 
                n_total, n_trans, ifelse(n_total>0, n_trans/n_total*100, 0)))
    
    if (n_trans == 0) {
      warning(sprintf("蛋白 %s 没有找到trans-pQTL", protein_name))
      return(NULL)
    }
    
 
    trans_pqtl$protein_name <- protein_name
    
    return(trans_pqtl)
    
  }, error = function(e) {
    warning(sprintf("读取文件 %s 时出错: %s", rsid_folder, e$message))
    return(NULL)
  })
}

results <- list()
success_count <- 0

for (i in 1:nrow(cis_info)) {
  protein <- cis_info$protein[i]
  chr <- cis_info$chromosome_name[i]
  start <- cis_info$cis_start[i]
  end <- cis_info$cis_end[i]
  
  trans_result <- filter_trans_pqtl(protein, chr, start, end)
  
  if (!is.null(trans_result)) {
    results[[protein]] <- trans_result
    
    output_file <- file.path(trans_folder, sprintf("%s_trans_pQTL_with_RSID.csv", protein))
    write_csv(trans_result, output_file)
    
    success_count <- success_count + 1
    cat(sprintf("  结果已保存: %s\n", basename(output_file)))
  }
  
  cat("---\n")
}

library(tidyverse)
library(ggpubr)
library(car)         
library(multcomp)     
library(rstatix)      
library(HH)       
library(broom)
library(emmeans)     



# 读取数据
df <- read.csv("pathway/validation/protein_val.csv",na.strings = "NaN")
df_cov_time <- read.csv("pathway/validation/tgap_table_pro_val.csv",na.strings = "NaN")

all_cols <- colnames(df)
start_physio <- which(all_cols == "A1BG")
end_physio <- which(all_cols == "ZPR1")

all_cols <- colnames(df)

start_covars <- which(all_cols == "Batch")
physio_columns <- all_cols[start_physio:end_physio]
covariate_columns <- c(all_cols[start_covars:length(all_cols)])


COV_data <- df[start_covars:length(all_cols)]

missing_values <- colSums(is.na(COV_data))

for (col in colnames(COV_data)) {
  median_value <- median(COV_data[[col]], na.rm = TRUE) # 计算列的均值，忽略NA
  COV_data[[col]][is.na(COV_data[[col]])] <- median_value # 填补缺失值
}


df[start_covars:length(all_cols)] <- COV_data

df <- df %>% 
  mutate(subtype = factor(subtype),
         sex = factor(sex),
         ethic = factor(ethic),
         Batch = factor(Batch),
         site = factor(site)) %>% 
  mutate(subtype = relevel(subtype, ref = "0"))  # 设置subtype0为参照

results <- data.frame()


for (i in 1:length(physio_columns)) {
  Y = physio_columns[i]
  col_i <- df_cov_time[Y]
  temp <- cbind(df, col_i)
  colnames(temp)[ncol(temp)] <- paste0(Y, "_tgap")

  idx <- which(!is.na(temp[Y]))

  temp <- temp[idx, ]
  covariate_columns <- c(colnames(temp)[start_covars:length(temp)])
  Independent_Item = paste(covariate_columns, collapse = " + ")

  formula_2 = paste(Y," ~ subtype + ",paste(Independent_Item))

  model <- lm(formula_2, data = temp)
  

  emmeans_model <- emmeans(model, ~ subtype)
  contrast_results <- contrast(emmeans_model, method = "trt.vs.ctrl", ref = "subtype0",adjust = 'none')
  contrast_df <- as.data.frame(contrast_results)

  mean_control = summary(emmeans_model)$emmean[1]
  

  contrast_df <- contrast_df %>%
    mutate(
      biomarker = Y,
      comparison = case_when(
        contrast == "subtype1 - subtype0" ~ "subtype1_vs_0",
        contrast == "subtype2 - subtype0" ~ "subtype2_vs_0", 
        contrast == "subtype3 - subtype0" ~ "subtype3_vs_0"
      ),
      mean_control = mean_control,
      cohens_d = estimate / sigma(model)
    )
  
  results <- bind_rows(results, contrast_df)
  # 打印进度
  if (i %% 100 == 0) {
    cat(sprintf("Processed %d out of %d proteins\n", i, length(physio_columns)))
  }
}

# 多重检验校正
results <- results %>%
  mutate(
    p_adj_fdr = p.adjust(p.value, method = "fdr"),
    p_adj_bonferroni = p.adjust(p.value, method = "bonferroni"),
    significance = case_when(
      p_adj_fdr < 0.001 ~ "***",
      p_adj_fdr < 0.01 ~ "**",
      p_adj_fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()

# 保存结果
write.csv(results, "pathway/validation/lm_ref_health_pro_val.csv", row.names = FALSE)



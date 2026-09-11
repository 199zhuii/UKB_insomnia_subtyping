library(tidyverse)
library(ggpubr)
library(car)
library(multcomp)
library(rstatix)
library(HH)
library(broom)
library(emmeans)
library(forestplot)

############################################  sleep behavior   ###############################################################
# 1. 数据准备 ----------------------------------------------------------------
df <- read_csv("pathway/sleep_behavior.csv")


all_cols <- colnames(df)
# 找到生理指标的开始和结束位置
start_physio <- which(all_cols == "BC_CCI")
end_physio <- which(all_cols == "Nap")

# 找到协变量的开始位置
start_covars <- which(all_cols == "sex")
physio_columns <- all_cols[start_physio:end_physio]
covariate_columns <- c(all_cols[start_covars:length(all_cols)])


# 处理协变量
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
         ethic = factor(ethic))%>% 
  mutate(subtype = relevel(subtype, ref = "0"))  # make the healthy as ref

# 组变量
subtype <-df$subtype

# 2. 计算前提假设 ----------------------------------------------------------------
attach(df)

results <- data.frame()

for (i in 1:length(physio_columns)) {
  Y = physio_columns[i]
  
  temp <- df %>% 
    dplyr::select(all_of(Y), all_of(covariate_columns), subtype) %>%
    filter(!is.na(.[[Y]]))
  
  Independent_Item = paste(covariate_columns, collapse = " + ")
  formula = paste(Y, " ~ ", paste(Independent_Item, " + subtype"))
  
  model <- lm(formula, data = temp)

  emmeans_model <- emmeans(model, ~ subtype)
  
  contrast_results <- contrast(emmeans_model, method = "trt.vs.ctrl", ref = "subtype0")
  contrast_df <- as.data.frame(contrast_results)
  
  mean_control <- summary(emmeans_model)$emmean[1]
  

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
    ) %>%
    rename(
      mean_diff = estimate,
      se_diff = SE,
      df_diff = df,
      t_value = t.ratio,
      p_value = p.value
    )
  
  results <- bind_rows(results, contrast_df)
}


results <- results %>%
  group_by(comparison) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr"),
    p_adj_bonferroni = p.adjust(p_value, method = "bonferroni"),
    significance_fdr = case_when(
      p_adj_fdr < 0.001 ~ "***",
      p_adj_fdr < 0.01 ~ "**",
      p_adj_fdr < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()
# 保存
write.csv(results, file = "pathway/lm_ref_health_sleep_behavior.csv", row.names = FALSE)




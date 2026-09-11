library(tidyverse)
library(ggpubr)
library(car)          
library(multcomp)     
library(rstatix)      
library(HH)       
library(broom)
library(emmeans)      


results <- data.frame()
# 1. 
df <- read_csv("pathway/cognitive_test.csv")
all_cols <- colnames(df)


start_physio <- which(all_cols == "Fluid_intelligence_score")
end_physio <- which(all_cols == "Matrix_pattern_completion")


start_covars <- which(all_cols == "sex")
physio_columns <- all_cols[start_physio:end_physio]
covariate_columns_1 <- c(all_cols[start_covars:length(all_cols)])


COV_data <- df[start_covars:length(all_cols)]
missing_values <- colSums(is.na(COV_data))
for (col in colnames(COV_data)) {
  median_value <- median(COV_data[[col]], na.rm = TRUE) 
  COV_data[[col]][is.na(COV_data[[col]])] <- median_value 
}


df[start_covars:length(all_cols)] <- COV_data
df <- df %>% 
  mutate(subtype = factor(subtype),
         sex = factor(sex),
         ethic = factor(ethic),
         EDU = factor(EDU))%>% 
  mutate(subtype = relevel(subtype, ref = "0"))  # 
# 组变量
subtype <-df$subtype


# 2. 
attach(df)

for (i in 1:length(physio_columns)) {
  Y = physio_columns[i]
  pattern <- paste0(".*", Y, ".*")
  matching_cols <- str_subset(covariate_columns_1, regex(pattern, ignore_case = TRUE))
  covariate_columns <- c(covariate_columns_1[c(1:5)],matching_cols)
  idx <- which(!is.na(df[Y]))
  temp <- df[idx, ]
  
  #
  Y_mean <- mean(temp[[Y]], na.rm = TRUE)
  Y_sd <- sd(temp[[Y]], na.rm = TRUE)
  
  # 
  Y_z <- paste0(Y, "_z")
  
  # 
  temp[[Y_z]] <- (temp[[Y]] - Y_mean) / Y_sd
  
  Independent_Item = paste(covariate_columns, collapse = " + ")
  # 
  formula = paste(Y_z," ~ ",paste(Independent_Item," + subtype"))
  # 
  model <- lm(formula, data = temp)
  
  # 
  emmeans_model <- emmeans(model, ~ subtype)
  
  # 
  contrast_results <- contrast(emmeans_model, method = "trt.vs.ctrl", ref = "subtype0",adjust = 'none')
  contrast_df <- as.data.frame(contrast_results)
  
  # 
  mean_control <- summary(emmeans_model)$emmean[1]
  
  # 
  contrast_df <- contrast_df %>%
    mutate(
      biomarker = Y,
      comparison = case_when(
        contrast == "subtype1 - subtype0" ~ "subtype1_vs_0",
        contrast == "subtype2 - subtype0" ~ "subtype2_vs_0", 
        contrast == "subtype3 - subtype0" ~ "subtype3_vs_0"
      ),
      mean_control = mean_control,
      # Cohen's d
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

# correct p value
results <- results %>%
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

# 
write.csv(results, "pathway/lm_ref_health_cog.csv", row.names = FALSE)


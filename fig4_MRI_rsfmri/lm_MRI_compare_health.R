library(tidyverse)
library(ggpubr)
library(car)          
library(multcomp)      
library(rstatix)      
library(HH)       
library(broom)
library(emmeans)      



############################################  all_the_brain   ###############################################################
results <- data.frame()
# 1. 
df <- read_csv("pathway/All_T1_MRI.csv")
# 
all_cols <- colnames(df)
start_physio <- which(all_cols == "lh_bankssts_volume")
end_physio <- which(all_cols == "Right_VentralDC")


start_covars <- which(all_cols == "instance2")
physio_columns <- all_cols[start_physio:end_physio]
covariate_columns <- c(all_cols[start_covars:length(all_cols)])


# COV
COV_data <- df[start_covars:length(all_cols)]
missing_values <- colSums(is.na(COV_data))
for (col in colnames(COV_data)) {
  median_value <- median(COV_data[[col]], na.rm = TRUE) # 
  COV_data[[col]][is.na(COV_data[[col]])] <- median_value # 
}

#
df[start_covars:length(all_cols)] <- COV_data
df <- df %>% 
  mutate(subtype = factor(subtype),
         sex = factor(sex),
         ethic = factor(ethic),
         MRI_site_instance2 = factor(MRI_site_instance2),
         EDU = factor(EDU))%>% 
  mutate(subtype = relevel(subtype, ref = "0"))  # 
subtype <-df$subtype


# 2. 
attach(df)

for (i in 1:length(physio_columns)) {
  Y = physio_columns[i]
  idx <- which(!is.na(df[Y]))
  temp <- df[idx, ]
  Y_mean <- mean(temp[[Y]], na.rm = TRUE)
  Y_sd <- sd(temp[[Y]], na.rm = TRUE)
  Y_z <- paste0(Y, "_z")
  
  temp[[Y_z]] <- (temp[[Y]] - Y_mean) / Y_sd
  Independent_Item = paste(covariate_columns, collapse = " + ")
  formula_2 = paste(Y_z," ~ ",paste(Independent_Item," + subtype"))
  
  reduced_model <- lm(formula_2, data = temp)
  model_summary <- tidy(reduced_model, conf.int = TRUE)
  
  emmeans_model <- emmeans(reduced_model, ~ subtype)
  contrast_results <- contrast(emmeans_model, method = "trt.vs.ctrl", ref = "subtype0",adjust = 'none')
  contrast_df <- as.data.frame(contrast_results)
  mean_control = summary(emmeans_model)$emmean[1]

  contrast_df <- contrast_df %>%
    mutate(
      biomarker = Y,
      comparison = case_when(
        contrast == "subtype1 - subtype0" ~ "subtype1_vs_0",
        contrast == "subtype2 - subtype0" ~ "subtype2_vs_0",
        contrast == "subtype3 - subtype0" ~ "subtype3_vs_0"))  %>%
    mutate(
      mean_control = mean_control,
      sigma_model = sigma(reduced_model),
      log_P = log(p.value),
      # Cohen's d
      cohens_d = estimate / sigma(reduced_model),
      CI_upper = estimate + 1.96*SE,
      CI_lower = estimate - 1.96*SE
    )
  results <- bind_rows(results, contrast_df)
}

# Correct p value
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

# 
write.csv(results, "pathway/lm_ref_health_all_brain_T1.csv", row.names = FALSE)


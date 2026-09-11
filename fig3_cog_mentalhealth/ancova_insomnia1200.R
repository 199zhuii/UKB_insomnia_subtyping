# 加载包
library(tidyverse)
library(MASS)
library(emmeans)
library(car)

# 1. 数据准备 ----------------------------------------------------------------
df <- read_csv("pathway/insomnia 1200/mergedData_insomnia.csv")


all_cols <- colnames(df)

start_covars <- which(all_cols == "instance0")

physio_columns <- all_cols[3]
covariate_columns <- c(all_cols[start_covars:length(all_cols)])


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

df$insomnia_factor <- factor(df[[physio_columns]], ordered = TRUE)

subtype <-df$subtype


# 2. 计算前提假设 ----------------------------------------------------------------
attach(df)

for (i in 1:length(physio_columns)) {
  Y = physio_columns[i]
  idx <- which(!is.na(df[Y]))
 
  temp <- df[idx, ]
  
 
  Independent_Item = paste(covariate_columns, collapse = " + ")
  formula_2 = paste("insomnia_factor ~ subtype + ",paste(Independent_Item))
  
  model_ord <- polr(formula_2, data = temp, method = "logistic", Hess = TRUE)
  
  formula_reduced <- as.formula(paste("insomnia_factor ~ ", paste(covariate_columns, collapse = " + ")))
  model_reduced <- polr(formula_reduced, data = temp, method = "logistic")
  overall_anova <- anova(model_reduced, model_ord, test = "Chisq")
  overall_p <- overall_anova$'Pr(Chi)'[2]
  
  emm <- emmeans(model_ord, ~ subtype)
  pairs_OR <- pairs(emm, adjust = "tukey", type = "response")
  posthoc_df <- as.data.frame(summary(pairs_OR))
  
  overall_result <- data.frame(Variable = physio_columns, Overall_p = overall_p)
}

# 保存
write.csv(overall_result, file = "pathway/insomnia 1200/insomnia_1200_ordinal_overall.csv", row.names = FALSE)
write.csv(posthoc_df, "pathway/insomnia 1200/insomnia_1200_ordinal_posthoc.csv", row.names = FALSE)





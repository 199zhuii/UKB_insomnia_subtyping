# 加载必要包
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyverse)
library(broom.mixed)
library(performance)
library(patchwork)


for (var_name in c("Alcohol_Use", "Anxiety", "Depression", "Mania","Well_being")) {
  
  df_data <- read_csv(paste0("pathway/mental health/",var_name,".csv"))

  all_cols <- colnames(df_data)
  
  base_score <- all_cols[3]
  follow_score <- all_cols[4]
  
  start_covars <- which(all_cols == "sex")
  covariate_columns <- c(all_cols[start_covars:length(all_cols)])
  
  COV_data <- df_data[start_covars:length(all_cols)]

  missing_values <- colSums(is.na(COV_data))

  for (col in colnames(COV_data)) {
    median_value <- median(COV_data[[col]], na.rm = TRUE) # 
    COV_data[[col]][is.na(COV_data[[col]])] <- median_value # 
  }
  
  df_data[start_covars:length(all_cols)] <- COV_data
  
  df <- df_data %>%
    pivot_longer(
      cols = c(all_of(base_score), all_of(follow_score)),
      names_to = "timepoint",
      values_to = "score"
    ) %>%
    mutate(
      time = ifelse(str_detect(timepoint, "baseline"), 0, 1))
  
  # 
  df <- df %>% 
    mutate(subtype = factor(subtype),
           sex = factor(sex),
           ethic = factor(ethic),
           EDU = factor(EDU))
  # 
  subtype <-df$subtype
  
  # 2. 
  ########################################## main anlysis###############################
  cat("ANALYZING COGNITIVE TEST:", var_name, "\n")

  # 
  Independent_Item = paste(covariate_columns, collapse = " + ")
  formula = paste("score ~ subtype * time + ", paste(Independent_Item, " + (1 | eid)"))
  model <- lmer(formula, data = df)
  
  model_summary <- summary(model)
  anova_results <- anova(model)
  
  #
  interaction_p <- anova_results["subtype:time", "Pr(>F)"]
  
  #
  simple_effects_df <- NULL
  posthoc_tukey_df <- NULL
  
  if(interaction_p < 0.05) {
    cat("\n--- SIMPLE EFFECTS ANALYSIS ---\n")
    
    #
    cat("Rate of Change by Subtype:\n")
    slopes <- emtrends(model, ~subtype, var = "time")
    slopes_contrasts <- pairs(slopes, adjust = "fdr")
    
    slopes_summary <- as.data.frame(slopes)
    slope_rows <- data.frame(slopes_summary)
    
    # 
    cat("\nGroup Differences at Each Timepoint:\n")
    time_effects <- emmeans(model, ~subtype | time)
    time_contrasts <- contrast(time_effects, "pairwise", adjust = "fdr")
    
    simple_effects <- list(
      slopes_contrasts = as.data.frame(slopes_contrasts),
      time_contrasts = as.data.frame(time_contrasts))
    # 
    slopes_df <- simple_effects$slopes_contrasts %>%
      mutate(
        test = xx,
        contrast_type = "slope_difference",
        adjustment_method = "FDR"
      )
    # 
    time_contrasts_df <- simple_effects$time_contrasts %>%
      mutate(
        test = xx,
        contrast_type = "timepoint_difference",
        adjustment_method = "FDR"
      )
    simple_effects_df <- bind_rows(slopes_df, time_contrasts_df)
  }
  
  # 
  cat("Tukey-adjusted pairwise comparisons across all groups and timepoints:\n")
  all_pairwise <- emmeans(model, ~ subtype * time) %>%
    pairs(by = NULL, adjust = "tukey")  # by = NULL 
  
  print(all_pairwise)
  
  # 
  posthoc_tukey_df <- as.data.frame(all_pairwise) %>%
    mutate(
      test = xx,
      contrast_type = "pairwise_comparison",
      adjustment_method = "Tukey"
    )
  # 
  significant_tukey <- posthoc_tukey_df %>%
    filter(p.value < 0.05)
  
  # 
  model_estimates <- tidy(model, effects = "fixed")
  
  # 
  results_df <- data.frame(
    test = xx,
    interaction_F = anova_results["subtype:time", "F value"],
    interaction_df_num = anova_results["subtype:time", "NumDF"],
    interaction_df_den = anova_results["subtype:time", "DenDF"],
    interaction_p = interaction_p,
    sample_size = n_distinct(df$eid),
    n_observations = nrow(df),
    n_significant_tukey = ifelse(!is.null(posthoc_tukey_df), 
                                 sum(posthoc_tukey_df$p.value < 0.05, na.rm = TRUE), 0),
    stringsAsFactors = FALSE
  )
  
  # 
  coefficients_df <- tidy(model, effects = "fixed")
  coefficient_results <- coefficients_df %>%
    mutate(
      test = xx,
      term_clean = case_when(
        term == "(Intercept)" ~ "Intercept",
        str_detect(term, "subtype") & str_detect(term, "time") ~ "Group_Time_Interaction",
        str_detect(term, "subtype") ~ "Group_Effect",
        term == "time" ~ "Time_Effect",
        TRUE ~ term
      )
    ) %>%
    select(test, term, term_clean, estimate, std.error, statistic, p.value)
  

  write_csv(paste0("pathway/base_follow/main_", var_name, ".csv"))
  write_csv(coefficient_results, paste0("pathway/base_follow/coefficients_", var_name, ".csv"))
  write_csv(slope_rows, paste0("pathway/base_follow/slope_", var_name, ".csv"))
  
 
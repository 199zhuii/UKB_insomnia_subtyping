# 1. 加载必要的包
library(tidyverse)
library(lubridate)
library(gtsummary)
library(forestplot)
library(readxl)
library(broom)
library(car)  # 用于VIF检查
library(performance)  # 用于模型诊断
library(patchwork)  # 用于图形组合
library(logistf)
library(dplyr)
library(tidyr)

disoer_id = c('Affective-Disorders','All-Cause-Dementia','Anxiety-and-Related-Disorders','Cardiovascular-disease','Sleep-Disorders','Parkinsonism','Autonomic-nervous-system-disorder')

for (i in 1:length(disoer_id)){
  # 正确读取和处理数据
  disorder_name = disoer_id[i]
  data <- read_excel(paste0("pathway/ICD/ICD merge/",disorder_name,".xlsx")) 

  covariate_cols <- c("sex", "ethic", "sleep_scale_age", "smoking", "alcohol", 
                      "BMI", "EDU", "TDI", "PCA1", "PCA2", "PCA3", "PCA4", 
                      "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10", 
                      "PCA11", "PCA12", "PCA13", "PCA14", "PCA15", "PCA16", 
                      "PCA17", "PCA18", "PCA19", "PCA20")
  
  data <- data %>%
    mutate(across(all_of(covariate_cols), 
                  ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
  
# 3. 数据清洗与预处理
data <- data %>%
  mutate(
    # 定义失眠分组因子
    subtype = factor(subtype),
    # 定义性别因子
    sex = factor(sex),
    # 定义民族因子（根据UKB编码）
    ethic = factor(ethic),
    # 定义吸烟状态因子
    smoking = factor(smoking),
    # 定义饮酒状态因子，不是
    #alcohol = factor(alcohol),
    # 定义教育水平因子
    EDU = factor(EDU),
    has_disease = ifelse((disorder == 1), 1, 0),
    # 计算诊断到问卷的时间间隔（年）
    years_diagnosis_to_questionnaire = ifelse(!is.na(disorder_age),
                                              sleep_scale_age - disorder_age,
                                              NA_real_))


# 4. 患病率分析
prevalence_by_group <- data %>%
  group_by(subtype) %>%
  summarise(
    n_total = n(),
    n_cases = sum(has_disease),
    prevalence = n_cases / n_total * 100,
    .groups = "drop"
  )

cat("\n=== 各组疾病患病率 ===\n")
print(prevalence_by_group)



# 5. 多因素逻辑回归分析
multivariate_model <- glm(
  has_disease ~ subtype + sex + ethic + sleep_scale_age + 
    BMI + smoking + alcohol + EDU + TDI +
    PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+PCA11+
    PCA12+PCA13+PCA14+PCA15+PCA16+PCA17+PCA18+PCA19+PCA20,
  data = data,
  family = binomial()
)

# 模型摘要
cat("\n=== 多因素逻辑回归模型摘要 ===\n")
print(summary(multivariate_model))

# 提取OR和95% CI
multivariate_results <- tidy(multivariate_model, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(
    term_clean = case_when(
      str_detect(term, "subtype") ~ gsub("subtype", "", term),
      str_detect(term, "sex") ~ gsub("sex", "", term),
      str_detect(term, "ethic") ~ gsub("ethic", "", term),
      str_detect(term, "smoking") ~ gsub("smoking", "", term),
      str_detect(term, "alcohol") ~ gsub("alcohol", "", term),
      str_detect(term, "EDU") ~ gsub("EDU", "", term),
      TRUE ~ term
    ),
    OR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 提取关键结果
all_results <- tibble(
  term = multivariate_results$term,
  OR = multivariate_results$estimate,
  lower = multivariate_results$conf.low,
  upper = multivariate_results$conf.high,
  p_value = multivariate_results$p.value
) %>%
  filter(grepl("subtype", term)) %>%  # 只保留亚型比较
  mutate(
    comparison = case_when(
      term == "subtype1" ~ "Subtype 1 vs 0",
      term == "subtype2" ~ "Subtype 2 vs 0",
      term == "subtype3" ~ "Subtype 3 vs 0",
      TRUE ~ term)
  ) %>%
  select(comparison, OR, lower, upper, p_value)


  
  cat("\n学术论文表格:\n")
  print(all_results)
  
  savepath = 'pathway/ICD/stat_result/'
  # 保存结果
  write_csv(all_results, paste0(savepath,disorder_name,"_glm_results.csv"))
} 







disoer_id = c('Affective-Disorders','All-Cause-Dementia','Anxiety-and-Related-Disorders','Cardiovascular-disease','Sleep-Disorders','Parkinsonism','Autonomic-nervous-system-disorder')
for (i in 1:length(disoer_id)){
  # 正确读取和处理数据
  disorder_name = disoer_id[i]
  data <- read_excel(paste0("pathway/ICD/ICD merge/",disorder_name,".xlsx")) 
  
  covariate_cols <- c("sex", "ethic", "sleep_scale_age", "smoking", "alcohol", 
                      "BMI", "EDU", "TDI", "PCA1", "PCA2", "PCA3", "PCA4", 
                      "PCA5", "PCA6", "PCA7", "PCA8", "PCA9", "PCA10", 
                      "PCA11", "PCA12", "PCA13", "PCA14", "PCA15", "PCA16", 
                      "PCA17", "PCA18", "PCA19", "PCA20")
  
  data <- data %>%
    mutate(across(all_of(covariate_cols), 
                  ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
  # 3. 数据清洗与预处理
  data <- data %>%
    mutate(
      # 定义失眠分组因子
      subtype = factor(subtype),
      # 定义性别因子
      sex = factor(sex),
      # 定义民族因子（根据UKB编码）
      ethic = factor(ethic),
      # 定义吸烟状态因子
      smoking = factor(smoking),
      # 定义饮酒状态因子，不是
      #alcohol = factor(alcohol),
      # 定义教育水平因子
      EDU = factor(EDU),
      has_disease = ifelse((disorder == 1), 1, 0),
      # 计算诊断到问卷的时间间隔（年）
      years_diagnosis_to_questionnaire = ifelse(!is.na(disorder_age),
                                                sleep_scale_age - disorder_age,
                                                NA_real_))

  
  
  # 4. 多因素逻辑回归分析
  firth_model <- logistf(
    has_disease ~ subtype + sex + ethic + sleep_scale_age + BMI + smoking + alcohol + EDU + TDI +
      PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+PCA11+
      PCA12+PCA13+PCA14+PCA15+PCA16+PCA17+PCA18+PCA19+PCA20,
    data = data,
    control = logistf.control(
      maxit = 5000,      # 增加最大迭代次数
      maxstep = 0.1,     # 减少步长以提高稳定性
    ),
    plcontrol = logistpl.control(
      maxit = 10000,     # 显著增加PL置信区间的迭代次数
      maxstep = 0.05,    # 减少步长
    )
  )
  
  # 模型摘要
  cat("\n=== 多因素逻辑回归模型摘要 ===\n")
  print(summary(firth_model))
  
  # 提取关键结果
  firth_results <- tibble(
    term = names(firth_model$coefficients),
    OR = exp(firth_model$coefficients),
    lower = exp(firth_model$ci.lower),
    upper = exp(firth_model$ci.upper),
    p_value = firth_model$prob
  ) %>%
    filter(grepl("subtype", term)) %>%  # 只保留亚型比较
    mutate(
      comparison = case_when(
        term == "subtype1" ~ "Subtype 1 vs 0",
        term == "subtype2" ~ "Subtype 2 vs 0",
        term == "subtype3" ~ "Subtype 3 vs 0",
        TRUE ~ term)
    )
  
  # 提取关键结果
  all_results <- tibble(
    term = firth_results$term,
    OR = firth_results$OR,
    lower = firth_results$lower,
    upper = firth_results$upper,
    p_value = firth_results$p_value
  ) %>%
    filter(grepl("subtype", term)) %>%  # 只保留亚型比较
    mutate(
      comparison = case_when(
        term == "subtype1" ~ "Subtype 1 vs 0",
        term == "subtype2" ~ "Subtype 2 vs 0",
        term == "subtype3" ~ "Subtype 3 vs 0",
        TRUE ~ term)
    ) %>%
    select(comparison, OR, lower, upper, p_value)
  
  
  
  cat("\n学术论文表格:\n")
  print(all_results)
  
  savepath = 'pathway/ICD/stat_result/'
  # 保存结果
  write_csv(all_results, paste0(savepath,disorder_name,"_firth_results.csv"))
} 




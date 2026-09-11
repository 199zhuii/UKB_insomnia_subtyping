# ============================================================
# BCH法评估LCA各条目对分类的贡献
# 需要：后验概率（第28-30列）和硬分类结果（第31列）
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(effectsize)
library(purrr)

# 读取数据
df <- read_csv("pathway/insomnia_with_prob_for_BCH.csv")  # 替换为你的文件名

# 重命名关键列（根据你的实际列名调整）
eid_col <- names(df)[1]           # eid列
item_cols <- names(df)[2:29]      # 第2-27列：26个LCA条目
prob_cols <- names(df)[30:32]     # 第28-30列：后验概率
class_col <- names(df)[33]        # 第31列：所属类别


df_clean <- df %>%
  rename(
    eid = all_of(eid_col),
    class = all_of(class_col),
    prob_c1 = all_of(prob_cols[1]),
    prob_c2 = all_of(prob_cols[2]),
    prob_c3 = all_of(prob_cols[3])
  )


# ============================================================
# BCH加权均值比较
# ============================================================

bch_weighted_comparison <- function(data, item_name, prob_vars) {
  
 
  item_values <- data[[item_name]]
  probs <- as.matrix(data[, prob_vars])
  
 
  probs <- probs / rowSums(probs)
  
  n_classes <- ncol(probs)
  n <- nrow(data)
  
  # === 步骤1：计算BCH加权均值 ===
  bch_means <- numeric(n_classes)
  for (k in 1:n_classes) {
    weights <- probs[, k]
    weights_norm <- weights / sum(weights)
    bch_means[k] <- sum(weights_norm * item_values, na.rm = TRUE)
  }
  
  # === 步骤2：计算BCH加权方差和标准误 ===
  bch_var <- numeric(n_classes)
  bch_se <- numeric(n_classes)
  
  for (k in 1:n_classes) {
    weights <- probs[, k]
    weights_norm <- weights / sum(weights)
    
    weighted_var <- sum(weights_norm * (item_values - bch_means[k])^2, na.rm = TRUE)
    n_eff <- sum(weights)^2 / sum(weights^2)
    bch_var[k] <- weighted_var / n_eff
    bch_se[k] <- sqrt(bch_var[k])
  }
  
  # === 步骤3：整体Wald检验 ===
  # 使用加权方差分析
  grand_mean <- mean(item_values, na.rm = TRUE)
  
  ss_between <- 0
  for (k in 1:n_classes) {
    ss_between <- ss_between + sum(probs[, k]) * (bch_means[k] - grand_mean)^2
  }
  
  ss_within <- 0
  for (k in 1:n_classes) {
    weights <- probs[, k]
    ss_within <- ss_within + sum(weights * (item_values - bch_means[k])^2, na.rm = TRUE)
  }
  
  df_between <- n_classes - 1
  df_within <- n - n_classes
  
  ms_between <- ss_between / df_between
  ms_within <- ss_within / df_within
  F_stat <- ms_between / ms_within
  
  p_value <- pf(F_stat, df_between, df_within, lower.tail = FALSE)
  
  # === 步骤4：计算效应量（加权η²） ===
  eta_sq <- ss_between / (ss_between + ss_within)
  
  # === 步骤5：两两比较的Cohen's d ===
  # 使用加权均值计算各类别间的d值
  pair_d <- matrix(NA, n_classes, n_classes)
  colnames(pair_d) <- paste0("Class", 1:n_classes)
  rownames(pair_d) <- paste0("Class", 1:n_classes)
  
  for (k1 in 1:n_classes) {
    for (k2 in 1:n_classes) {
      if (k1 < k2) {
        pooled_sd <- sqrt((bch_var[k1] * sum(probs[,k1]) + bch_var[k2] * sum(probs[,k2])) / 
                            (sum(probs[,k1]) + sum(probs[,k2])))
        pair_d[k1, k2] <- (bch_means[k1] - bch_means[k2]) / pooled_sd
        pair_d[k2, k1] <- -pair_d[k1, k2]
      }
    }
  }
  
  # === 返回结果 ===
  result <- list(
    item = item_name,
    bch_means = bch_means,
    bch_se = bch_se,
    F_statistic = F_stat,
    df1 = df_between,
    df2 = df_within,
    p_value = p_value,
    eta_squared = eta_sq,
    cohens_d_matrix = pair_d,
    max_cohens_d = max(abs(pair_d[upper.tri(pair_d)]))
  )
  
  return(result)
}

# ============================================================
# 对所有条目执行BCH分析
# ============================================================

cat("\n========== 开始BCH分析 ==========\n")
cat("分析条目数:", length(item_cols), "\n")


bch_results <- list()
pb <- txtProgressBar(min = 0, max = length(item_cols), style = 3)

for (i in seq_along(item_cols)) {
  bch_results[[i]] <- bch_weighted_comparison(
    data = df_clean,
    item_name = item_cols[i],
    prob_vars = c("prob_c1", "prob_c2", "prob_c3")
  )
  setTxtProgressBar(pb, i)
}
close(pb)

# ============================================================
# 提取并汇总结果
# ============================================================

summary_df <- data.frame(
  item = sapply(bch_results, `[[`, "item"),
  F_statistic = sapply(bch_results, `[[`, "F_statistic"),
  df1 = sapply(bch_results, `[[`, "df1"),
  df2 = sapply(bch_results, `[[`, "df2"),
  p_value = sapply(bch_results, `[[`, "p_value"),
  eta_squared = sapply(bch_results, `[[`, "eta_squared"),
  max_cohens_d = sapply(bch_results, `[[`, "max_cohens_d"),
  stringsAsFactors = FALSE
)


summary_df$p_adjusted <- p.adjust(summary_df$p_value, method = "bonferroni")
summary_df$p_fdr <- p.adjust(summary_df$p_value, method = "fdr")


summary_df$significance <- ifelse(summary_df$p_fdr < 0.001, "***",
                                  ifelse(summary_df$p_fdr < 0.01, "**",
                                         ifelse(summary_df$p_fdr < 0.05, "*", "ns")))

summary_df <- summary_df[order(-summary_df$eta_squared), ]

# 添加排名
summary_df$rank <- 1:nrow(summary_df)

cat("\n========== BCH分析结果（按贡献排序） ==========\n")
print(summary_df[, c("rank", "item", "eta_squared", "F_statistic", 
                     "p_value", "p_fdr", "max_cohens_d", "significance")], 
      digits = 4)

# ============================================================
# 提取各类别的BCH加权均值
# ============================================================

means_matrix <- do.call(rbind, lapply(bch_results, function(res) {
  data.frame(
    item = res$item,
    Class1_mean = res$bch_means[1],
    Class2_mean = res$bch_means[2],
    Class3_mean = res$bch_means[3],
    Class1_se = res$bch_se[1],
    Class2_se = res$bch_se[2],
    Class3_se = res$bch_se[3]
  )
}))

# 合并到主结果中
final_results <- merge(summary_df, means_matrix, by = "item")
final_results <- final_results[order(-final_results$eta_squared), ]

# ============================================================
# 保存结果
# ============================================================

write.csv(final_results, "pathway/BCH_results_complete.csv", row.names = FALSE)


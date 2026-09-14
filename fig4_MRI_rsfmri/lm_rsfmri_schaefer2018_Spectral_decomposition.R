library(tidyverse)
library(ggpubr)
library(car)        
library(multcomp)     
library(rstatix)     
library(HH)       
library(broom)
library(emmeans)      
library(data.table)
library(tidyverse)
library(matrixStats)



############################################  schaefer2018  Spectral Decomposition ###############################################################
 
df <- read_csv("pathway/rsfmri/Schaefer2018_416_node.csv")
n_subjects <- nrow(df)

all_cols <- colnames(df)

start_physio <- which(all_cols == "N1_N2")
end_physio <- which(all_cols == "N415_N416")

conn_matrix <- df[start_physio:end_physio]
conn_matrix <- as.matrix(conn_matrix)
if (any(is.na(conn_matrix))) {
  is_na <- is.na(conn_matrix)
  rows <- unique(row(conn_matrix)[is_na])
  df <- df[-rows,]
}

n_analysis <- nrow(df) 


start_covars <- which(all_cols == "instance2")
physio_columns <- all_cols[start_physio:end_physio]

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
         ethnic = factor(ethnic),
         MRI_site_instance2 = factor(MRI_site_instance2),
         EDU = factor(EDU))%>% 
  mutate(subtype = relevel(subtype, ref = "0")) 


design_formula <- ~ 0 + instance2 + sex  + MRI_site_instance2 + EDU + p24441_i2 + ethnic + sleep_scale_age + subtype 

X_design <- model.matrix(design_formula, data = df)

X <- cbind(
  Intercept = 1,  
  X_design)        # model.matrix

conn_matrix <- df[start_physio:end_physio]
conn_matrix <- as.matrix(conn_matrix)

cat(sprintf("get the conn_matrix"))

rm(df_filter)
rm(X_design)
rm(is_na)
gc()


lambda <- 1e-8 * mean(diag(crossprod(X)))

XtX_reg <- crossprod(X) + lambda * diag(ncol(X))
XtX_inv <- solve(XtX_reg)  

XtY <- crossprod(X, conn_matrix)

beta <- XtX_inv %*% XtY  
subtype_cols <- grep("^subtype[123]$", colnames(X), value = FALSE)
beta_subtype <- beta[subtype_cols, ]  

cat(sprintf("beta_subtype right"))


fitted_values <- X %*% beta
residuals <- conn_matrix - fitted_values


n_predictors <- ncol(X)
df_resid <- n_analysis - n_predictors

mse <- colSums(residuals^2, na.rm = TRUE) / df_resid

XtX_inv_diag <- diag(XtX_inv)
se_scale_subtype <- XtX_inv_diag[subtype_cols]

se_subtype <- sqrt(se_scale_subtype %*% t(mse))
t_stats <- beta_subtype / se_subtype


p_values <- 2 * pt(-abs(t_stats), df = df_resid)

rm(fitted_values)
rm(residuals)
rm(matching_rows)
gc()


X_scaled <- scale(conn_matrix, center = TRUE, scale = TRUE)  

# 2. Meff

R_small <- tcrossprod(X_scaled)
eigen_small <- eigen(R_small, symmetric = TRUE, only.values = TRUE)$values
eigen_small <- eigen_small[eigen_small > 1e-10]

m <- nrow(X_scaled)
eig_Rvar <- eigen_small / (m - 1)

f_li_ji <- function(x) {
  if (x >= 1) return(1 + (x - floor(x))) else return(0)
}

Meff <- floor(sum(sapply(abs(eig_Rvar), f_li_ji)))
cat(sprintf("Meff: %.0f\n", Meff))


alpha <- 0.05
alpha_liji <- 1 - (1 - alpha)^(1 / Meff)
p_adj_liji <- 1 - (1 - p_values)^Meff

write.csv(p_adj_liji, 'pathway/adjust_p_Scheafer2018_416_node_liji.csv', row.names = TRUE)


significance_meff <- matrix("ns", nrow = nrow(p_adj_liji), ncol = ncol(p_adj_liji))
significance_meff[p_adj_liji < 0.001] <- "***"
significance_meff[p_adj_liji >= 0.001 & p_adj_liji < 0.01] <- "**"
significance_meff[p_adj_liji >= 0.01 & p_adj_liji < alpha] <- "*"

connection_labels <- all_cols[start_physio:end_physio]


significant_results <- list()
for (i in 1:nrow(p_adj_liji)) {
  comp_name <- rownames(p_values)[i]  
  sig_idx <- which(p_adj_liji[i, ] < alpha)
  if (length(sig_idx) > 0) {

    result_df <- data.frame(
      comparison = comp_name,
      connection_index = sig_idx,
      connection_label = connection_labels[sig_idx],
      df = df_resid,
      beta = beta_subtype[i, sig_idx],
      se = se_subtype[i, sig_idx],
      t_stat = t_stats[i, sig_idx],
      p_raw = p_values[i, sig_idx],
      p_adj_liji = p_adj_liji[i, sig_idx],  
      significance = significance_meff[i, sig_idx],
      effect_size = abs(beta_subtype[i, sig_idx]) / sqrt(mse[sig_idx]),
      mse = mse[sig_idx]
    )
    significant_results[[comp_name]] <- result_df
  }
}


if (length(significant_results) > 0) {
  all_sig_results <- do.call(rbind, significant_results)

  rownames(all_sig_results) <- NULL
  output_file <- sprintf("pathway/rsfmri/significant_spectrum_Scheafer2018_416_node.csv")
  write.csv(all_sig_results, output_file, row.names = FALSE)
} else {
  all_sig_results <- data.frame()
}


output_file <- sprintf("pathway/rsfmri/t_stats_Scheafer2018_416_node.csv")
mat_transposed <- t(t_stats)
write.csv(mat_transposed, output_file, row.names = TRUE)


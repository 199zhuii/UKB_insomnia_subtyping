# Load R packages
library(gprofiler2)
library(data.table)


background_df <- read.csv("pathway/pro_name.csv")  # all Background gene
background_symbols <- background_df$protein_name
# 
background_genes <- gconvert(query = background_symbols, 
                         target = "ENSG",
                         organism = "hsapiens")


# Read Protein results
df_sub1 <- read.csv("pathway/protein_fdr_subtype1.csv")
df_sub1 <- gconvert(query = df_sub1$protein_name, 
                    target = "ENSG",
                    organism = "hsapiens")
df_sub2 <- read.csv("pathway/protein_fdr_subtype2.csv")
df_sub2 <- gconvert(query = df_sub2$protein_name, 
                    target = "ENSG",
                    organism = "hsapiens")
df_sub3 <- read.csv("pathway/protein_fdr_subtype3.csv")
df_sub3 <- gconvert(query = df_sub3$protein_name, 
                    target = "ENSG",
                    organism = "hsapiens")
#version information
get_version_info(organism = "hsapiens")

############################### subtype1 ################################################
sub_1 <- gost(query = df_sub1$target,
              organism = 'hsapiens', ordered_query = FALSE,
              domain_scope = "custom",custom_bg = background_genes$target,
              multi_query = FALSE, exclude_iea = FALSE,
              measure_underrepresentation = FALSE, evcodes = FALSE,
              sources = c('GO:BP','GO:MF','GO:CC','KEGG','REAC'),
              significant = TRUE, user_threshold = 0.05,correction_method = "fdr")

# save results
sub_1_result <- sub_1$result
write.csv(sub_1_result[,1:13],file = 'pathway/protein/fdr_subtype1_pathway.csv',row.names = F)


############################### subtype2 ################################################
sub_2 <- gost(query = df_sub2$target,
              organism = 'hsapiens', ordered_query = FALSE,
              domain_scope = "custom",custom_bg = background_genes$target,
              multi_query = FALSE, exclude_iea = FALSE,
              measure_underrepresentation = FALSE, evcodes = FALSE,
              sources = c('GO:BP','GO:MF','GO:CC','KEGG','REAC'),
              significant = TRUE, user_threshold = 0.05,correction_method = "fdr")

# save results
sub_2_result <- sub_2$result
write.csv(sub_2_result[,1:13],file = 'pathway/protein/fdr_subtype2_pathway.csv',row.names = F)



############################### subtype3 ################################################
sub_3 <- gost(query = df_sub3$target,
              organism = 'hsapiens', ordered_query = FALSE,
              domain_scope = "custom",custom_bg = background_genes$target,
              multi_query = FALSE, exclude_iea = FALSE,
              measure_underrepresentation = FALSE, evcodes = FALSE,
              sources = c('GO:BP','GO:MF','GO:CC','KEGG','REAC'),
              significant = TRUE, user_threshold = 0.05,correction_method = "fdr")

# save results
sub_3_result <- sub_3$result
write.csv(sub_3_result[,1:13],file = 'pathway/protein/fdr_subtype3_pathway.csv',row.names = F)

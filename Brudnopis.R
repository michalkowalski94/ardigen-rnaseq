# Installation of BiocManager to handle Bioconductor ----

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Installation of necessary libraries ----
cran_packages <- c("tidyverse", "pheatmap", "ggrepel", "RColorBrewer", "UpSetR")
bioc_packages <- c("edgeR", "limma", "clusterProfiler",
                   "org.Hs.eg.db", "EDASeq")

# Install CRAN packages if not already installed
missing_cran <- cran_packages[!(cran_packages %in%
                                  installed.packages()[,"Package"])]
if (length(missing_cran) > 0) {
  cat("Instalowanie pakietów z CRAN:", paste(missing_cran, collapse=", "), "\n")
  install.packages(missing_cran)
}

# Install Bioconductor packages if not already installed
missing_bioc <- bioc_packages[!(bioc_packages %in%
                                  installed.packages()[,"Package"])]
if (length(missing_bioc) > 0) {
  cat("Instalowanie pakietów z Bioconductor:",
      paste(missing_bioc, collapse=", "), "\n")
  BiocManager::install(missing_bioc)
}

# Import of necessary libraries ----
library(tidyverse)
library(edgeR)
library(limma)
library(UpSetR)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EDASeq)


# Data reading and pre-processing ----
# This code reads in RNA-Seq data and metadata, checks for missing values, and prepares the data for analysis.

## Import of data ----
counts_data <- read.delim2("./expression_data.tsv", header = TRUE, row.names = 1,                       check.names = FALSE)
# View(counts_data)

## Import of metadata ----
metadata <- read.delim2("./metadata.tsv", header = TRUE, row.names = 1,
                       check.names = FALSE)
# View(metadata)

## Check missing observations and standardize column names ----
# Filter out missing patients, by comparing colnames of count_data with
# rownames in metadata
metadata_set <- rownames(metadata)
counts_data_set <- colnames(counts_data)[6:23]
# Function to remove whitespaces in rownames of metadata
remove_whitespace <- function(x) {
  gsub("\\s+", "", x)
}

#Removal of whitespaces in metadata
metadata[,"tissue"] <- remove_whitespace(metadata[,"tissue"])
metadata[,"smoker"] <- remove_whitespace(metadata[,"smoker"])

# Vector operation on the metadata_set to remove whitespaces
metadata_set <- remove_whitespace(metadata_set)

stopifnot(all(metadata_set %in% counts_data_set))
stopifnot(all(counts_data_set %in% metadata_set))
rownames(metadata) <- metadata_set

# Checking completeness of the records and finding records with missing values
# in counts_data

missing_records <- is.na(counts_data)
# Countring true values
missing_records_count <- colSums(missing_records)
# Displaying the records with missing values
missing_records_display <- counts_data[, missing_records_count > 0]
if (ncol(missing_records_display) > 0) {
  print("Records with missing values:")
  print(missing_records_display)
} else {
  print("No records with missing values found.")
}

# Extracting the raw counts data with preservation of rownames
raw_counts <- counts_data[, metadata_set]
# Displaying the first few rows of the raw counts data
print("Raw counts data:")
print(head(raw_counts))

# Exploratory Data Analysis ----

## Exploration of raw counts data ----

### Construction of DGE List ----
dge <- DGEList(counts = raw_counts, 
               samples = metadata, 
               genes = data.frame(GeneSymbol = rownames(raw_counts)))

# Displaying the DGEList object
print("DGEList object:")
print(dge)

### Filtering lowly expressed genes ----
# Filtering out genes with low counts
keep <- filterByExpr(dge, group = dge$samples$tissue)
dge_filtered <- dge[keep, , keep.lib.sizes=FALSE]
cat("Number of genes before filtering:", nrow(dge), "\n")
cat("Number of genes after filtering:", nrow(dge_filtered), "\n")

### Distribution of library sizes
plot_data <- dge_filtered$samples %>%
  rownames_to_column("Sample")
mean_lib_size_mil <- mean(plot_data$lib.size) / 1e6

ggplot(plot_data, aes(x = reorder(Sample, lib.size), y = lib.size / 1e6)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = mean_lib_size_mil, linetype = "dashed",
             color = "red", size = 1) +
  geom_text(aes(x = -Inf, y = mean_lib_size_mil,
                label = paste("Average =", round(mean_lib_size_mil, 2), "million")),
            hjust = -0.1, vjust = -0.5, color = "red", size = 4) +
  labs(
    title = "Sample library size",
    x = "Sample",
    y = "Library size (million)" #Adjusted for better readability
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Counts distribution within sample, before normalization
logcpm <- cpm(dge_filtered, log=TRUE)
logcpm_long <- as.data.frame(logcpm) %>%
  rownames_to_column("GeneSymbol") %>%
  pivot_longer(
    cols = -GeneSymbol,
    names_to = "Sample",
    values_to = "logCPM"
  ) %>%
  left_join(dge_filtered$samples %>% rownames_to_column("Sample"), by = "Sample")

ggplot(logcpm_long, aes(x = Sample, y = logCPM, fill = tissue)) +
  geom_boxplot(outlier.shape = NA) + 
  labs(
    title = "Distribution of log-CPM values before normalization",
    x = "Sample",
    y = "log-CPM"
  ) +
  # coord_cartesian(ylim = c(0, 20)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

### MA-plot on psuedo-reference ----
sample_index <- 2
sample_name <- colnames(logcpm)[sample_index]
normal_samples_mask <- dge_filtered$samples$tissue == "normal"
pseudo_reference <- rowMeans(logcpm[, normal_samples_mask])

M <- logcpm[, sample_index] - pseudo_reference
A <- (logcpm[, sample_index] + pseudo_reference) / 2
ma_data <- data.frame(A = A, M = M)

ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue", size = 1) +
  labs(
    title = paste("MA-plot for sample:", sample_name),
    subtitle = "Comparison to the pseudo-reference from 'normal' tissue cohort",
    x = "A (Average log-CPM)",
    y = "M (Difference in log-CPM)"
  ) +
  theme_bw(base_size = 8)

### Plotting RLE (Relative Log Expression) before normalization ----
order_df <- dge_filtered$samples %>%
  rownames_to_column("Sample") %>%
  arrange(tissue, Sample)
ordered_counts <- dge_filtered$counts[, order_df$Sample]
ordered_colors <- as.numeric(factor(order_df$tissue))

plotRLE(as.matrix(ordered_counts), 
        outline = FALSE, 
        las = 2, 
        col = ordered_colors,
        main = "RLE before normalization",
        colour_by = "tissue")
legend("topright", 
       legend = levels(factor(order_df$tissue)), 
       fill = 1:nlevels(factor(order_df$tissue)), 
       title = "Tissue")

## Normaalization using TMM ----
# Normalization using TMM (Trimmed Mean of M-values)
dge_norm <- calcNormFactors(dge_filtered, method = "TMM")
temp_design <- model.matrix(~tissue, data = dge_norm$samples)
v_check <- voom(dge_norm, temp_design, plot = TRUE)

### Searching for the source of variability in the normalized data ----

##### Plotting MDS to visualize the variability ----
plotMDS(dge_norm, 
        col = as.numeric(factor(dge_norm$samples$tissue)),
        labels = dge_norm$samples$smoker,
        dim.plot = c(1,2))
title("MDS Plot (counts normalized with TMM")

##### Plotting PCA to visualize the variability ----
logcpm_norm <- cpm(dge_norm, log=TRUE)
pca_res <- prcomp(t(logcpm_norm), scale. = TRUE)
pca_data <- as.data.frame(pca_res$x) %>%
  rownames_to_column("Sample") %>%
  left_join(dge_norm$samples %>% rownames_to_column("Sample"), by = "Sample")

# Computing percentage of variance explained by each principal component
percent_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# Plotting PCA results: PC1 vs PC2
ggplot(pca_data, aes(x = PC1, y = PC2, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA Plot (counts normalized with TMM): PC1 vs PC2",
    x = paste0("PC1: ", round(percent_var[1] * 100), "% variance explained"),
    y = paste0("PC2: ", round(percent_var[2] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

# Plotting the PCA results: PC2 vs PC3
ggplot(pca_data, aes(x = PC2, y = PC3, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA Plot (counts normalized with TMM): PC2 vs PC3",
    x = paste0("PC2: ", round(percent_var[2] * 100), "% variance explained"),
    y = paste0("PC3: ", round(percent_var[3] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

# Plotting the PCA results: PC1 vs PC3
ggplot(pca_data, aes(x = PC1, y = PC3, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA Plot (counts normalized with TMM): PC1 vs PC3",
    x = paste0("PC1: ", round(percent_var[1] * 100), "% variance explained"),
    y = paste0("PC3: ", round(percent_var[3] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

#There is clearly a cluster missing - data is probably being clustered by sex,
#and that should be taken into account before the main analysis.

# Upstream & Differential Expression Analysis ----

## Add missing sex data, based on chromosome Y genes expression ----
gene_info <- counts_data[, "Chr", drop = FALSE]
y_genes_mask <- sapply(strsplit(as.character(gene_info$Chr), ";"), function(x) {
  # Usunięcie pustych wpisów, które mogą powstać z NA
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(FALSE)
  all(x == "chrY")
})

# Extracting the Y chromosome genes from the raw counts data
y_gene_counts <- raw_counts[y_genes_mask, ]
# Creating the log2 transformed expression values for Y chromosome genes
mean_y_expr <- colMeans(log2(y_gene_counts + 1))
# Creating a data frame for plotting the average Y chromosome genes expression
sex_plot_data <- data.frame(
  sample = names(mean_y_expr),
  mean_y_expr = mean_y_expr
)
# Plotting the average Y chromosome genes expression
ggplot(sex_plot_data, aes(x = reorder(sample, mean_y_expr), y = mean_y_expr)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(
    title = "Average Y chromosome genes expression",
    x = "Sample",
    y = "Average expression log2(counts + 1)"
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Settling log2 transformation threshold to 1
threshold <- 1
metadata$sex <- ifelse(mean_y_expr > threshold, "male", "female")

## Construct DGEList ----
dge <- DGEList(counts = raw_counts, 
               samples = metadata, 
               genes = data.frame(GeneSymbol = rownames(raw_counts)))

# Filtering out genes with low counts... again
keep <- filterByExpr(dge, group = dge$samples$tissue)
dge_filtered <- dge[keep, , keep.lib.sizes=FALSE]
cat("Number of genes before filtering:", nrow(dge), "\n")
cat("Number of genes after filtering:", nrow(dge_filtered), "\n")

# Normalization using TMM (Trimmed Mean of M-values)
dge_norm <- calcNormFactors(dge_filtered, method = "TMM")
temp_design <- model.matrix(~tissue, data = dge_norm$samples)
v_check <- voom(dge_norm, temp_design, plot = TRUE)

## Plotting PCA to visualize the variability (including sex) ----
## Preparing separete data for both sexes

male_samples <- dge_norm$samples$sex == "male"
female_samples <- dge_norm$samples$sex == "female"

logcpm_male <- logcpm_norm[, male_samples]
logcpm_female <- logcpm_norm[, female_samples]

### PCA for male cohort ----
pca_res_male <- prcomp(t(logcpm_male), scale. = TRUE)
pca_data_male <- as.data.frame(pca_res_male$x) %>%
  rownames_to_column("Sample") %>%
  left_join(dge_norm$samples %>% rownames_to_column("Sample"), by = "Sample")
percent_var_male <- pca_res_male$sdev^2 / sum(pca_res_male$sdev^2)

ggplot(pca_data_male, aes(x = PC1, y = PC2, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue), type = 't', linetype = "dashed") +
  labs(
    title = "PCA for male cohort: PC1 vs PC2",
    x = paste0("PC1: ", round(percent_var_male[1] * 100), "% variance explained"),
    y = paste0("PC2: ", round(percent_var_male[2] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

ggplot(pca_data_male, aes(x = PC2, y = PC3, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue),type = 't', linetype = "dashed") +
  labs(
    title = "PCA for male cohort: PC2 vs PC3",
    x = paste0("PC2: ", round(percent_var_male[2] * 100), "% variance explained"),
    y = paste0("PC3: ", round(percent_var_male[3] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

ggplot(pca_data_male, aes(x = PC1, y = PC3, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue), type = 't', linetype = "dashed") +
  labs(
    title = "PCA for male cohort: PC1 vs PC3",
    x = paste0("PC1: ", round(percent_var_male[1] * 100), "% variance explained"),
    y = paste0("PC3: ", round(percent_var_male[3] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

### PCA for female cohort ----
pca_res_female <- prcomp(t(logcpm_female), scale. = TRUE)
pca_data_female <- as.data.frame(pca_res_female$x) %>%
  rownames_to_column("Sample") %>%
  left_join(dge_norm$samples %>% rownames_to_column("Sample"), by = "Sample")
percent_var_female <- pca_res_female$sdev^2 / sum(pca_res_female$sdev^2)

ggplot(pca_data_female, aes(x = PC1, y = PC2, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue), type = 't', linetype = "dashed") +
  labs(
    title = "PCA for female cohort: PC1 vs PC2",
    x = paste0("PC1: ", round(percent_var_male[1] * 100), "% variance explained"),
    y = paste0("PC2: ", round(percent_var_male[2] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

ggplot(pca_data_female, aes(x = PC2, y = PC3, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue),type = 't', linetype = "dashed") +
  labs(
    title = "PCA for female cohort: PC2 vs PC3",
    x = paste0("PC2: ", round(percent_var_male[2] * 100), "% variance explained"),
    y = paste0("PC3: ", round(percent_var_male[3] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

ggplot(pca_data_female, aes(x = PC1, y = PC3, color = tissue, shape = smoker)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue), type = 't', linetype = "dashed") +
  labs(
    title = "PCA for female cohort: PC1 vs PC3",
    x = paste0("PC1: ", round(percent_var_male[1] * 100), "% variance explained"),
    y = paste0("PC3: ", round(percent_var_male[3] * 100), "% variance explained")
  ) +
  theme_bw(base_size = 8)

## The female cohort demonstrates huge clustering of tumor samples against
## the normal samples. It clearly shows that some genes are different from the
## reference tissue, but as the normal tissue clusters quite tightly, the tumor
## tissue seems to be quite dispersed on PC1 vs PC2 plot. The PC2 plot vs PC3 plot
## shows same history. Sex was definetely the noisy factor in here, thus separate
## analyses will be performed on male and female cohorts

## Differential expression analysis & Downstream ----

### Task 1 & 2 ----

# Setting levels for factors in the DGEList object
dge_norm$samples$tissue <- factor(dge_norm$samples$tissue,
                                  levels = c("normal", "tumor"))
dge_norm$samples$smoker <- factor(dge_norm$samples$smoker,
                                  levels = c("never", "former", "current"))

#### Task 1: Difference between tumor tissue and normal tissue ----

tissue_design <- model.matrix(~0+tissue, data = dge_norm$samples)
colnames(tissue_design) <- gsub("tissuetumor", "tumor", colnames(tissue_design))
colnames(tissue_design) <- gsub("tissuenormal", "normal", colnames(tissue_design))

tissue_design

v <- voom(dge_norm, tissue_design, plot = TRUE)
fit <- lmFit(v, tissue_design)

cont.matrix <- makeContrasts(
  tumor_vs_normal = tumor - normal,
  levels = colnames(tissue_design)
)

cont.matrix

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2, p.value = 0.05, lfc = 1))

#Extracting results for the tumor vs normal comparison
tumor_vs_normal_results <- topTable(fit2, coef = "tumor_vs_normal",
                                    number = Inf, sort.by = "none",
                                    adjust.method = "BH",
                                    confint = T) %>%
  arrange(adj.P.Val) %>%
  mutate(Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Displaying the top results
View(tumor_vs_normal_results)

#### Task 2: Difference between smoking status ----
# Creating a new design matrix for smoking status
smoking_design <- model.matrix(~0+smoker, data = dge_norm$samples)
colnames(smoking_design) <- gsub("smokercurrent", "current", colnames(smoking_design))
colnames(smoking_design) <- gsub("smokerformer", "former", colnames(smoking_design))
colnames(smoking_design) <- gsub("smokernever", "never", colnames(smoking_design))

smoking_design

v_smoking <- voom(dge_norm, smoking_design, plot = TRUE)
fit_smoking <- lmFit(v_smoking, smoking_design)
cont.matrix_smoking <- makeContrasts(
  current_vs_never = current - never,
  former_vs_never = former - never,
  current_vs_former = current - former,
  levels = colnames(smoking_design)
)

cont.matrix_smoking

fit2_smoking <- contrasts.fit(fit_smoking, cont.matrix_smoking)
fit2_smoking <- eBayes(fit2_smoking)
summary(decideTests(fit2_smoking, p.value = 0.05, lfc = 1))

#Extracting results for the current vs never comparison
current_vs_never_results <- topTable(fit2_smoking, coef = "current_vs_never",
                                     number = Inf, sort.by = "none",
                                     adjust.method = "BH", confint = T) %>%
  arrange(adj.P.Val) %>%
  mutate(Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))
# Displaying the top results
View(current_vs_never_results)

former_vs_never_results <- topTable(fit2_smoking, coef = "former_vs_never",
                                    number = Inf, sort.by = "none",
                                    adjust.method = "BH", confint = T) %>%
  arrange(adj.P.Val) %>%
  mutate(Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))
# Displaying the top results
View(former_vs_never_results)

# Extracting results for the current vs former comparison
current_vs_former_results <- topTable(fit2_smoking, coef = "current_vs_former",
                                      number = Inf, sort.by = "none",
                                      adjust.method = "BH", confint = T) %>%
  arrange(adj.P.Val) %>%
  mutate(Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))
# Displaying the top results
View(current_vs_former_results)


#### Task 3 finding interactions between smokers and tumor ----

# Adding smoking status to the DGEList object
dge_norm$samples$smoking_status <- ifelse(dge_norm$samples$smoker == 'never', 'no', 'yes')
dge_norm$samples$smoking_status <- factor(dge_norm$samples$smoking_status, levels = c("no", "yes"))

# Creating a new design matrix including smoking status
design_interaction <- model.matrix(~0 + tissue * smoking_status, data = dge_norm$samples)
colnames(design_interaction) <- gsub("tissuetumor:smoking_statusyes", "smokertumor", colnames(design_interaction))
print(design_interaction)

# Performing voom transformation
v_int <- voom(dge_norm, design_interaction, plot = TRUE)
fit_int <- lmFit(v_int, design_interaction)

cont.matrix_int <- makeContrasts(
  smoking_tumor_vs_normal = smokertumor - tissuenormal,
  levels = colnames(design_interaction)
)

cont.matrix_int
fit2_int <- contrasts.fit(fit_int, cont.matrix_int)

fit2_int <- eBayes(fit2_int)

summary(decideTests(fit2_int, p.value = 0.05, lfc = 1))

# Extracting results for the smoking_tumor_vs_normal comparison
smoking_tumor_vs_normal_results <- topTable(fit2_int, coef = "smoking_tumor_vs_normal",
                                            number = Inf, sort.by = "none",
                                            adjust.method = "BH", confint = T) %>%
  arrange(adj.P.Val) %>%
  mutate(Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))
# Displaying the top results
View(smoking_tumor_vs_normal_results)

#### MA plot for tumor vs normal comparison ----
tumor_vs_normal_results$ggplotSignificant <- ifelse(
  abs(tumor_vs_normal_results$logFC) >= 2 & 
    tumor_vs_normal_results$adj.P.Val < 0.05,
  "Yes",
  "No"
)

tumor_normal_ggplot <- ggplot(tumor_vs_normal_results, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = ggplotSignificant), alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  labs(
    title = "MA Plot: Tumor vs Normal",
    x = "Average Expression",
    y = "Log Fold Change"
  ) +
  theme_bw(base_size = 8)

tumor_normal_ggplot

#### MA plot for smoking_tumor vs normal comparison ----

smoking_tumor_vs_normal_results$ggplotSignificant <- ifelse(
  abs(smoking_tumor_vs_normal_results$logFC) >= 4 & 
    smoking_tumor_vs_normal_results$adj.P.Val < 0.05,
  "Yes",
  "No"
)

smoking_tumor_normal_ggplot <- ggplot(smoking_tumor_vs_normal_results, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = ggplotSignificant), alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  labs(
    title = "MA Plot: Smoking Tumor vs Normal",
    x = "Average Expression",
    y = "Log Fold Change"
  ) +
  theme_bw(base_size = 8)

smoking_tumor_normal_ggplot

#### MA plots for interactions between types of smokers ----

##### MA plot for current vs never comparison ----

current_vs_never_results$ggplotSignificant <- ifelse(
  abs(current_vs_never_results$logFC) >= 2 & 
    current_vs_never_results$P.Value < 0.05,
  "Yes",
  "No"
)

current_never_ggplot <- ggplot(current_vs_never_results, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = ggplotSignificant), alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  labs(
    title = "MA Plot: Current vs Never",
    x = "Average Expression",
    y = "Log Fold Change"
  ) +
  theme_bw(base_size = 8)

current_never_ggplot

##### MA plot for former vs never comparison ----

former_vs_never_results$ggplotSignificant <- ifelse(
  abs(former_vs_never_results$logFC) >= 1 & 
    former_vs_never_results$P.Value < 0.05,
  "Yes",
  "No"
)

former_never_ggplot <- ggplot(former_vs_never_results, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = ggplotSignificant), alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  labs(
    title = "MA Plot: Former vs Never",
    x = "Average Expression",
    y = "Log Fold Change"
  ) +
  theme_bw(base_size = 8)

former_never_ggplot

##### MA plot for current vs former comparison ----

current_vs_former_results$ggplotSignificant <- ifelse(
  abs(current_vs_former_results$logFC) >= 2 & 
    current_vs_former_results$P.Value < 0.05,
  "Yes",
  "No"
)

current_former_ggplot <- ggplot(current_vs_former_results, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = ggplotSignificant), alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  labs(
    title = "MA Plot: Current vs Former",
    x = "Average Expression",
    y = "Log Fold Change"
  ) +
  theme_bw(base_size = 8)

current_former_ggplot

#### Definition of function for enriched GO

run_go_enrichment <- function(de_results, p_value_col = "adj.P.Val", p_cutoff = 0.05, fc_cutoff = 1) {
  # Filter for significant genes based on specified cutoffs
  sig_genes <- de_results %>%
    filter(.data[[p_value_col]] < p_cutoff, abs(logFC) > fc_cutoff) %>%
    pull(GeneSymbol)
  
  if (length(sig_genes) == 0) {
    cat("No significant genes found for the given criteria. Skipping enrichment analysis.\n")
    return(NULL)
  }
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = sig_genes,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Perform GO enrichment analysis for Biological Process
  ego <- enrichGO(gene = na.omit(entrez_ids),
                  OrgDb = org.Hs.eg.db,
                  keyType = 'ENTREZID',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.02)
  
  return(ego)
}

#### Pathway enrichment analysis using clusterProfiler for tumor vs normal----

go_tumor_vs_normal <- run_go_enrichment(tumor_vs_normal_results,
                                        p_value_col = "adj.P.Val",
                                        p_cutoff = 0.05,
                                        fc_cutoff = 2)

# Visualize the results if any were found
if (!is.null(go_tumor_vs_normal) && nrow(as.data.frame(go_tumor_vs_normal)) > 0) {
  print(dotplot(go_tumor_vs_normal, showCategory = 10) + 
          labs(title = "GO Enrichment: Tumor vs Normal"))
}

#### Pathway enrichment analysis using clusterProfiler for smoking_tumor vs normal ----

go_smoking_tumor_vs_normal <- run_go_enrichment(smoking_tumor_vs_normal_results,
                                                  p_value_col = "adj.P.Val",
                                                  p_cutoff = 0.05,
                                                  fc_cutoff = 6)
# Visualize the results if any were found
if (!is.null(go_smoking_tumor_vs_normal) && nrow(as.data.frame(go_smoking_tumor_vs_normal)) > 0) {
  print(dotplot(go_smoking_tumor_vs_normal, showCategory = 30) + 
          labs(title = "GO Enrichment: Smoking Tumor vs Normal"))
}

#### Pathway enrichment analysis using clusterProfiler for current vs never comparison ----

go_current_vs_never <- run_go_enrichment(current_vs_never_results,
                                           p_value_col = "P.Value",
                                           p_cutoff = 0.05,
                                           fc_cutoff = 2)
# Visualize the results if any were found
if (!is.null(go_current_vs_never) && nrow(as.data.frame(go_current_vs_never)) > 0) {
  print(dotplot(go_current_vs_never, showCategory = 10) + 
          labs(title = "GO Enrichment: Current vs Never"))
}

# Notice that this result is not statistically significant after p-value adjustment
# hence they can only be treated as exploratory results.

#### Pathway enrichment analysis using clusterProfiler for former vs never comparison ----

go_former_vs_never <- run_go_enrichment(former_vs_never_results,
                                          p_value_col = "P.Value",
                                          p_cutoff = 0.05,
                                          fc_cutoff = 1)
# Visualize the results if any were found
if (!is.null(go_former_vs_never) && nrow(as.data.frame(go_former_vs_never)) > 0) {
  print(dotplot(go_former_vs_never, showCategory = 10) + 
          labs(title = "GO Enrichment: Former vs Never"))
}

# Notice that this result is not statistically significant after p-value adjustment
# hence they can only be treated as exploratory results.

#### Pathway enrichment analysis using clusterProfiler for current vs former comparison ----

go_current_vs_former <- run_go_enrichment(current_vs_former_results,
                                             p_value_col = "P.Value",
                                             p_cutoff = 0.05,
                                             fc_cutoff = 2)
# Visualize the results if any were found
if (!is.null(go_current_vs_former) && nrow(as.data.frame(go_current_vs_former)) > 0) {
  print(dotplot(go_current_vs_former, showCategory = 10) + 
          labs(title = "GO Enrichment: Current vs Former"))
}
# Notice that this result is not statistically significant after p-value adjustment
# hence they can only be treated as exploratory results.

#### Heatmaps ----
##### Function for generation of heatmaps ----

# Function takes as an argument the DGEList object and GO enrichment results
# and generates heatmaps for the differentially expressed genes in each significant category
plot_go_heatmaps <- function(voom_obj, go_results, de_results, annotation_col, 
                             n_top = 5, p_value_col = "adj.P.Val", p_cutoff = 0.05, fc_cutoff = 1) {
  if (is.null(go_results) || nrow(as.data.frame(go_results)) == 0) {
    cat("No GO results to plot.\n")
    return(invisible(NULL))
  }
  
  # Get top n GO terms
  top_go_terms <- head(as.data.frame(go_results), n_top)
  
  # Get the expression matrix
  expr_matrix <- voom_obj$E
  
  # Get the list of all significant genes from the DE results
  all_sig_genes <- de_results %>%
    filter(.data[[p_value_col]] < p_cutoff, abs(logFC) > fc_cutoff) %>%
    pull(GeneSymbol)
  
  # Loop through top GO terms and create a heatmap for each
  for (i in 1:nrow(top_go_terms)) {
    term_description <- top_go_terms$Description[i]
    gene_ids <- top_go_terms$geneID[i]
    
    # Convert gene IDs string to a vector of Entrez IDs
    genes_in_term_entrez <- str_split(gene_ids, "/")[[1]]
    
    # Convert Entrez IDs back to Gene Symbols
    symbols_in_term <- mapIds(org.Hs.eg.db,
                              keys = genes_in_term_entrez,
                              column = "SYMBOL",
                              keytype = "ENTREZID",
                              multiVals = "first")
    
    symbols_in_term <- na.omit(symbols_in_term)
    
    # Find which of these genes are BOTH in the GO term AND significant
    genes_to_plot <- intersect(symbols_in_term, all_sig_genes)
    
    # Final check to ensure these genes are in our expression matrix
    genes_to_plot <- intersect(genes_to_plot, rownames(expr_matrix))
    
    if (length(genes_to_plot) < 2) {
      cat("Skipping heatmap for '", term_description, "' (fewer than 2 significant genes found in data).\n", sep="")
      next
    }
    
    # Subset the expression matrix
    mat_subset <- expr_matrix[genes_to_plot, ]
    
    # Create the heatmap
    pheatmap(mat_subset,
             main = paste("Heatmap for:", term_description),
             annotation_col = voom_obj$targets[, annotation_col, drop = FALSE],
             scale = "row", # Scale genes to have mean 0 and sd 1
             show_colnames = FALSE,
             fontsize_row = 8)
  }
}

##### Heatmaps for tumor vs normal comparison ----
plot_go_heatmaps(voom_obj = v, 
                 go_results = go_tumor_vs_normal,
                 de_results = tumor_vs_normal_results,
                 annotation_col = "tissue", 
                 n_top = 10,
                 fc_cutoff = 2,
                 p_value_col = "adj.P.Val",
                 p_cutoff = 0.05)
##### Heatmaps for smoking_tumor vs normal comparison ----
v_int$targets$smoker_tumor <- ifelse(v_int$targets$tissue == "tumor" & 
                                      v_int$targets$smoking_status == "yes", 
                                      "yes", "no")
plot_go_heatmaps(voom_obj = v_int,
                 go_results = go_smoking_tumor_vs_normal,
                 de_results = smoking_tumor_vs_normal_results,
                 annotation_col = "smoker_tumor", 
                 n_top = 10,
                 fc_cutoff = 8,
                 p_value_col = "adj.P.Val",
                 p_cutoff = 0.05)

##### Heatmaps for current vs never comparison ----
plot_go_heatmaps(voom_obj = v_smoking,
                 go_results = go_current_vs_never,
                 de_results = current_vs_never_results,
                 annotation_col = "smoker", 
                 n_top = 10,
                 fc_cutoff = 2,
                 p_value_col = "P.Value",
                 p_cutoff = 0.05)

##### Heatmaps for former vs never comparison ----
plot_go_heatmaps(voom_obj = v_smoking,
                 go_results = go_former_vs_never,
                 de_results = former_vs_never_results,
                 annotation_col = "smoker", 
                 n_top = 10,
                 fc_cutoff = 1,
                 p_value_col = "P.Value",
                 p_cutoff = 0.05)

##### Heatmaps for current vs former comparison ----
plot_go_heatmaps(voom_obj = v_smoking,
                 go_results = go_current_vs_former,
                 de_results = current_vs_former_results,
                 annotation_col = "smoker", 
                 n_top = 10,
                 fc_cutoff = 2,
                 p_value_col = "P.Value",
                 p_cutoff = 0.05)

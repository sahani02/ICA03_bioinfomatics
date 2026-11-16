library(GEOquery)
library(limma)
library(ggplot2)
library(tidyverse)
library(pheatmap)


# Download GEO dataset
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/matrix/GSE42568_series_matrix.txt.gz",
  destfile = "GSE42568_series_matrix.txt.gz",
  mode = "wb",
  method = "libcurl"
)

# Load dataset
gse <- getGEO(
  filename = "GSE42568_series_matrix.txt.gz"
)

# Extract expression data
expression_data <- exprs(gse)

# Clean and prepare data
# Extract sample information
metadata <- pData(gse)

# Extract sample titles
sample_titles <- metadata$title

# Identify sample groups
group <- ifelse(
  grepl("normal", sample_titles, ignore.case = TRUE),
  "Normal",
  "Cancer"
)

# Convert group into factor
group <- factor(group)

# Display sample counts
table(group)

# View sample titles
head(metadata$title)

# Check group counts
table(group)

# Check missing values
sum(is.na(expression_data))

#Check duplicated sample names
anyDuplicated(colnames(expression_data))

#Exploratory Data Analysis (EDA)
#Sample distribution
group_data <- data.frame(
  Group = group
)

ggplot(group_data, aes(x = Group, fill = Group)) +
  geom_bar() +
  labs(
    title = "Distribution of Samples",
    x = "Sample Type",
    y = "Number of Samples"
  ) +
  theme_minimal()

#Gene expression distribution
boxplot(
  expression_data[, 1:10],
  las = 2,
  main = "Gene Expression Distribution",
  ylab = "Expression Level"
)


# Check dimensions
dim(expression_data)

# Check missing values
sum(is.na(expression_data))

# Check genes with zero variance
gene_variance <- apply(
  expression_data,
  1,
  var,
  na.rm = TRUE
)

sum(gene_variance == 0, na.rm = TRUE)

# Remove genes with zero variance

nonzero_genes <- gene_variance > 0 &
  !is.na(gene_variance)

expression_filtered <- expression_data[
  nonzero_genes,
]

# Perform PCA

pca <- prcomp(
  t(expression_filtered),
  scale. = TRUE
)

# Create PCA dataframe

pca_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group
)

# Create PCA plot

ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Breast Cancer Samples",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()

#Main analysis – Differential Gene Expression
# Create design matrix
design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

print(design)
# Fit linear model
fit_model <- lmFit(expression_data, design)

# Compare cancer and normal samples
contrast <- makeContrasts(
  Cancer - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit_model, contrast)

fit2 <- eBayes(fit2)

results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

head(results)

class(fit_model)

head(results)

dim(results)
# Extract results
results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)

head(results)

significant_genes <- results %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

# Display significant genes
head(significant_genes)

# Count significant genes
nrow(significant_genes)

#save significant genes
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  significant_genes,
  "results/tables/significant_genes.csv",
  row.names = TRUE
)

#Create a volcano plot

results$Significance <- ifelse(
  results$adj.P.Val < 0.05 &
    abs(results$logFC) > 1,
  "Significant",
  "Not Significant"
)

ggplot(
  results,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Significance
  )
) +
  geom_point(alpha = 0.6) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme_minimal()

#Perform clustering and create a heatmap
gene_variance <- apply(
  expression_data,
  1,
  var
)

top_genes <- names(
  sort(gene_variance, decreasing = TRUE)
)[1:30]

heatmap_data <- expression_data[
  top_genes,
]

#create heatmap
annotation <- data.frame(
  Group = group
)

rownames(annotation) <- colnames(heatmap_data)

pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation,
  show_rownames = FALSE,
  main = "Breast Cancer Gene Expression Heatmap"
)
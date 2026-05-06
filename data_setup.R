###########################################################
# Step 0. Install and load packages
###########################################################
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager", repos = "https://cloud.r-project.org")
# }
# 
# Install required Bioconductor packages
# BiocManager::install(c("GEOquery", "DESeq2"), force = TRUE)

# Load libraries
library(GEOquery)
library(DESeq2)

###########################################################
# Step 1. Download GSM-level supplementary files
###########################################################
gsm_list <- c(
  "GSM8346402","GSM8346403","GSM8346404","GSM8346405",  # control
  "GSM8346417","GSM8346418","GSM8346419","GSM8346420"   # drought
)

# Download files for all GSMs
lapply(gsm_list, getGEOSuppFiles)

###########################################################
# Step 2. Locate abundance.tsv.gz files
###########################################################
all_files <- list.files(".", pattern = "abundance\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)

# Keep only the 8 WT root files
wt_root_files <- all_files[sapply(all_files, function(f) {
  any(sapply(gsm_list, grepl, f))
})]

length(wt_root_files)  # should be 8
wt_root_files

###########################################################
# Step 3. Extract est_counts and build count matrix
###########################################################
counts_list <- lapply(wt_root_files, function(f) {
  df <- read.table(gzfile(f), header = TRUE, row.names = 1)
  est_counts <- df$est_counts
  names(est_counts) <- rownames(df)  # keep gene IDs!
  return(est_counts)
})

df <- read.table(gzfile("./GSM8346402/GSM8346402_S01.abundance.tsv.gz"), header = TRUE, row.names = 1)

# Combine into a single matrix
counts <- do.call(cbind, counts_list)

# Remove version suffix (e.g., ".1", ".2", etc.)
gene_ids <- sub("\\.\\d+$", "", rownames(counts))

# Aggregate counts across transcripts for each gene
counts_gene <- rowsum(counts, group = gene_ids)

# Check result
head(rownames(counts_gene))
dim(counts_gene)

# Assign column names
colnames(counts_gene) <- c("Control1","Control2","Control3","Control4",
                           "Drought1","Drought2","Drought3","Drought4")


###########################################################
# Step 4. Create metadata
###########################################################
condition <- factor(c(rep("control", 4), rep("drought", 4)))
metadata <- data.frame(condition = condition)
rownames(metadata) <- colnames(counts_gene)

all(rownames(metadata) == colnames(counts_gene))  # should be TRUE

###########################################################
# Step 5. Round
###########################################################
# Round counts to integers
counts_int <- round(counts_gene)

# Confirm integers
all(counts_int == floor(counts_int))  # should be TRUE

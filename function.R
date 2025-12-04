#' Data Pre-processing
#'
#' @description
#' From TCGA clinical, count, and survival TSV files, this function:
#' 1) aligns common samples,
#' 2) maps Ensembl IDs to gene symbols (via biomaRt),
#' 3) cleans clinical info + creates survival variables,
#' 4) selects top variable genes and Z-score normalizes expression,
#' 5) merges clinical + expression into one dataset,
#' 6) writes the merged dataset to a CSV file.
#'
#' @param clinical_file Character. Path to clinical TSV file.
#' @param count_file Character. Path to count TSV file.
#' @param survival_file Character. Path to survival TSV file.
#' @param top_n_genes Integer. Number of most variable genes to keep. Default 2000.
#'
#' @return A list with:
#'   \item{merged_data}{Final merged data frame (also written to CSV).}
#'   \item{clinical_clean}{Cleaned clinical data.}
#'   \item{expr_z}{Z-score normalized expression matrix (genes x samples).}
#'   \item{common_samples}{Vector of sample IDs used.}
#'
#' @export

# Load Library
library(readr)
library(dplyr)
library(biomaRt)

build_tcga <- function(clinical_file,
                       count_file,
                       survival_file,
                       top_n_genes = 2000) {
  
  # Read data 
  luad_clinical <- read_tsv(clinical_file, show_col_types = FALSE)
  luad_count <- read_tsv(count_file, show_col_types = FALSE)
  luad_survival <- read_tsv(survival_file, show_col_types = FALSE)
  
  # Map sample IDs across 3 datasets
  clinical_sample <- luad_clinical$sample
  count_sample <- colnames(luad_count)[-1]
  survival_sample <- luad_survival$sample
  
  common_samples <- Reduce(
    intersect,
    list(clinical_sample, count_sample, survival_sample)
  )
  cat("Number of common samples across 3 datasets:", length(common_samples), "\n")
  
  clinical_idx <- match(common_samples, clinical_sample)
  luad_clinical <- luad_clinical[clinical_idx, ]
  
  survival_idx <- match(common_samples, survival_sample)
  luad_survival <- luad_survival[survival_idx, ]
  
  count_idx <- match(common_samples, count_sample)
  luad_count <- luad_count[, c(1, count_idx + 1)]  # +1 since the first column is Ensembl_ID
  
  # Ensembl to gene symbol using biomaRt
  # expression matrix (genes x samples)
  expr_mat <- as.matrix(luad_count[, -1])
  # strip version from Ensembl IDs
  rownames(expr_mat) <- sub("\\..*$", "", luad_count$Ensembl_ID)
  
  # connect to Ensembl
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # get annotation for the Ensembl IDs in the matrix
  annotation <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = unique(rownames(expr_mat)),
    mart = ensembl
  )
  
  # map Ensembl to gene symbol
  gene_symbols <- annotation$external_gene_name[
    match(rownames(expr_mat), annotation$ensembl_gene_id)
  ]
  
  # report Ensembl IDs not in annotation
  missing_ids <- setdiff(rownames(expr_mat), annotation$ensembl_gene_id)
  cat("Number of Ensembl IDs not found in biomaRt:", length(missing_ids), "\n")
  
  # keep genes with non-empty symbols
  keep_nonempty <- !is.na(gene_symbols) & gene_symbols != ""
  expr_mat <- expr_mat[keep_nonempty, , drop = FALSE]
  gene_symbols  <- gene_symbols[keep_nonempty]
  
  # remove duplicated gene symbols
  dup_indices <- duplicated(gene_symbols)
  dup_genes <- unique(gene_symbols[dup_indices])
  cat("Number of duplicated gene symbols removed:", length(dup_genes), "\n")
  
  expr_mat <- expr_mat[!dup_indices, , drop = FALSE]
  gene_symbols <- gene_symbols[!dup_indices]
  
  # final rownames = gene symbols
  rownames(expr_mat) <- gene_symbols
  cat("Number of genes after mapping and filtering:", nrow(expr_mat), "\n")
  
  # clean clinical and add survival
  luad_clinical$OS <- luad_survival$OS
  luad_clinical$OS.time <- luad_survival$OS.time
  
  luad_clinical_clean <- luad_clinical %>%
    dplyr::mutate(
      sample_full = sample,
      patient_id  = substr(sample, 1, 12),
      race = trimws(tolower(race.demographic)),
      race = ifelse(
        is.na(race) | race == "" | race == "not reported",
        NA_character_,
        race
      )
    ) %>%
    dplyr::transmute(
      sample       = sample_full,
      patient_id   = patient_id,
      disease_type = disease_type,
      primary_site = primary_site,
      race         = race,
      gender       = gender.demographic,
      age          = as.numeric(age_at_index.demographic),
      sample_type_code = substr(sample_full, 14, 15),
      tissue_status = dplyr::case_when(
        sample_type_code %in% c("01", "02", "03", "06") ~ "Primary_tumor",
        sample_type_code %in% c("11") ~ "Normal",
        TRUE ~ "Other"
      ),
      # survival outcome
      event = OS,
      time  = OS.time
    ) %>%
    dplyr::filter(!is.na(race))  # drop samples with missing race
  
  # Normalization, select top variable genes and Z-score
  expr_log <- expr_mat
  gene_var <- apply(expr_log, 1, var, na.rm = TRUE)
  
  top_n   <- min(top_n_genes, nrow(expr_log))
  top_idx <- order(gene_var, decreasing = TRUE)[1:top_n]
  expr_log_filt <- expr_log[top_idx, , drop = FALSE]
  
  expr_z <- t(scale(t(expr_log_filt)))  # genes x samples
  
  # Merge clinical and expression data
  expr_df <- as.data.frame(t(expr_z))   # samples x genes
  expr_df$sample <- rownames(expr_df)
  
  merged_data <- merge(
    luad_clinical_clean,
    expr_df,
    by  = "sample",
    all = FALSE
  )
  
  clinical_cols <- c(
    "sample", "patient_id", "tissue_status",
    "event", "time",
    "age", "gender", "race",
    "primary_site", "disease_type",
    "sample_type_code"
  )
  
  gene_cols <- setdiff(colnames(merged_data), clinical_cols)
  merged_data <- merged_data[, c(clinical_cols, gene_cols)]
  
  merged_data$event <- factor(merged_data$event, levels = c(0, 1), labels = c("Alive", "Dead"))
  merged_data$race <- factor(merged_data$race, levels = c("white", "black or african american", "asian", "american indian or alaska native"),
                             labels = c("White", "Black or African American", "Asian", "American Indian or Alaska Native"))
  merged_data$gender <- factor(merged_data$gender, levels = c("male", "female"), labels = c("Male", "Female"))
  merged_data$tissue_status <- factor(merged_data$tissue_status, levels = c("Primary_tumor", "Normal"),
                                      labels = c("Primary Tumor", "Normal"))
  list(
    merged_data    = merged_data,
    clinical_clean = luad_clinical_clean,
    expr_z         = expr_z,
    common_samples = common_samples
  )
}
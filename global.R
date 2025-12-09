suppressPackageStartupMessages({
  library(shiny)
  library(shinycssloaders)
  library(tidyverse)
  library(ggplot2)
  library(plotly)
  library(pheatmap)
  library(survival)
  library(randomForestSRC)
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(pROC)
})

app_config <- list(
  data_dir   = "top2000",
  cache_dir  = "cache",
  top_n_genes = 2000,
  use_cache   = TRUE
)

cancer_codes <- list(
  "Thyroid" = "thca",
  "Liver" = "lihc",
  "Kidney" = "kirc",
  "Breast" = "brca",
  "Colon" = "coad",
  "Bronchus and Lung" = "luad"
)

clinical_fixed_cols <- c(
  "sample", "patient_id", "tissue_status", "event", "time",
  "age", "gender", "race", "primary_site", "disease_type",
  "sample_type_code"
)

cache_path <- function(ct, type) {
  code <- cancer_codes[[ct]]
  
  dir <- app_config$cache_dir
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  file.path(app_config$cache_dir, paste0(type, "_", code, ".rds"))
}

cache_read <- function(path) {
  if (app_config$use_cache && file.exists(path)) {
    readRDS(path)
  } else {
    NULL
  }
}

cache_write <- function(object, path) {
  if (app_config$use_cache) {
    saveRDS(object, path)
  }
}

load_expr_log <- function(ct) {
  cp <- cache_path(ct, "expr_log")
  cached <- cache_read(cp)
  if (!is.null(cached)) return(cached)
  
  message("Loading expression matrix for: ", ct)
  code <- cancer_codes[[ct]]
  
  file <- file.path(app_config$data_dir,
                    paste0("TCGA_", toupper(code), "_merged_2000genes.csv"))
  
  df <- data.table::fread(file)
  data.table::setDF(df)
  
  gene_cols <- setdiff(
    names(df)[sapply(df, is.numeric)],
    c("age", "time", "sample_type_code")
  )
  
  expr_mat <- as.matrix(t(df[, gene_cols, drop = FALSE]))
  rownames(expr_mat) <- gene_cols
  colnames(expr_mat) <- df$sample
  
  expr_log <- log2(expr_mat + 1)
  
  cache_write(expr_log, cp)
  return(expr_log)
}

get_expr_log <- function(ct) load_expr_log(ct)


load_merged_data <- function(ct) {
  cp <- cache_path(ct, "merged")
  cached <- cache_read(cp)
  if (!is.null(cached)) return(cached)
  
  message("Loading merged clinical data for: ", ct)
  code <- cancer_codes[[ct]]
  
  file <- file.path(app_config$data_dir,
                    paste0("TCGA_", toupper(code), "_merged_2000genes.csv"))
  
  df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))
  
  clinical_cols_present <- intersect(clinical_fixed_cols, names(df))
  
  gene_cols <- setdiff(
    names(df)[sapply(df, is.numeric)],
    c("age", "time", "sample_type_code")
  )
  
  vars <- apply(df[, gene_cols, drop = FALSE], 2, var, na.rm = TRUE)
  top_n <- min(app_config$top_n_genes, length(vars))
  top_genes <- names(sort(vars, decreasing = TRUE))[1:top_n]
  
  df_final <- cbind(
    df[, clinical_cols_present, drop = FALSE],
    df[, top_genes, drop = FALSE]
  )
  
  cache_write(df_final, cp)
  return(df_final)
}

get_merged_data <- function(ct) load_merged_data(ct)


load_surv_df <- function(ct) {
  cp <- cache_path(ct, "surv_df")
  cached <- cache_read(cp)
  if (!is.null(cached)) return(cached)
  
  message("Generating survival-cleaned data for: ", ct)
  df <- load_merged_data(ct)
  
  df2 <- df %>%
    filter(tissue_status == "Primary Tumor") %>%
    mutate(
      event = case_when(
        event == "Dead"  ~ 1,
        event == "Alive" ~ 0,
        TRUE ~ NA_real_
      ),
      time = as.numeric(time),
      age = as.numeric(age),
      gender = as.factor(gender)
    ) %>%
    drop_na(time, event, age, gender)
  
  cache_write(df2, cp)
  return(df2)
}

get_surv_df <- function(ct) load_surv_df(ct)


get_gene_list <- function(ct) {
  expr <- get_expr_log(ct)
  gene_var <- apply(expr, 1, var, na.rm = TRUE)
  names(sort(gene_var, decreasing = TRUE))[1:10]
}

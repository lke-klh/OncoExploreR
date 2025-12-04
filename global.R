suppressPackageStartupMessages({
  library(shiny)
  library(shinycssloaders)
  library(tidyverse)
  library(ggplot2)
  library(plotly)
  library(pheatmap)
  library(biomaRt)
})

source("function.R")

app_config <- list(
  data_dir = "data",
  cache_dir = "cache",
  top_n_genes = 2000,
  use_cache = TRUE
)

cancer_codes <- list(
  "Thyroid" = "thca",
  "Liver" = "lihc",
  "Kidney" = "kirc",
  "Breast" = "brca",
  "Colon" = "coad",
  "Bronchus and Lung" = "luad"
)

.cancer_env <- new.env(parent = emptyenv())

load_cancer_data <- function(
    cancer_type,
    data_dir = app_config$data_dir,
    cache_dir = app_config$cache_dir,
    top_n_genes = app_config$top_n_genes,
    use_cache = app_config$use_cache
) {
  
  # warnings for mismatch
  cancer_code <- cancer_codes[[cancer_type]]
  if (is.null(cancer_code)) {
    warning("Unknown cancer type: ", cancer_type)
    return(NULL)
  }
  
  # need cache for dealing with big data
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  cache_file <- file.path(cache_dir, paste0(cancer_code, "_processed.rds"))
  
  if (use_cache && file.exists(cache_file)) {
    msg <- paste0("Reading cached data for ", cancer_type, " from ", cache_file)
    message(msg)
    out <- tryCatch(
      readRDS(cache_file),
      error = function(e) {
        warning("Failed to read cache file for ", cancer_type, ": ", e$message)
        NULL
      }
    )
    if (!is.null(out)) return(out)
  }
  
  clinical_file <- file.path(data_dir, paste0("TCGA-", toupper(cancer_code), ".clinical.tsv.gz"))
  count_file <- file.path(data_dir, paste0("TCGA-", toupper(cancer_code), ".star_counts.tsv.gz"))
  survival_file <- file.path(data_dir, paste0("TCGA-", toupper(cancer_code), ".survival.tsv.gz"))
  
  if (!file.exists(clinical_file) ||
      !file.exists(count_file)    ||
      !file.exists(survival_file)) {
    warning("One or more TCGA files missing for ", cancer_type,
            ". Expected in directory: ", normalizePath(data_dir, winslash = "/"))
    return(NULL)
  }
  
  message("Building TCGA object for ", cancer_type, " which may take a while...")
  
  result <- tryCatch(
    build_tcga(
      clinical_file = clinical_file,
      count_file = count_file,
      survival_file = survival_file,
      top_n_genes = top_n_genes
    ),
    error = function(e) {
      warning("Failed to build TCGA object for ", cancer_type, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(result)) return(NULL)
  
  if (use_cache) {
    tryCatch(
      saveRDS(result, cache_file),
      error = function(e) {
        warning("Failed to save cache file for ", cancer_type, ": ", e$message)
      }
    )
  }
  
  result
}

get_cancer_object <- function(cancer_type) {
  if (is.null(cancer_type) || !nzchar(cancer_type)) return(NULL)

  if (exists(cancer_type, envir = .cancer_env, inherits = FALSE)) {
    return(get(cancer_type, envir = .cancer_env, inherits = FALSE))
  }
  
  obj <- load_cancer_data(cancer_type)
  if (is.null(obj)) {
    message("No data available for cancer type: ", cancer_type)
    return(NULL)
  }
  
  assign(cancer_type, obj, envir = .cancer_env)
  obj
}

get_cancer_data <- function(cancer_type) {
  obj <- get_cancer_object(cancer_type)
  if (is.null(obj)) return(NULL)
  obj$merged_data
}

get_cancer_expr <- function(cancer_type) {
  obj <- get_cancer_object(cancer_type)
  if (is.null(obj)) return(NULL)
  obj$expr_z
}

get_available_genes <- function(cancer_type) {
  expr_data <- get_cancer_expr(cancer_type)
  if (is.null(expr_data)) return(character(0))
  sort(rownames(expr_data))
}
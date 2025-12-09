OncoExploreR: An Interactive Shiny Application for Cancer Data
Exploration
================

**Authors:** Liuhan Ke · Ziyi Ou · Hailin Zhang

**Date:** 2025-12-06

<p align="center">

<img src="www/user_interface.png" alt="UI of OncoExploreR" width="700"><br>
<em>Project website:
<a href="https://lke-klh.shinyapps.io/OncoExploreR/">
https://lke-klh.shinyapps.io/OncoExploreR/</a></em>
</p>

# 1. Overview 🌟

OncoExploreR (Oncology Explorer) is an interactive Shiny application
designed to help users explore cancer genomics data from The Cancer
Genome Atlas (TCGA). The app integrates exploratory visualization,
differential expression gene detection benchmarking, and survival
modeling into a single, user-friendly interface aimed at learners and
beginners in cancer data analysis.

OncoExploreR allows users to select a cancer of interest out of 6
different types and provides:

- a body-map based entry point for visual and demographic exploration of
  each cancer,
- a simulation-based benchmark of four popular differential expression
  (DE) gene detection methods,
- a survival analysis module that illustrates how gene expression may
  affect predicted survival among certain demographic groups,
- direct access to cleaned TCGA-derived datasets for further offline
  analysis.

The underlying data include clinical variables (age, sex, tissue status,
race), survival outcomes (time-to-event and event indicators), and
preprocessed gene expression matrices (top variable genes based on
log2-transformed STAR counts).

The following table summarizes the cancer types and datasets included in
OncoExploreR:

| Organ / System    | Cancer Type                              |
|-------------------|------------------------------------------|
| Bronchus and Lung | Lung Adenocarcinoma (LUAD)               |
| Colon             | Colon Adenocarcinoma (COAD)              |
| Breast            | Breast Invasive Carcinoma (BRCA)         |
| Kidney            | Kidney Renal Clear Cell Carcinoma (KIRC) |
| Liver             | Liver Hepatocellular Carcinoma (LIHC)    |
| Thyroid           | Thyroid Carcinoma (THCA)                 |

# 2. Application Components 🧩

The application is structured into four distinct functional modules:
Body Map Explorer, DE Gene Detection Benchmark, Survival Analysis, and
Data Access.

## 2.1 Body Map Explorer 🩻

This module provides a visual and intuitive starting point for
exploratory data analysis. Users start with an interactive human body
map where they can click on specific organs to select the cancer type of
interest. Upon selection, the app displays:

- a brief textual overview of the disease,
- demographic summaries: distributions of age, sex, and race among
  patients
- top 20 variable genes heatmap.

This component is designed to help users situate each cancer in its
anatomical and population context before moving to more technical
analyses.

## 2.2 DE Gene Detection Benchmark 🧬

This module allows users to simulate differential expression scenarios
and benchmark the performance of four popular DE gene detection methods:
`DESeq2`, `edgeR`, `limma-trend`, and `limma-voom`. Users can customize
simulation parameters including: sample size per group, percentage of DE
genes, dispersion level and log$_2$ fold-change values.

Each simulated dataset is analyzed by all four methods, and the
application reports 3 plots to help users evaluate model performances:
ROC curves, Sensitivity–FDR trade-offs and Sensitivity vs logFC.

## 2.3 Survival Analysis 📉

This module uses TCGA primary tumor samples to study how gene expression
and clinical factors relate to predicted survival.

Users will input patient characteristics (age, sex) and a selected gene
with expression “treatment” value to compare against a baseline. The app
fits a Random Survival Forest (RSF) model to estimate survival curves
under both scenarios.

To incorporate biological and inter-individual variability, the module
simulates many individuals under each scenario and averages their
predicted survival curves.

This module outputs 2 plots: one is the average survival curves under
baseline vs treatment, with shaded ribbons reflecting confidence
intervals gained by simulation; the other is the distribution of
survival probability differences at 3, 5, and 10 years.

## 2.4 Data Access 🗃️

To facilitate reproducible research, OncoExploreR provides a data
download panel that exposes the processed datasets used in the
application. Users can download integrated tables containing clinical
metadata, survival data, and gene expression matrices for all six
supported cancer types.

# 3. Methodology 🔬

## 3.1 Differential Expression Modeling and Simulation 🔢

To simulate differential expression scenarios, we model raw counts using
a Negative Binomial (NB) distribution:
$$Y_{gi} \sim \mathrm{NB}(\mu_{gi}, \phi), \qquad
  \log(\mu_{gi}) = \beta_{0g} + \beta_{1g} X_i,$$ where:

- $Y_{gi}$ is the count of gene $g$ for sample $i$,
- $X_i$ is the group indicator (healthy vs. simulated patient),
- $\phi$ is the dispersion parameter,
- $\mu_{gi}$ is the mean expression level.

Baseline means $\mu_g$ are estimated from healthy samples. Differential
expression is introduced by multiplying the means with gene-specific
fold-changes:

- if $g$ is up-DE:  
  $\mu_{g,\text{patient}} = μ_g * 2^{\text{logFC}}$
- if $g$ is down-DE:  
  $\mu_{g,\text{patient}} = μ_g / 2^{\text{logFC}}$
- if $g$ is non-DE:  
  $\mu_{g,\text{patient}} = μ_g$

Users can select:

- sample size per group
- percentage of DE genes
- dispersion parameter

The benchmark evaluates four widely used DE tools:

- `DESeq2`: NB GLM with empirical shrinkage for dispersion and
  fold-change estimates.

- `edgeR`: NB model with empirical Bayes dispersion estimation and exact
  test.

- `limma-trend`: Linear modeling on log2(count+1) with mean–variance
  trend adjustment.

- `limma-voom`: Precision weights computed via voom transform, enabling
  linear modeling.

Each method returns statistics, p-values, and adjusted p-values for
every gene.

## 3.2 Survival Analysis with Random Survival Forests ⏳

To model patient-level outcomes, we fit a Random Survival Forest on TCGA
primary tumor samples using clinical variables and a selected gene:
$$ \mathrm{Surv}(T_i, \delta_i) \sim \mathrm{RSF}(\text{age}, \text{sex}, \text{gene expression}),$$
where $T_i$ is the survival time and $\delta_i$ is the event indicator.

The forest consists of $B = 200$ trees, each producing a survival
function $S_b(t \mid x)$ for covariate vector $x$. The ensemble
estimator is given by

$$ \hat{S}(t \mid x) = \frac{1}{B} \sum_{b=1}^B S_b(t \mid x).$$

Users can select:

- age
- sex
- gene name and its a gene expression value (in $\log_2$ scale) as
  treatment

We sample from the population to get the baseline value.

To reflect biological variability, we simulate 500 individuals for each
scenario: $$E_i \sim \mathcal{N}(E_{\text{mean}}, \sigma^2)$$ For each
individual, RSF produces a survival curve, and we obtain averaged curves
and distributions of survival probability differences at 3, 5, and 10
years.

This simulation framework illustrates the expected improvement in
survival probability given altered gene expression, while incorporating
cross-patient variability.

# 4. Installation 🚀

## 4.1 Use the Online Version 🌐

You can directly access the app via the deployed Shiny server:

    https://lke-klh.shinyapps.io/OncoExploreR/

## 4.2 Run Locally 🖥️

### Step 1 - Clone the repository

In your terminal:

``` bash
git clone https://github.com/Ziyi-Ou/OncoExploreR.git
cd OncoExploreR
```

### Step 2 — Restore the environment

Open R / RStudio inside the project directory:

``` r
install.packages("renv")

# Restore all package dependencies from renv.lock
renv::restore()
```

If you prefer not to use renv, you may manually install all required
packages:

``` r
install.packages(c(
  "shiny", "tidyverse", "survival", "randomForestSRC",
  "DESeq2", "edgeR", "limma", "pROC", "ggplot2",
  "plotly", "pheatmap", "tidyverse"
))
```

### Step 3 — Run the Shiny App

``` r
library(shiny)
shiny::runApp(".")
```

# 5. Acknowledgments 💡

We gratefully acknowledge the publicly available genomic data resources and open-source software that made this project possible.

### Data Resources
- **TCGA**: The Cancer Genome Atlas (https://portal.gdc.cancer.gov/)
- **UCSC Xena Browser**: (https://xenabrowser.net/datapages/)

### Key References

- Feng, H., Meng, G., Lin, T., Parikh, H., Pan, Y., Li, Z., ... & Li, Q. (2023). *Islet: individual-specific reference panel recovery improves cell-type-specific inference.* Genome Biology, 24(1), 174.

- Li, D., Zand, M. S., Dye, T. D., Goniewicz, M. L., Rahman, I., & Xie, Z. (2022). *An evaluation of RNA-seq differential analysis methods.* PLoS One, 17(9), e0264246.

- Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 15(12), 550.

- Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). *voom: Precision weights unlock linear model analysis tools for RNA-seq read counts.* Genome Biology, 15(2), R29.

- Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). *Random survival forests.*

- Ritchie, M. E., Phipson, B., Wu, D. I., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies.* Nucleic Acids Research, 43(7), e47.

- Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.* Bioinformatics, 26(1), 139-140.

### Software and R Packages

We used multiple open-source R packages to conduct statistical analysis, develop the user interface, and evaluate model performance, including:

- **DESeq2**: https://bioconductor.org/packages/DESeq2  
- **edgeR**: https://bioconductor.org/packages/edgeR  
- **limma**: https://bioconductor.org/packages/release/bioc/html/limma.html  
- **randomForestSRC**: https://cran.r-project.org/package=randomForestSRC  
- **shiny**: https://cran.r-project.org/web/packages/shiny/index.html  
- **pROC**: https://cran.r-project.org/package=pROC  
- **tidyverse**: https://www.tidyverse.org/packages/  

We sincerely thank the authors and developers of these tools for their contributions to open science and computational biology.


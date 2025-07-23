# RNA-Seq Analysis for Recruitment Task

## Project Description

This repository contains a complete analysis of RNA-Seq data, covering pre-processing, exploratory data analysis (EDA), differential gene expression (DE) analysis, and functional enrichment analysis. The project was completed as part of a recruitment task.

The data is derived from paired tumor and normal lung tissue samples from patients with different smoking statuses (`never`, `former`, `current`).

### Main Analysis Goals:

1.  **Identify differences in gene expression** between tumor and normal tissue.
2.  **Investigate differences in gene expression** between smokers, former smokers, and never-smokers.
3.  **Quantitatively assess the impact of smoking** on gene expression changes leading to carcinogenesis (interaction analysis).

### Methods Used:

* **Pre-processing:** Filtering of lowly expressed genes, TMM normalization.
* **Exploratory Data Analysis (EDA):** PCA, MDS, MA-plot, RLE plots.
* **Inference of Missing Metadata:** Sex identification based on the expression of Y-chromosome genes.
* **Differential Expression Analysis:** The `limma-voom` pipeline.
* **Enrichment Analysis:** Gene Ontology (GO) analysis using the `clusterProfiler` package.
* **Visualization:** Volcano plots, heatmaps, and UpSet plots.

## Repository Structure

* `Ardigen.Rmd`: The main R Markdown file containing all the code, visualizations, and commentary for the analysis.
* `Ardigen.html`: The rendered HTML report from the `Ardigen.Rmd` file, presenting the results in an accessible format.
* `Brudnopis.R`: An R script containing the analysis code (draft/development version).
* `renv.lock`: The `renv` lockfile ensuring full reproducibility of the analysis environment.

## How to Reproduce the Analysis

### 1. Download the Data

The input data files are required to run the analysis. Please download them and place them in the root directory of the project.

* **expression_data.tsv**: [Download from Google Drive](https://drive.google.com/file/d/1HzA2pgNVWmzYtqlWCL7PDx2UdwazicJs/view)
* **metadata.tsv**: [Download from Google Drive](https://drive.google.com/file/d/1T_1jMcS7C-3-aZ48JgdKoyRy_xp0iXAO/view)

### 2. Clone the Repository

```bash
git clone [https://github.com/michalkowalski94/ardigen-rnaseq.git](https://github.com/michalkowalski94/ardigen-rnaseq.git)
cd ardigen-rnaseq
```

### 3. Restore the `renv` Environment

This project uses `renv` for dependency management. To install all required packages in the exact versions used in the original analysis, open the project in R/RStudio and run the following command in the R console:

```r
renv::restore()
```

### 4. Generate the Report

To generate the HTML report, open the `Ardigen.Rmd` file in RStudio and click the "Knit" button, or run the following command in the R console:

```r
rmarkdown::render("Ardigen.Rmd")
```

Alternatively, you can run the plain R script:

```r
source("Brudnopis.R")

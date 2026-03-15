Two-Sample Mendelian Randomization: Inflammatory Proteins, Immune Cells, and Head and Neck Cancer


This repository contains the reproducible code for the two-sample Mendelian Randomization (MR) study: "Investigating the causal relationships between circulating inflammatory proteins, immune cell traits, and head and neck squamous cell carcinoma (HNSCC)".

Overview
This study systematically investigates causal relationships between 91 inflammatory proteins, 731 immune cell traits, and HNSCC risk using bidirectional two-sample MR. It employs a multi-threshold validation approach to identify robust, high-confidence associations.

Data Sources
All analyses use publicly available GWAS summary statistics:
HNSCC: IEU OpenGWAS project (ebi-a-GCST012235)
Circulating inflammatory proteins (91): GWAS Catalog (GCST90274758-GCST90274848)
Immune cell traits (731): GWAS Catalog (GCST90001908-GCST90002121)

Quick Start
To replicate the main analyses:
Clone the repository and install dependencies:
git clone https://github.com/KEERYREEK/Mendelian_HNC_IF_IC.git
cd Mendelian_HNC_IF_IC

# Install required R packages
Rscript -e "install.packages(c('TwoSampleMR', 'data.table', 'dplyr', 'ggplot2', 'forestploter', 'yaml'))"
Run the analysis pipeline:
The main workflow is executed by the master script:
Rscript scripts/run_all_analysis.R
This performs: data harmonization → IV selection → MR analysis → sensitivity tests → visualization.
Project Structure
  scripts/          # R analysis scripts
  data/             # Processed datasets (not included, see Data Availability)
  results/          # Output tables and figures (see our article =)
  config/           # Analysis parameters (p-value thresholds, LD clumping settings)
Key Outputs
The code generates:
MR Results: Causal estimates (OR, 95% CI, P-value) for all exposure-outcome pairs.
High-Confidence Associations: A filtered list meeting pre-specified criteria (significant across thresholds, passes sensitivity analyses, ≥3 IVs).
Complete Instrument Lists: Detailed SNP information for each analysis, including p-value thresholds, F-statistics, and LD clumping parameters.
Forest plots and scatter plots​ for visualization.
Data and Code Availability
In line with journal requirements for full reproducibility:
Harmonized Datasets & Full Instrument Lists: Uploading.
Analysis Code: This repository contains all R scripts. 
Please see the published article's Data Availability Statement for the final DOIs.

Citation
If you use this code or findings, please cite our publication:
not yet [Paper Title]. [Journal Name]. [Year]. DOI: [DOI]

Contact
For questions regarding the analysis, please open an issue on GitHub or contact the corresponding author.

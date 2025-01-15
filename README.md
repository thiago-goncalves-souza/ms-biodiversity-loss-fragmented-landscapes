
# Species turnover does not rescue biodiversity in fragmented landscapes
## Table of Contents

1. [Introduction](#introduction)
2. [Directory Structure](#directory-structure)
3. [Requirements](#requirements)
4. [Usage](#usage)
5. [Scripts Description](#scripts-description)
6. [Figures and Results](#figures-and-results)
7. [Acknowledgments](#acknowledgments)
8. [License](#license)

## Introduction

This is the repository with data and scripts to replicate all analyses run in the paper "Species turnover does not rescue biodiversity in fragmented landscapes"

## Directory Structure

The project is organized as follows:

```
- figures/                  # Contains the outputs of Figures 2 and 3
- processed_data/           # Folder for processed datasets after running analyses such as rarefaction
- renv/                     # Directory for R environment configurations
- results/                  # Summary and relevant results files
- 01_study_design.R         # Script for estimating alpha diversity based on study design
- 02_pairwise_div_all.R     # Script for estimating pairwise diversity (psd) using all plot pairs
- 03_pairwise_div_near.R    # Script for estimating pairwise diversity (psd) using the nearest plot pairs
- 04_pairwise_div_raref.R   # Script for estimating pairwise diversity (psd) using q = 0 and 2 (all pairs and nearest pairs)
- 05_pairwise_div_raref_samp_cov.R # Further rarefaction analysis to estimate pairwise diversity with a coverage based approach
- 06_figures.R              # Script for generating figures 2 and 3
- 07_glmm.R                 # Script for GLMMs to test the effects of landscape type on species diversity
- 08_meta_analysis.R        # Script for performing meta-analysis to test the effects of landscape type on species diversity
- Engel_novel_code.R        # Relevant functional to calculate beta diversity (coverage based) from Engel et al. 2021, Ecosphere, https://doi.org/10.1002/ecs2.3745
- biod_loss_frag.Rproj      # R project file 
- renv.lock                 # Lock file for R environment
- utility_functions.R       # Utility functions used across scripts (it was specified when needed across the scripts 01 to 08)
```

## Requirements

To run these scripts, you'll need:

- R and RStudio installed on your computer
- Packages defined in the `renv.lock` file
- Any additional software or dependencies required should be installed
- You must run the scripts from 01 to 08 to make sure all required files are created.

## Usage

You may want to download this repository (and associated R Project) to fully reproduce the scripts. All required packages and explanations for the usage can be run on each scripts, number from 01 (01_study_design.R) to 08 (08_meta_analysis.R). All required data is available in the folders data and processed_data. 

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/yourproject.git
   cd yourproject
   ```

2. Open `biod_loss_frag.Rproj` in RStudio to activate the R project environment.

3. Ensure all dependencies are installed, typically handled via the `renv` package:

   ```R
   renv::restore()
   ```

4. Run the analysis scripts in order as follows:

   - `source("01_study_design.R")`
   - `source("02_pairwise_div_all.R")`
   - (Repeat for other scripts...)

## Scripts Description

- **01_study_design.R**: Estimates alpha diversity based on study design
- **02_pairwise_div_all.R**: Estimates pairwise diversity (psd) using all plot pairs
- **03_pairwise_div_near.R**: Estimates pairwise diversity (psd) using the nearest plot pairs
- **04_pairwise_div_raref.R**: Estimates pairwise diversity (psd) using q = 0 and 2 (all pairs and nearest pairs)
- **05_pairwise_div_raref_samp_cov.R**: Further rarefaction analysis to estimate pairwise diversity with a coverage based approach
- **06_figures.R**: Generates figures 2 and 3
- **07_glmm.R**: Runs GLMMs to test the effects of landscape type on species diversity
- **08_meta_analysis.R**: Runs meta-analysis to test the effects of landscape type on species diversity

## Figures and Results

Images and results generated from the `06_figures.R` and other relevant scripts can be found in the `figures/` and `results/` directories.

## License

CC-BY-4.0 License
---

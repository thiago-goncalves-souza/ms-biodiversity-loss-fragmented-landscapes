
# Project Title

A brief description of what this project does and who it's for.

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

This is the script to replicate all analyses run in the paper "Increasing species turnover does not alleviate biodiversity loss in fragmented landscapes"

## Directory Structure

The project is organized as follows:

```
- figures/                  # Contains figure files for the project
- processed_data/           # Folder for processed datasets
- renv/                     # Directory for R environment configurations
- results/                  # Generated results files
- 01_study_design.R         # Script for setting up the study design
- 02_pairwise_div_all.R     # Script for pairwise diversity calculations (all)
- 03_pairwise_div_near.R    # Script for pairwise diversity calculations (near)
- 04_pairwise_div_raref.R   # Script for rarefaction diversity calculations
- 05_pairwise_div_raref_samp_cov.R # Further diversity analysis with sampling coverage
- 06_figures.R              # Script for generating figures
- 07_glmm.R                 # Script for Generalized Linear Mixed Models
- 08_meta_analysis.R        # Script for performing meta-analysis
- Engel_novel_code.R        # Novel code specific to Engel study
- biod_loss_frag.Rproj      # R project file for the study
- renv.lock                 # Lock file for R environment
- utility_functions.R       # Utility functions used across scripts
```

## Requirements

To run these scripts, you'll need:

- R and RStudio installed on your computer
- Packages defined in the `renv.lock` file
- Any additional software or dependencies required should be installed

## Usage

Instructions on how to replicate the results, including any installation commands or setup steps necessary for reproducing the analysis.

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

- **01_study_design.R**: Sets up the experimental design and initial parameters.
- **02_pairwise_div_all.R**: Computes pairwise diversity across the entire dataset.
- **03_pairwise_div_near.R**: Calculates pairwise diversity for nearby elements.
- **04_pairwise_div_raref.R**: Conducts rarefaction analyses.
- **05_pairwise_div_raref_samp_cov.R**: Incorporates sampling coverage in rarefaction.
- **06_figures.R**: Generates and saves plots and visualizations.
- **07_glmm.R**: Performs analyses using generalized linear mixed models.
- **08_meta_analysis.R**: Conducts meta-analysis using compiled data.
- **Engel_novel_code.R**: Contains novel methodologies developed during the Engel study.
- **utility_functions.R**: Supports all R scripts with reusable functions.

## Figures and Results

Images and results generated from the `06_figures.R` and other relevant scripts can be found in the `figures/` and `results/` directories.

## Acknowledgments

Acknowledge any collaborators, institutions, or repositories you leveraged during this project.

## License

Specify the license for the project, e.g., MIT License. Provide a link to the full license text if applicable.

---

This is a basic template for your `README.md`. Be sure to fill in specific details for each section to make the documentation useful to other users and collaborators.

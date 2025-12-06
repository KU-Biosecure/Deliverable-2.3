#  Biosecure Project – Hotspot Mapping Pipeline for Highly Pathogenic Avian Influenza (HPAI) and African Swine Fever (ASF)

This repository contains the full R-based pipeline for mapping and analyzing Avian Influenza (AIV) risk across NUTS2 regions in the EU. It integrates biosecurity metrics, environmental drivers, and disease case data to estimate regional risk scores and visualize hotspots. 

**Note the flow or process in the HPAI is the same ASF code**
---

## Project Structure

- **Data Sources**:
  - Biocheck data (2023 & 2024) for poultry types:
    - Free-range layers
    - Free-range broilers
    - Broilers
    - Laying hens
  - Bird density and land cover data
  - AIV poultry case counts
  - NUTS2 shapefiles and metadata

- **Main Outputs**:
  - Cleaned and merged datasets
  - Risk scores and disease introduction probabilities
  - Spatial visualizations of risk and scenario comparisons
  - Statistical modeling results (Negative Binomial regression)

---
## FAIR Data Principles
This repository follows the **FAIR data principles** to ensure that all resources are **Findable, Accessible, Interoperable, and Reusable**. The open-source R code is fully documented and version-controlled, and the harmonised datasets are structured to support reproducibility and integration with other tools and workflows. By aligning with FAIR standards, this project promotes transparency, encourages collaboration, and supports the broader research community in advancing spatial risk mapping for disease surveillance.

---
##  Pipeline Overview

### 1. Data Cleaning & Aggregation
- Reads and processes Biocheck Excel sheets
- Aggregates mean scores by region and question ID
- Merges 2023 and 2024 scores for each poultry type

### 2. Driver Integration
- Combines bird density and land cover with Biocheck scores
- Joins AIV case data to form the final modeling dataset

### 3. Modeling
- Fits a Negative Binomial model using `glmmTMB`
- Estimates dispersion parameter (theta) via profile likelihood
- Computes risk scores and disease introduction probabilities

### 4. Visualization
- Maps risk scores and disease introduction using `ggplot2` and `sf`
- Compares baseline vs. biosecurity improvement/decrease scenarios
- Highlights top 15 regions with disease introduction
---

## Dependencies

Make sure the following R packages are installed:


library(sf)
library(tidyverse)
library(glmmTMB)
library(readxl)
library(giscoR)
library(eurostat)
library(ggplot2)
library(gganimate)
library(gstat)
library(sp)
library(terra)
library(raster)
library(lattice)
library(httr)
library(utils)
library(cowplot)
library(stringr)
library(spdep)




 ## Author

This R code was developed and written by **Ofosuhene Okofrobour Apenteng**,for mapping biosecurity risks of African Swine Fever (ASF) and Highly Pathogenic Avian Influenza (HPAI) across EU regions. Developed for the Biosecure project (EU Horizon 101083923).



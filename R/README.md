# IncucyteFlow

<!-- badges: start -->
[![R-CMD-check](https://github.com/YOUR_USERNAME/IncucyteFlow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YOUR_USERNAME/IncucyteFlow/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Convert Incucyte microscopy object-level data to FCS format for analysis with flow cytometry tools.

## Overview

**IncucyteFlow** provides a pipeline to:

1. Parse Incucyte `.platemap` XML files to extract experimental metadata
2. Read and tidy Incucyte object-level CSV exports
3. Convert to FCS files compatible with FlowJo, flowCore, and other flow cytometry software
4. Load into GatingSet objects for downstream analysis with flowWorkspace and CytoExploreR

## Installation

### Dependencies

IncucyteFlow requires several Bioconductor packages. Install them first:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("flowCore", "flowWorkspace", "Biobase"))

# Optional but recommended:
BiocManager::install(c("CytoExploreR", "flowVS"))

# IncucyteFlow <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/bjreisman/IncucyteFlow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bjreisman/IncucyteFlow/actions/workflows/R-CMD-check.yaml)
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

# Recommended for downstream analysis:
BiocManager::install("CytoExploreR")
```

### Install IncucyteFlow

```r
# install.packages("devtools")
devtools::install_github("bjreisman/IncucyteFlow")
```

## Quick Start

```r
library(IncucyteFlow)
library(CytoExploreR)

# Step 1: Read platemap
platemap <- read_incucyte_platemap("experiment.PlateMap")

# Step 2: Read and tidy Incucyte data
mydata <- read_incucyte("Data")
mydata.tidy <- tidy_incucyte(mydata, platemap)

# Step 3: Create FCS files
fcs_files <- incucyte_to_fcs(
  data = mydata.tidy,
  output_dir = "incucyte_fcs",
  channel_desc = create_channel_desc(FLR1 = "GFP", FLR2 = "TMRM", FLR3 = "Sytox"),
  fluor_channels = c("FLR1MeanIntensity", "FLR2MeanIntensity", "FLR3MeanIntensity"),
  filename_prefix = "P493"
)

# Step 4: Load into GatingSet
gs <- load_incucyte_fcs("incucyte_fcs", platemap, prefix = "P493")

# Step 5: Compensate (if needed)
gs_comp <- cyto_compensate(gs, spillover = "Spillover-Matrix.csv")

# Step 6: Transform using logicle
channels <- c("FLR1MeanIntensity", "FLR2MeanIntensity", "FLR3MeanIntensity")
trans <- cyto_transformer_logicle(gs_comp, channels = channels, plot = FALSE)
gs_trans <- cyto_transform(gs_comp, trans = trans, plot = FALSE)

# Step 7: Gate and analyze with CytoExploreR
cyto_gate_draw(gs_trans, 
               parent = "root", 
               alias = "cells", 
               channels = c("Area", "Eccentricity"), 
               gatingTemplate = "gatingTemplate.csv")

cyto_gate_draw(gs_trans, 
               parent = "cells", 
               alias = "viable", 
               channels = "FLR3MeanIntensity", 
               gatingTemplate = "gatingTemplate.csv")
```

## Functions

| Function | Description |
|----------|-------------|
| `read_incucyte_platemap()` | Parse Incucyte .platemap XML files |
| `read_incucyte()` | Read Incucyte CSV data files |
| `tidy_incucyte()` | Tidy data and join platemap metadata |
| `create_channel_desc()` | Create channel description mapping |
| `incucyte_to_fcs()` | Convert to FCS files |
| `load_incucyte_fcs()` | Load FCS files into a GatingSet |
| `write_plater_csv()` | Export platemap to well-list CSV |
| `convert_to_plater_format()` | Convert to plater plate-shaped format |

## Expected Data Format

### Incucyte CSV Filenames

The `tidy_incucyte()` function expects CSV filenames in this format:

```
{dye}_{day}_{time}_{suffix}.csv
```

For example: `Phase_d0_0h00_Objects.csv`, `GFP_d1_24h00_Objects.csv`

### Platemap XML

Standard Incucyte `.platemap` files exported from the Incucyte software.

## Recommended Downstream Analysis

This package pairs well with [CytoExploreR](https://github.com/DillonHammill/CytoExploreR) for:

- Interactive gating (`cyto_gate_draw()`)
- Transformation (`cyto_transformer_logicle()`)
- Compensation (`cyto_compensate()`)
- Visualization and statistics

## License

MIT

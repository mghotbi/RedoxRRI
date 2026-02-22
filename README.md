# RedoxRRI  
**Holobiont Redox Resilience Index for Integrative Stress Biology**

---

## Overview

**RedoxRRI** is an R package for computing and visualizing a **Holobiont Redox Resilience Index (RRI)** by integrating **plant physiology**, **soil redox chemistry**, and **microbial resilience** into a unified, directionally identifiable framework.

The package is designed for applications in **redox ecology**, **plant–soil–microbe interactions**, and **holobiont resilience research**, with explicit support for spatial structure, temporal dynamics, and the statistical challenges typical of ecological data.

---

## Key Features

### Multidomain integration
- Plant physiological stress buffering (e.g. ROS-related traits)
- Soil redox chemistry and electron-acceptor stability proxies
- Microbial resilience derived from abundance data and/or ecological networks

### Directionally identifiable latent scores
- Dimension reduction via PCA, FA, NMF, WGCNA, and related approaches
- Explicit sign control with optional anchor variables to ensure reproducibility

### Robustness to ecological data realities
- Median imputation within spatial strata
- Tolerance to MNAR missingness, zero inflation, and collinearity

### Complementary outputs
- **Absolute domain scores**, scaled to \([0,1]\), suitable for statistical modeling
- **Compositional domain contributions**, summing to 1, designed for ternary visualization

---

## Conceptual Framework

Redox resilience at the holobiont scale emerges from **interacting biological subsystems**, rather than any single measurable trait.

**RedoxRRI** formalizes this concept by quantifying three domains:

- **Physio** — plant oxidative buffering and stress response  
- **Soil** — redox chemistry and electron-acceptor stability  
- **Micro** — microbial functional capacity and/or network resilience  

Each domain is summarized using a **latent variable** derived from multivariate indicators. Domain scores are explicitly oriented so that higher values consistently indicate greater resilience, scaled to \([0,1]\), and integrated into a composite **Redox Resilience Index (RRI)**.

---

## Design Principles

- **Modular** — domains can be added, removed, or reweighted  
- **Transparent** — all transformations are inspectable and reproducible  
- **Flexible** — supports multiple dimension-reduction strategies  
- **Mechanistic** — designed for inference, not black-box prediction  

---

## Methodological Components

### Dimension reduction
- Principal Component Analysis (PCA)
- Factor analysis
- Nonlinear embeddings (e.g. UMAP)
- Network-based summaries (e.g. WGCNA)

### Microbial resilience metrics
- Abundance- or function-based representations
- Optional network topology metrics derived from ecological graphs

### Aggregation
- User-defined weighting across biological domains
- Optional coupling terms to capture cross-domain coherence

---

## Typical Workflow

1. **Prepare domain-specific data**
   - plant physiological traits  
   - soil redox or chemical indicators  
   - microbial abundance, functional profiles, or networks  

2. **Derive latent scores** for each domain

3. **Integrate domains** into a composite Redox Resilience Index

4. **Evaluate and compare**
   - alternative weighting schemes  
   - stress or disturbance scenarios  
   - ecological or experimental contexts  

---

## Installation

```r
install.packages("remotes")
remotes::install_github("mghotbi/RedoxRRI")

# with vignettes
install.packages("remotes")
install.packages("BiocManager")

BiocManager::install(c(
  "BiocStyle",
  "rmarkdown",
  "knitr"))

remotes::install_github(
  "mghotbi/RedoxRRI",
  build_vignettes = TRUE,
  dependencies = TRUE)

```

## Example
```r
library(RedoxRRI)

# Simulate a holobiont redox system
sim <- simulate_redox_holobiont(seed = 1)
str(sim)

# Compute the Redox Resilience Index
res <- rri_pipeline_st(
  ROS_flux     = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data   = sim$micro_data,
  id           = sim$id,

  # Direction anchoring 
  direction_phys = "auto",
  direction_anchor_phys = "FvFm",
  direction_soil = "auto",
  direction_anchor_soil = "Eh"
)

# Absolute domain scores (for analysis)
head(res$row_scores)

# Compositional domain contributions (for ternary)
head(res$row_scores_comp)

# Ternary visualization
plot_RRI_ternary(res$row_scores_comp)


```

<img src="https://github.com/user-attachments/assets/abf56f94-582a-49d3-a28d-83087a01dd61"
     alt="Conceptual workflow of the Holobiont Redox Resilience Index (RRI)"
     width="900">

*Conceptual overview of the RedoxRRI framework. Multivariate indicators from plant physiology,
soil redox chemistry, and microbial systems are summarized into domain-level latent scores,
directionally aligned, scaled, and integrated into a unified holobiont-scale Redox
Resilience Index (RRI).*



## Output Structure

The primary output is an object of class RRI, containing:
**row_scores**
Absolute domain scores and per-sample RRI (all scaled to [0,1])
**row_scores_comp**
Compositional domain contributions (Physio + Soil + Micro = 1) and RRI
**meta**
Metadata including model settings and the system-level RRI index (rri_index)

The compositional scores are directly suitable for ternary visualization.

## Intended Use

RedoxRRI is designed for hypothesis-driven analysis of stress resilience in complex biological systems, with a focus on interpretability and mechanistic insight.
The package is particularly well suited for:
Ecological and environmental research
Quantifying resilience across abiotic stress gradients
Integrative stress biology
Soil–plant–microbiome interaction studies
Comparative and sensitivity analyses
Method development and exploratory modeling

## Design Philosophy

**RedoxRRI is not a black-box predictive model.**
Instead, it prioritizes:
explicit modeling choices
traceable transformations from raw data to index values
comparison of alternative hypotheses rather than automated optimization 
This makes RedoxRRI especially suitable for research contexts where understanding why resilience changes is as important as measuring how much.


#### Citation

If you use RedoxRRI in your research, please cite:

```r
citation("RedoxRRI")
```

Ghotbi, M. et al. RedoxRRI: A framework for holobiont-level redox resilience
(manuscript in preparation)


#### License

GPL-3 © 2025 Mitra Ghotbi

### Contact

**Maintainer:** Mitra Ghotbi  
📧 mitra.ghotbi@gmail.com  
🔗 ORCID: [0000-0001-9185-9993](https://orcid.org/0000-0001-9185-9993)




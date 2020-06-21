# SymSim
SymSim (Synthetic model of multiple variability factors for Simulation) is an R package for simulation of single cell RNA-Seq data. 

### Install from Github
This package can be installed with R package devtools. First, pull the package with git clone to your working directory. Make sure that you have installed the packages listed in the DESCRIPTION file.

The required Bioconductor packages can be installed as follows in R:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IRanges")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("S4Vectors")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```
To install SymSim, run:
```{r}
library("devtools")
devtools::install_github("YosefLab/SymSim")
```

### Generating datasets with SymSim

Refer to the [SymSim vignette](https://github.com/YosefLab/SymSim/blob/master/vignettes/SymSimTutorial.Rmd) for examples of using SymSim to simulate datasets.

### References

Xiuwei Zhang &ast;, Chenling Xu &ast;, Nir Yosef. **Simulating multiple faceted variability in Single Cell RNA sequencing**. _Nature Communications_, 10:2611, 2019. (https://www.nature.com/articles/s41467-019-10500-w).

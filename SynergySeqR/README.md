# SynergySeqR
An R package for SynergySeq functions

## Setup 
SynergySeqR is an R package that was developed and tested in [RStudio 1.1.442](https://www.rstudio.com/products/rstudio/) using [R 3.4.2](https://cran.rstudio.com/).

## Installing
To install the latest SynergySeqR package version from GitHub, use these commands:
```R
install.packages("devtools")
library(devtools)
devtools::install_github(“beronicao/SynergySeq", subdir = "SynergySeqR")
library(SynergySeqR)
```

## Uninstalling
```R
remove.packages("SynergySeqR")
```

## Updating
```R
library(devtools)
update_packages("SynergySeqR")
```

## Demo
```R
res <- Drugs_SigsR()

res[[2]]$preRenderHook(p=res[[2]])  # view plot

```

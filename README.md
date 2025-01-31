<p style="text-align: center;">
  <img src="man/figures/logo.png" alt="Example Image" width="70">
</p>

# vmrseq: Detecting variably methylated regions (VMRs) from single-cell bisulfite sequencing

The R package `vmrseq` is a novel computational tool developed for pinpointing variably methylated regions (VMRs) in scBS-seq data without prior knowledge on size or location. High-throughput single-cell measurements of DNA methylation allows studying inter-cellular epigenetic heterogeneity, but this task faces the challenges of sparsity and noise. vmrseq overcomes these challenges and identifies variably methylated regions accurately and robustly. 

![](man/figures/method.png)


## Installation

You can install the development version of `vmrseq` in R from Bioconductor (recommended) or GitHub with:

``` r
### Install stable version from Bioconductor
 if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("vmrseq")

## Or development version from Github
# install.packages("remotes")
remotes::install_github("nshen7/vmrseq")
```

## Online Vignette

- An online vignette of 'Get Started' on the `vmrseq` package can be found at 
[https://rpubs.com/nshen7/vmrseq-vignette](https://rpubs.com/nshen7/vmrseq-vignette).
- An example workflow on scBS-seq data analysis using `vmrseq` can be found at
[https://github.com/nshen7/vmrseq-workflow-vignette](https://github.com/nshen7/vmrseq-workflow-vignette).

## Docker Image

We provide a Docker image for robust setup and use of this package. The Docker image includes all necessary dependencies for the package and vignettes. To pull the Docker image from Docker Hub, use the following command in bash:

``` bash
docker pull nshen7/vmrseq-bioc-3.19:latest
```


## Citation

Shen, N., Korthauer, K. vmrseq: probabilistic modeling of single-cell methylation heterogeneity. Genome Biol 25, 321 (2024). https://doi.org/10.1186/s13059-024-03457-7

## License/Copyright

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 
This package is made available under an MIT license.  

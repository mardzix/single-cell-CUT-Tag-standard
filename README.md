# Project Title

This is a basic single-cell CUT&Tag preprocessing pipeline

## Description

The pipeline is based on snakemake workflow manager.
The pipeline takes output from cellranger-atac as input and performs following:
- Peak calling using MACS2 instead of standard cellranger algorithm
- Generate a bulk bigwig track (all cells) for QC purposes
- Performs cell selection more suited for scCUT&Tag data
- Creates seurat object with a number of various feature matrices:
  - bins of various sizes (5kb-250kb)
  - gene body and promoter matrix
  - promoter matrix
  - peak matrix
- Attempts dimensionality reduction (LSI/UMAP) using Seurat/Signac
## Getting Started

Prepare conda environment. 
```
conda env create -f envs/environment.yaml
```

### Dependencies

TODO - add from environment file

TODO - add R packages 


### Installing

Clone this repo
```
git clone https://github.com/mardzix/single-cell-CUT-Tag-standard.git
```


### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help


## Authors

Marek Bartosovic
[@marekbartosovic](https://twitter.com/marekbartosovic)

## Version History


## License
<!--- 
(This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details)
--->

## Acknowledgments

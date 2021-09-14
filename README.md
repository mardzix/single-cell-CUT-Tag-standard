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

### Dependencies

- python3
- 


### Installing

* How/where to download your program
* Any modifications needed to be made to files/folders

### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
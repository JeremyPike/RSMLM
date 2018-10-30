# RSMLM

A R package for pointillist analysis of single molecule localization microscopy (SMLM) data. This package includes methods for persistence based clustering (ToMATo), alongside DBSCAN, Voronoi tessellation and Ripley's K based clustering. There is also the capacity to simulate dSTORM data.

## Tutorials

A set of Binder ready tutorials and demos showcasing this package can be found within the [SMLM-tutorials](https://github.com/JeremyPike/RSMLM-tutorials) repository. We recommend looking at these before installing the package.  

## Installation
This package is under active devlopment and we havn't released a version on CRAN yet.

The **devtools** package can be used to install the development version directly from github:

```r
install_github("jeremypike/RSMLM")
```

To run this command you will need [devtools](https://cran.r-project.org/package=devtools). Windows users will also need [Rtools](http://www.murdoch-sutherland.com/Rtools/).

## Citation

If RSMLM is useful please cite our pre-print:

Topological data analysis quantifies biological nano-structure from single molecule localization microscopy\
Jeremy A Pike, Abdullah O Khan, Chiara Pallini, Steven G Thomas, Markus Mund, Jonas Ries, Natalie S Poulter, Iain B Styles\
bioRxiv 400275; doi: [https://doi.org/10.1101/400275](https://doi.org/10.1101/400275)

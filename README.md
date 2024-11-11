# RClass2024

<!--   주석처리  toc:table of content    
author: "Sung Young Kim, MD, PhD^[mailto:palelamp@gmail.com]</center>"
date: "03 June 2020"
output: 
    html_document:
        toc: false  
        toc_float: 
           collapsed: false
        toc_depth: 3
        
output: 
    md_document:
        variant: markdown_github
-->
<!--
This page provides the supplementary R code and data to reproduce the experiments in the following paper : **Accurate prediction of acquired EGFR TKIs resistance using a pathway-based individualized machine learning approach**  
-->

------------------------------------------------------------------------

**Main Dataset**
----------------

-   The main method function R file can be downloaded from [here](http://centromics.org/info/142sup/mainFunctions.R)
-   Preprocessed dataset can be downloaded from [here](http://centromics.org/info/142sup/EGFRTKIs_8set.RData)
-   The pathways used for model construction can be downloaded from [here](http://centromics.org/info/142sup/p.KEGG.PID.BioCarta.RData)

------------------------------------------------------------------------

**Package Download**
--------------------

-   Package source file can be downloaded from [here](http://centromics.org/info/142sup/mainFunctions.R)

------------------------------------------------------------------------

**Install Dependencies**
------------------------

-   If don't have the R software installed in our computer, download and install it (check out the [R home page](http://www.r-project.org/))
-   Open the R command line interface, and install package dependencies (if they have not been installed yet):
-   Dependencies: R (&gt;= 3.0.0), shiny (&gt;= 0.8.0), WGCNA, igraph, shinyBS, RColorBrewer, Hmisc, psych, RJSONIO, whisker, yaml, pheatmap, preprocessCore, GO.db, AnnotationDbi, impute, and ggplot2

``` r
if (!requireNamespace("BiocManager", quietly = T)) install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")
install.packages(c("ggplot2", "ggrepel", "WGCNA", "igraph", "Hmisc"))
```

-   Please, install versions 0.8.0 for shiny. <!--We are working to make the package compatible with the new versions of the packages as soon as possible.--> To install the recommended versions for shiny, just type the following commands on the R command-line:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_version("shiny", "0.8.0")
```

------------------------------------------------------------------------

**Install package via Github (Recommended)**
--------------------------------------------

To install the latest version of package via Github, run following commands in R:

``` r
devtools::install_github("malcogene/RClass")
```

------------------------------------------------------------------------

**Install package from the source**
-----------------------------------

-   **Linux/Mac OS**
    -   Download the package \*.tar.gz.
    -   Open a command prompt. then:

``` r
install.packages(path_of_the_downloaded_file, repos = NULL, type="source")
```

-   **Windows**
    -   Start by reviewing the section on Windows packages in [the R Installation and Administration manual](https://cran.r-project.org/doc/manuals/R-admin.html), then carefully follow the instructions from [The Windows toolset appendix](https://cran.r-project.org/doc/manuals/R-admin.html#The-Windows-toolset).
    -   Download the package \*.tar.gz.
    -   Make sure you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed.
    -   Open a command prompt. then:

``` r
install.packages(path_of_the_downloaded_file, repos = NULL, type="source")
```

------------------------------------------------------------------------

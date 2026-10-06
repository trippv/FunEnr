# FunEnr

<!-- badges: start -->
<!-- badges: end -->

## Overview

**FunEnr** is an R package for performing Gene Ontology (GO) functional enrichment analysis using the [`topGO`](https://bioconductor.org/packages/topGO/) framework.

The package provides a simple workflow for:

- performing GO enrichment analysis for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC);
- preparing gene-to-GO annotation tables;
- automatically detecting comma- or semicolon-separated GO identifiers;
- identifying candidate genes associated with enriched GO terms; and
- reducing redundant GO terms based on semantic similarity.

The main function, `FunEnr_Topgo()`, is designed to simplify the use of `topGO` when working with gene lists and custom gene-to-GO annotations.

## Installation

You can install the development version of **FunEnr** from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("trippv/FunEnr")
```

The package depends on several Bioconductor packages. If they are not already installed, they can be installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "topGO",
  "GOSemSim",
  "org.Hs.eg.db"
))
```

## Basic workflow

Load the package:

```r
library(FunEnr)
```

Prepare a gene-to-GO annotation table. The table must contain two columns: one with gene identifiers and one with GO identifiers.

For example:

```r
background <- data.frame(
  genes = c("gene1", "gene2", "gene3", "gene4"),
  GO = c(
    "GO:0008150,GO:0009987",
    "GO:0008150",
    "GO:0003674,GO:0008150",
    "GO:0009987"
  )
)
```

Define the genes of interest:

```r
genelist <- c("gene1", "gene2")
```

Run the enrichment analysis:

```r
results <- FunEnr_Topgo(
  genelist = genelist,
  background = background,
  ontology = "BP"
)
```

The resulting data frame contains the enriched GO terms together with their enrichment statistics, adjusted p-values, and the candidate genes associated with each term.

## GO ontologies

The `ontology` argument allows enrichment analysis for the three main Gene Ontology domains:

```r
results_BP <- FunEnr_Topgo(
  genelist,
  background,
  ontology = "BP"
)

results_MF <- FunEnr_Topgo(
  genelist,
  background,
  ontology = "MF"
)

results_CC <- FunEnr_Topgo(
  genelist,
  background,
  ontology = "CC"
)
```

## Reducing redundant GO terms

GO enrichment analyses can return multiple related or redundant terms. FunEnr can optionally reduce these terms using semantic similarity:

```r
results_reduced <- FunEnr_Topgo(
  genelist = genelist,
  background = background,
  ontology = "BP",
  reduce_terms = TRUE
)
```

When `reduce_terms = TRUE`, semantic similarity is calculated using `GOSemSim` and redundant terms are reduced using `rrvgo`.

The current implementation uses `org.Hs.eg.db` for semantic similarity calculations and therefore uses Homo sapiens annotation data for this step.

## Input format

The background annotation must contain exactly two columns:

```text
gene        GO
gene1       GO:0008150,GO:0009987
gene2       GO:0008150
gene3       GO:0003674
```

GO identifiers may be separated by either commas or semicolons:

```text
GO:0008150,GO:0009987
```

or

```text
GO:0008150;GO:0009987
```

The order of the columns does not matter. FunEnr automatically identifies the column containing GO identifiers.

## Documentation

A detailed introduction to the package and examples of its workflow are available in the package vignette:

```r
vignette("FunEnr")
```

## Development

FunEnr is under active development. Contributions, bug reports, and suggestions are welcome through the GitHub repository.

## License

FunEnr is released under the MIT License.

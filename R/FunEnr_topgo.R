#' Perform GO Enrichment Analysis with topGO
#'
#' This function performs Gene Ontology (GO) enrichment analysis using the `topGO` package for all GO ontologies (BP, MF, CC). It can optionally reduce redundant GO terms.
#'
#' @param genelist A character vector of differentially expressed genes to be tested for enrichment.
#' @param background A data frame with two columns: the first column contains gene names, and the second column contains corresponding GO IDs.
#' @param ontology A character string specifying the ontology to use for analysis. Valid options are "BP" (Biological Process), "MF" (Molecular Function), and "CC" (Cellular Component). Defaults to "BP".
#' @param reduce_terms Logical, whether to reduce redundant GO terms in the results. Defaults to `FALSE`.
#' @param topnode An integer specifying the number of top GO terms (nodes) to report. Defaults to 30.
#' @return A data frame containing the enrichment analysis results, including GO terms, p-values, and scores.
#' @details This function uses the `topGO` package to perform the analysis. It takes a gene list, compares it against the background GO terms, and identifies enriched GO categories based on the given ontology.
#' @export

FunEnr_Topgo <- function(genelist,
                         background,
                         ontology = "BP",
                         topnode = 30,
                         reduce_terms = FALSE) {

  # Load required libraries
  require(topGO)
  require(rrvgo)
  require(dplyr)
  require(tidyr)
  require(GOSemSim)
  require(org.Hs.eg.db)

  # Validate inputs
  if (!is.character(genelist) || length(genelist) == 0) {
    stop("genelist must be a non-empty character vector.")
  }

  # Check if the background data has exactly two columns
  if (ncol(background) != 2) {
    stop("Background data is not in the correct format. It must have exactly two columns.")

  } else {
    # Check for "GO:" in the first or second column
    if (any(str_detect(background[, 1], "GO:"))) {
      print("GO terms detected in the first column")

      # Change column names: "GO" should be in the second column
      colnames(background)[1] <- "GO"
      colnames(background)[2] <- "genes"

    } else if (any(str_detect(background[, 2], "GO:"))) {
      print("GO terms detected in the second column")

      # Change column names: "GO" should be in the second column
      colnames(background)[2] <- "GO"
      colnames(background)[1] <- "genes"

    } else {
      stop("No GO terms detected in either column.")
    }
  }

  # Validate ontology
  if (!ontology %in% c("BP", "MF", "CC")) {
    stop("Invalid ontology. Choose from 'BP', 'MF', or 'CC'.")
  }

  # Detect separator for GO terms
  first_go_string <- background$GO[1]
  separator <- detect_separator(first_go_string)
  print(paste0("Separator detected: '", separator, "'"))

  # Set ontology
  message(paste0("Ontology set to: ", ontology))


  # Prepare gene-to-GO mapping
  gene_2_go <- background |>
    separate_rows(GO, sep = separator)

  # unstack GO terms for each gene
  gene_2_go <- unstack(gene_2_go[, c("GO", "genes")])

  # Filter candidate genes to keep only those present in the background
  candidate_list <- genelist[genelist %in% background$genes]

  # Create factor for the gene list
  bg_genes <- as.character(background$genes)
  geneList <- factor(as.integer(bg_genes %in% candidate_list))
  names(geneList) <- bg_genes

  # Create TopGO data object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = gene_2_go)

  # Run enrichment analysis
  weight_fisher_result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  # Extract GO results
  allGO <- usedGO(GOdata)
  all_res <- GenTable(GOdata, weightFisher = weight_fisher_result, topNodes = topnode)


  # Memory management: Clear unnecessary objects
  rm(weight_fisher_result)
  rm(gene_2_go)
  gc()  # Garbage collection

  # Extract significant genes for each GO term
  GOnames <- as.vector(all_res$GO.ID)
  allGenes <- genesInTerm(GOdata, GOnames)
  significantGenes <- list()

  for(x in 1:nrow(all_res)){

    significantGenes[[x]] <- allGenes[[x]][allGenes[[x]] %in% as.vector(candidate_list)]
  }
  names(significantGenes) <- all_res$Term


  # Add significant genes to the results table
  all_res$Genes <- sapply(significantGenes, function(g) paste(g, collapse = ","))

  # Adjust p-values for multiple testing
  all_res <- all_res %>%
    mutate(p.adj = p.adjust(weightFisher, method = "BH")) %>%
    arrange(p.adj)

  # Reduce GO terms using rrvgo, if required
  if (reduce_terms) {
    message("Reducing GO terms using rrvgo.")

    semdata <- GOSemSim::godata("org.Hs.eg.db", ont = ontology)
    simMatrix <- calculateSimMatrix(all_res$GO.ID,
                                    orgdb = "org.Hs.eg.db",
                                    ont = ontology,
                                    method = "Rel",
                                    semdata = semdata)

    scores <- setNames(-log10(all_res$p.adj), all_res$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.6, orgdb = "org.Hs.eg.db")

    # Filter only parent terms
    all_res <- all_res[all_res$GO.ID %in% reducedTerms$parent, ]

  }

  return(all_res)
}

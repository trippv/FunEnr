#' Perform Gene Ontology enrichment analysis using topGO
#'
#' Performs Gene Ontology (GO) enrichment analysis using the
#' \pkg{topGO} package. The function tests for enrichment of GO terms
#' among a set of candidate genes relative to a user-provided
#' gene-to-GO background annotation.
#'
#' The analysis can be performed for Biological Process (BP),
#' Molecular Function (MF), or Cellular Component (CC) ontologies.
#' When requested, redundant GO terms can be reduced using semantic
#' similarity with \pkg{rrvgo} and \pkg{GOSemSim}.
#'
#' @param genelist A character vector containing the genes of interest,
#'   typically a set of differentially expressed genes.
#' @param background A data frame containing two columns representing
#'   a gene-to-GO annotation. One column must contain gene identifiers
#'   and the other GO identifiers. GO identifiers can be provided as
#'   comma- or semicolon-separated strings.
#' @param ontology A character string specifying the GO ontology to
#'   analyze. One of `"BP"` (Biological Process), `"MF"` (Molecular
#'   Function), or `"CC"` (Cellular Component). Defaults to `"BP"`.
#' @param topnode Integer specifying the maximum number of GO terms
#'   returned by \code{\link[topGO]{GenTable}}. Defaults to 30.
#' @param reduce_terms Logical indicating whether redundant GO terms
#'   should be reduced using semantic similarity. Defaults to `FALSE`.
#'
#' @return A data frame containing the GO enrichment results. The
#'   returned table includes GO identifiers, GO term descriptions,
#'   enrichment statistics, adjusted p-values, and the genes from the
#'   candidate list associated with each GO term.
#'
#' @details
#' The function first identifies the gene and GO columns in the
#' supplied background annotation and automatically detects whether
#' multiple GO identifiers are separated by commas or semicolons.
#'
#' Candidate genes are restricted to genes present in the supplied
#' background annotation. A \pkg{topGO} data object is then constructed
#' using the gene-to-GO mapping and enrichment is tested using the
#' `weight01` algorithm with Fisher's exact test.
#'
#' P-values are adjusted for multiple testing using the Benjamini-Hochberg
#' method.
#'
#' If `reduce_terms = TRUE`, semantically similar GO terms are evaluated
#' using \pkg{GOSemSim} and redundant terms are reduced with
#' \pkg{rrvgo}. Semantic similarity is calculated using the
#' \pkg{org.Hs.eg.db} annotation database and is therefore based on
#' Homo sapiens GO annotations.
#'
#' @section Input format:
#' The background annotation must contain exactly two columns. For
#' example:
#'
#' \preformatted{
#' gene        GO
#' gene_001    GO:0008150,GO:0009987
#' gene_002    GO:0003674
#' gene_003    GO:0005575;GO:0008150
#' }
#'
#' The order of the two columns is not important as long as one column
#' contains GO identifiers beginning with `"GO:"`.
#'
#' @examples
#' # Example gene-to-GO annotation
#' background <- data.frame(
#'   genes = c("gene1", "gene2", "gene3", "gene4", "gene5"),
#'   GO = c(
#'     "GO:0008150,GO:0009987",
#'     "GO:0008150",
#'     "GO:0003674,GO:0008150",
#'     "GO:0009987",
#'     "GO:0003674"
#'   )
#' )
#'
#' # Candidate genes
#' genes <- c("gene1", "gene2")
#'
#' # Run enrichment analysis
#' \dontrun{
#' results <- FunEnr_Topgo(
#'   genelist = genes,
#'   background = background,
#'   ontology = "BP"
#' )
#' }
#'
#' @seealso
#' \code{\link[topGO]{topGOdata}},
#' \code{\link[topGO]{runTest}},
#' \code{\link[topGO]{GenTable}},
#' \code{\link{detect_separator}}
#'
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

#' Detect the Separator in a GO Term String
#'
#' This function detects the separator used in a string of Gene Ontology (GO) terms. It checks whether the terms are separated by commas (`,`), semicolons (`;`), or if no valid separator is found.
#'
#' @param go_string A character string containing GO terms.
#' @return A character value indicating the detected separator (either `","` or `";"`).
#' @examples
#' # Example: GO terms separated by commas
#' detect_separator("GO:0008150,GO:0003674,GO:0005575")
#'
#' # Example: GO terms separated by semicolons
#' detect_separator("GO:0008150;GO:0003674;GO:0005575")
#'
#' # This will throw an error:
#' \dontrun{
#' detect_separator("GO:0008150 GO:0003674 GO:0005575")
#' }
#' @export
detect_separator <- function(go_string) {
  if (stringr::str_detect(go_string, ",")) {
    return(",")
  } else if (stringr::str_detect(go_string, ";")) {
    return(";")
  } else {
    stop("No valid separator detected in GO terms")
  }
}

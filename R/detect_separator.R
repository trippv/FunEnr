#' Detect the separator used between GO identifiers
#'
#' Detects whether multiple Gene Ontology (GO) identifiers in a
#' character string are separated by commas or semicolons.
#'
#' @param go_string A character string containing one or more GO
#'   identifiers.
#'
#' @return A character string containing the detected separator:
#'   `","` for comma-separated GO identifiers or `";"` for
#'   semicolon-separated GO identifiers.
#'
#' @details
#' The function searches for a comma first and then for a semicolon.
#' If neither separator is present, the function returns an error.
#'
#' @examples
#' detect_separator("GO:0008150,GO:0003674,GO:0005575")
#'
#' detect_separator("GO:0008150;GO:0003674;GO:0005575")
#'
#' \dontrun{
#' detect_separator("GO:0008150 GO:0003674 GO:0005575")
#' }
#'
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

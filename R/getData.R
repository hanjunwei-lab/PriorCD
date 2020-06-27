##' This function is used to get example data.
##'
##' @title getData
##' @param exampleData String. These example data are included: mRNA_path, microRNA_path, drug.ic50, drug.r, drug.fdr, drug.info, drug.edgelist, breast_cancer and brc_candidates.
##' @export
getData <- function(exampleData) {
  if (!exists("envData")) {
    envData <- initializePriorCD()
  }

  if (exampleData == "drug.ic50") {
    dataset <- get("drug.ic50", envir = envData)
    return(dataset)
  }

  if (exampleData == "mRNA_path") {
    dataset <- get("mRNA_path", envir = envData)
    return(dataset)
  }

  if (exampleData == "microRNA_path") {
    dataset <- get("microRNA_path", envir = envData)
    return(dataset)
  }

  if (exampleData == "drug.edgelist") {
    dataset <- get("drug.edgelist", envir = envData)
    return(dataset)
  }

  if (exampleData == "drug.info") {
    dataset <- get("drug.info", envir = envData)
    return(dataset)
  }

  if (exampleData == "breast_cancer") {
    dataset <- get("breast_cancer", envir = envData)
    return(dataset)
  }

  if (exampleData == "drug.r") {
    dataset <- get("drug.r", envir = envData)
    return(dataset)
  }

  if (exampleData == "drug.fdr") {
    dataset <- get("drug.fdr", envir = envData)
    return(dataset)
  }

  if (exampleData == "brc_candidates") {
    dataset <- get("brc_candidates", envir = envData)
    return(dataset)
  }

  if (exampleData == "priorlist") {
    dataset <- get("priorlist", envir = envData)
    return(dataset)
  }
}

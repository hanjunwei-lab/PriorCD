initializePriorCD <- function() {
  utils::data("envData", package = "PriorCD")
}
getenvir <- function(envData) {
  if (!exists("envData"))
    initializePriorCD()
  return(get(envData, envir = envData))
}

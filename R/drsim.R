##' This function is used to construct a binary adjacency matrix of drug similarity where 1 means strong similarity and 0 means weak similarity.
##'
##' @title drsim
##' @param r.mat The input matrix of drug correlations.
##' @param p.mat The input matrix of probability values(p-value) of drug correlations.
##' @param top A value to measure drug similarity. It's a threshold of correlation, top=0.005(default) means that top 0.005 of drugs for each row are considered as strong similarity.
##' @param r.thres A value to measure drug similarity. It's a threshold of correlation, r.thres=0.7(default) means that the similarity between drugs are strong when r greater than 0.7.
##' @param p.thres A value to measure the significance level of drug similarity. It's a threshold of probability values, p.thres=0.01(default) means that the similarity between drugs are significant when p less than 0.01.
##' @importFrom grDevices dev.off
##' @importFrom grDevices pdf
##' @importFrom graphics abline
##' @importFrom graphics legend
##' @return A binary adjacency matrix of drug similarity.
##' @examples
##' r <- getData("drug.r")
##' fdr <- getData("drug.fdr")
##' m <- drsim(r, fdr, top = 0.5)
##' @export
drsim <- function(r.mat, p.mat, top = 0.005, r.thres = 0.7, p.thres = 0.01) {
  drug.adjm <- matrix(data = 0, nrow = nrow(r.mat), ncol = ncol(r.mat))
  colnames(drug.adjm) <- colnames(r.mat)
  rownames(drug.adjm) <- rownames(r.mat)
  rank <- ceiling(ncol(r.mat) * (top))

  pb <- txtProgressBar(min = 0, max = nrow(r.mat), style = 3)

  for (i in 1:nrow(r.mat)) {
    value <- sort(r.mat[i, ], decreasing = T)[rank]
    drug.adjm[i, which(r.mat[i, ] >= as.numeric(value))] <- 1

    setTxtProgressBar(pb, i)
  }
  close(pb)

  drug.adjm[r.mat < r.thres] <- 0
  drug.adjm[p.mat > p.thres] <- 0
  t.drug.adjm <- t(drug.adjm)
  drug.adjm <- t.drug.adjm + drug.adjm
  drug.adjm[drug.adjm == 2] <- 1
  diag(drug.adjm) <- 0
  return(drug.adjm)
}

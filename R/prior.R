##' This function is used to generate drug prioritizing result.
##'
##' @title prior
##' @param drug.el A edge list of drugs, which is a two-column matrix, each row defines one edge. Numbers in the edge list represent NSC-ID of drugs.
##' @param p0 A vector of approved drugs' NSC-ID of interested cancer.
##' @param gamma gamma = 0.7(default). A probability of losing when doing Random Walk. On the contray, there is a probability of 1-gamma left to itself. The range of this value is (0, 1).
##' @param times times = 100(default). Loop times when getting p-values.
##' @return Detailed information about drug prioritizing, which contain NSC-id, name, prioritizing score, p-value, FDR, status and MOA(mechanism of action) of drugs.
##' @importFrom stats p.adjust
##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @importFrom igraph graph_from_edgelist
##' @importFrom igraph as_adj
##' @importFrom dplyr %>%
##' @importFrom dplyr arrange
##' @importFrom dplyr desc
##' @examples
##' e <- getData("drug.edgelist")
##' brc <- getData("breast_cancer")
##' \donttest{result <- prior(e, brc,time=20)}
##' @export

prior <- function(drug.el, p0, gamma = 0.7, times = 100) {
  drug.el <- apply(drug.el, 2, as.character)
  g <- graph_from_edgelist(drug.el, directed = F)
  iNet <- as_adj(g, type = "both", names = T, sparse = F)

  new.p0 <- matrix(data = 0, ncol = 1, nrow = nrow(iNet))
  rownames(new.p0) <- row.names(iNet)
  new.p0[which(rownames(new.p0) %in% p0), 1] <- 1

  real.pt <- rw(as.matrix(iNet), new.p0[, 1], gamma)
  real.pt[,1] <- (as.data.frame(real.pt) %>% arrange(as.numeric(rownames(real.pt))))[,1]
  rownames(real.pt) <- sort(as.numeric(rownames(real.pt)))

  print("Generating table...")
  rp <- function(times) {
    rp.pt <- matrix(data = 0, nrow = nrow(iNet), ncol = times)
    rownames(rp.pt) <- sort(as.numeric(rownames(iNet)))

    pb <- txtProgressBar(min = 0, max = times, style = 3)

    for (i in 1:times) {
      names <- sample(rownames(iNet), nrow(iNet), replace = F)
      rownames(iNet) <- names
      colnames(iNet) <- names
      rp.p0 <- matrix(data = 0, ncol = 1, nrow = nrow(iNet))
      rownames(rp.p0) <- names
      rp.p0[which(rownames(rp.p0) %in% p0), ] <- 1
      index <- rw(as.matrix(iNet), rp.p0[, 1], gamma)
      rp.pt[, i] <- (as.data.frame(index) %>% arrange(as.numeric(rownames(index))))[, 1]

      setTxtProgressBar(pb, i)
    }
    close(pb)

    rp.p <- vector(length = nrow(iNet), mode = "numeric")
    for (i in 1:length(rp.p)) {
      rp.p[i] <- length(which(rp.pt[i, ] >= real.pt[i, ]))/as.numeric(times)
    }

    return(rp.p)
  }
  p <- rp(times)
  fdr <- p.adjust(p, method = "BH", n = length(p))

  m <- data.frame(NSC.ID = 1:nrow(iNet))
  drug.info <- getData("drug.info")
  m$NSC.ID <- sort(as.numeric(rownames(iNet)))
  m$DRUG.NAME <- drug.info[match(m$NSC.ID, drug.info[,1]), 2]
  m$PRI.SCORE <- real.pt[,1]
  m$P.VALUE <- p
  m$FDR <- fdr
  m$STATUS <- drug.info[match(m$NSC.ID, drug.info[,1]), 3]
  m$MOA <- drug.info[match(m$NSC.ID, drug.info[,1]), 4]
  m <- m[-which(m$NSC.ID %in% p0), ]
  m <- m %>% arrange(m$FDR, desc(m$PRI.SCORE))

  return(m)
}

rw <- function(drug.adjm, p0, gamma) {
  p0 <- t(p0)
  p0 <- p0/sum(p0)
  PT <- p0
  k <- 0
  delta <- 1
  Ng <- dim(drug.adjm)[2]
  for (i in 1:Ng) {
    sumr <- sum(drug.adjm[i, ])
    if (sumr == 0) {
      drug.adjm[i, ] <- numeric(length = length(drug.adjm[i, ]))
    }
    if (sumr > 0) {
      drug.adjm[i, ] <- drug.adjm[i, ]/sum(drug.adjm[i, ])
    }
  }
  drug.adjm <- as.matrix(drug.adjm)
  drug.adjm <- t(drug.adjm)

  while (delta > 1e-10) {
    PT1 <- (1 - gamma) * drug.adjm
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma * p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
  }
  PT <- t(PT)
  rownames(PT) <- rownames(drug.adjm)
  return(PT)
}

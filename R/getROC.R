##' This function is used to plot ROC.
##'
##' @title getROC
##' @param drug.el A edge list of drugs, which is a two-column matrix, each row defines one edge. Numbers in the edge list represent NSC-ID of drugs.
##' @param p0 A vector of approved drugs' NSC-ID of interested cancer.
##' @param gamma gamma=0.7(default). A probability of losing when doing Random Walk. On the contray, there is a probability of 1-gamma left to itself. The range of this value is (0, 1).
##' @param filename filename = "ROC.pdf"(default). File name and path where to save the PDF. Filetype is decided by the extension in the path. Currently only .pdf formats are supported.
##' @return ROC
##' @import ROCR
##' @examples
##' e <- getData("drug.edgelist")
##' brc <- getData("breast_cancer")
##' \donttest{getROC(e, brc)}
##' @export
getROC <- function(drug.el, p0, gamma = 0.7, filename = "ROC.pdf") {
  drug.el <- apply(drug.el, 2, as.character)
  g <- graph_from_edgelist(drug.el, directed = F)
  iNet <- as_adj(g, type = "both", names = T, sparse = F)

  seed_p0 <- p0
  roc_input <- data.frame()
  all_p0 <- matrix(data = 0, ncol = 1, nrow = nrow(iNet))
  rownames(all_p0) <- row.names(iNet)

  pb <- txtProgressBar(min = 0, max = length(seed_p0), style = 3)

  for (i in 1:length(seed_p0)){
    all_p0[which(rownames(all_p0) %in% seed_p0),1] <- 1
    all_p0[which(rownames(all_p0)==seed_p0[i]),1] <- 0
    res <- rw(as.matrix(iNet),all_p0[,1],gamma)
    label <- data.frame(name=rownames(iNet),size=0)
    label[which(label$name==seed_p0[i]),2] <- 1
    temp_roc_input <- data.frame(score=res,label=label$size)
    roc_input <- rbind(roc_input,temp_roc_input)

    setTxtProgressBar(pb, i)
  }
  close(pb)

  pdf(filename,height=5,width=5)
  pred <- prediction(roc_input$score, roc_input$label)
  perf <- performance(pred,'tpr','fpr')
  auc <- performance(pred, "auc")@y.values[[1]]

  plot(perf,colorize=F,col="black",xlab="FPR",ylab="TPR",font.lab=1)
  abline(a=0,b=1,col="red",lty=3)
  legend("bottomright",paste("drrp (",round(auc,3),")",sep=""),lty=1)
  dev.off()
  print(paste("AUC is",auc))
}

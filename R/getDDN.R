##' This function is used to generate drug drug similarity network.
##'
##' @title getDDN
##' @param drug.el A edge list of drugs, which is a two-column matrix, each row defines one edge. Numbers in the edge list represent NSC-ID of drugs.
##' @param r.set A set of drugs that you used to prioritize candidates.
##' @param candidates A set of drugs that have been prioritized.
##' @param file file = "network.html"(default). File name and path where to save the HTML web page. Currently only .html formats are supported.
##'
##' @return A HTML web page within drug drug similarity network
##' @import igraph
##' @import visNetwork
##' @importFrom igraph degree
##' @examples
##' e <- getData("drug.edgelist")
##' brc <- getData("breast_cancer")
##' candidates <- getData("brc_candidates")
##' getDDN(e, brc, candidates)
##' @export
getDDN <- function(drug.el,r.set,candidates,file="network.html"){
  can <- drug.el[which((drug.el[,1] %in% candidates)|(drug.el[,2] %in% candidates)),]
  se <- drug.el[which((drug.el[,1] %in% r.set)|(drug.el[,2] %in% r.set)),]
  cs <- rbind(can,se)
  cs <- apply(cs,2,as.character)
  cs <- unique(cs)
  g <- graph_from_edgelist(cs,directed=F)

  deg <- degree(g)
  new.el <- as_edgelist(g)
  index <- deg[which(deg==1)]
  index <- as.data.frame(index)
  deg <- deg[-which(deg==1)]
  new.el <- new.el[-which((new.el[,1] %in% rownames(index)|(new.el[,2] %in% rownames(index)))),]
  g <- graph_from_edgelist(new.el,directed=F)

  d <- data.frame(V=V(g)$name)
  d$type <- 1
  d[which(d[,1] %in% candidates),2] <- 2
  d[which(d[,1] %in% r.set),2] <- 3

  nodes <- data.frame(id=V(g)$name)
  nodes$shape <- "circle"
  nodes$shadow <- FALSE
  nodes$label <- d$type
  nodes$color.background <- c("#cdd1d3","gold","tomato")[nodes$label]
  nodes$color.border <- "black"
  nodes$color.highlight.background <- "orange"
  nodes$color.highlight.border <- "darkred"
  #label
  new <- data.frame(id=nodes$id)
  new$nsc <- "NSC"
  new$nscid <- paste(new$nsc,new$id,sep="")
  nodes$label <- new$nscid
  #legend
  d <- data.frame(V=V(g)$name)
  d$type <- "Non-candidates"
  d[which(d[,1] %in% candidates),2] <- "Candidates"
  d[which(d[,1] %in% r.set),2] <- "Restart Set"
  nodes$group <- d$type


  links <- data.frame(from=new.el[,1],to=new.el[,2])
  links$color <- "gray"
  links$width <- 1.25


  visn <- visNetwork(nodes,links,width="1500px",height="750px") %>% visIgraphLayout(randomSeed=130,layout="layout_with_graphopt") #gem graphopt
  visn <- visGroups(visn,groupname="Non-candidates",shape="dot",color=list(background="#cdd1d3",border="black"))
  visn <- visGroups(visn,groupname="Candidates",shape="dot",color=list(background="gold",border="black"))
  visn <- visGroups(visn,groupname="Restart Set",shape="dot",color=list(background="tomato",border="black"))
  visLegend(visn,main=NULL,ncol=1,position="right",width=0.15) %>% visSave(file=file,selfcontained=F)
}

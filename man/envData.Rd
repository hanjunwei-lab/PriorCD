\name{envData}
\alias{envData}
\alias{drug.edgelist}
\alias{drug.info}
\alias{brc_candidates}
\alias{breast_cancer}
\alias{drug.fdr}
\alias{drug.ic50}
\alias{drug.r}
\alias{mRNA_path}
\alias{microRNA_path}
\alias{priorlist}
\docType{data}

\title{The variables in the environment include an example profile, a edgelist of our drug simiarity network, comprehensive drug information, restart drug set of breast cancer, candidate drugs of breast cancer, fdr of drug similarity network, correlation between drugs, mRNA and microRNA pathway activity profiles we've enriched.
}
\description{
Drug repurposing has become the focus of experts in drug development. In PriorCD, pathway activities and drug activities are combine to construct drug functional similarity network, and on which a global network propagation algorithm is applied. First, drug functional similarity network is constructed by the correlation and fdr of drug pairs. Then a global netwokr propagation (RWR) is performed on this network to prioritize candidates. Finally, ROC and network structure of the result can be browsed in PriorCD by getROC and getDDN functions.
}
\details{
The environment variable includes the variable \code{drug.edgelist}, \code{drug.info},\code{brc_candidates},\code{breast_cancer},\code{drug.fdr},\code{drug.r},\code{drug.ic50},\code{mRNA_path},\code{microRNA_path},\code{priorlist}
}

\format{An environment variable}
\author{Junwei Han\email{hanjunwei1981@163.com},Baotong Zheng\email{btzheng1116@163.com},Jieyi Di\email{dijy_0419@yeah.net}}

\keyword{data}

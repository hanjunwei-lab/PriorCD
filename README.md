# PriorCD
Prioritizing Cancer Drugs for Interested Cancer

> Prioritize candidate cancer drugs for drug repositioning based on the random walk with restart algorithm in a drug-drug functional similarity network. 1) We firstly constructed a drug-drug functional similarity network by integrating pathway activity and drug activity derived from the NCI-60 cancer cell lines. 2) Secondly, we calculated drug repurposing score according to a set of approved therapeutic drugs of interested cancer based on the random walk with restart algorithm in the drug-drug functional similarity network. 3) Finally, the permutation test was used to calculate the statistical significance level for the drug repurposing score.

### how to install

Install the **development version** from Github:

```R
Installation method：

1. library(devtools); 
   install_github("hanjunwei-lab/PriorCD")
2. install.packages("PriorCD")

Use：
library(PriorCD)
```

Please cite the following article when using `PriorCD`:

Di, J., B. Zheng, Q. Kong, Y. Jiang, S. Liu, Y. Yang, X. Han, Y. Sheng, Y. Zhang, L. Cheng, and J. Han, Prioritization of candidate cancer drugs based on a drug functional similarity network constructed by integrating pathway activities and drug activities. Mol Oncol, 2019. 13(10): p. 2259-2277.
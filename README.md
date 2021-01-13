# R_WGCNA
R scripts for RNAseq analysis unsing the WGCNA R package.

# Citing the WGCNA package
If you use WGCNA in published work, please cite it to properly credit people who have created it.

The WGCNA as an analysis method is described in
    Zhang B and Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834

The package implementation is described in the article
    Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559

If you use any q-value (FDR) calculations, please also cite at least one of the following articles:
    Storey JD. (2002) A direct approach to false discovery rates. Journal of the Royal Statistical Society, Series B, 64: 479-498.
    Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide experiments. Proceedings of the National Academy of Sciences, 100: 9440-9445.
    Storey JD, Taylor JE, and Siegmund D. (2004) Strong control, conservative point estimation, and simultaneous conservative consistency of false discovery rates: A unified approach. Journal of the Royal Statistical Society, Series B, 66: 187-205.

If you use the collapseRows function to summarize/convert probe-level data to gene-level data, please cite
    Miller JA, Cai C, Langfelder P, Geschwind DH, Kurian SM, Salomon DR, Horvath S (2011) Strategies for aggregating gene expression data: The collapseRows R function. BMC Bioinformatics 12:322

If you use module preservation calculations, please cite
    Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is my network module preserved and reproducible? PloS Comp Biol. 7(1): e1001057

If you use functions rgcolors.func, plotCor, plotMat, stat.bwss, or stat.diag.da, please also cite the article
    Sandrine Dudoit, Yee Hwa Yang, Matthew J. Callow, and Terence P. Speed, Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments, STATISTICA SINICA 2002, 12:111

# Acknowledgments
The core of the functions and other code was written by Peter Langfelder and Steve Horvath, partly based on older code written by Steve Horvath and Bin Zhang. Multiple people contributed additional code, most prominently Jeremy Miller, Chaochao (Ricky) Cai, Lin Song, Jun Dong, and Andy Yip. The package also contains code adapted from external packages that were either orphaned (such as package sma) or their development has made the code difficult to use in WGCNA (such as package qvalue). A big thanks goes out to people who continue report the many bugs in the package.

The package is currently maintained by Peter Langfelder.

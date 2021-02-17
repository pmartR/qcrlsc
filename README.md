# qcrlsc
Quality Control-based Robust LOESS Signal Correction

This package requires an omicsData object that can be created using the pmartR package. The sole functionality in this package is to perform normalization utilizing quality control (QC) samples, with the aim to reduce batch effects. The functionality is based on the quality control sample based robust LOESS (locally estimated scatterplot smoothing) signal correction (QC-RLSC) method described by Dunn \emph{et al}. (2011)

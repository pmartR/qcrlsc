# qcrlsc
Quality Control-based Robust LOESS Signal Correction

This package requires an omicsData object that can be created using the pmartR package. The sole functionality in this package is to perform normalization utilizing quality control (QC) samples, with the aim to reduce batch effects. The functionality is based on the quality control sample based robust LOESS (locally estimated scatterplot smoothing) signal correction (QC-RLSC) method described by Dunn et al. (2011).

Dunn,W.B., Broadhurst,D., Begley,P., Zelena,E., Francis-McIntyre,S., Anderson,N., Brown,M., Knowles,J.D., Halsall,A., Haselden,J.N. et al. (2011) Procedures for large-scale metabolic profiling of serum and plasma using gas chromatography and liquid chromatography coupled to mass spectrometry. Nat. Protoc., 6, 1060-1083

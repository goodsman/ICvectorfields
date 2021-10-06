---
title: "NEWS"
output: html_document
---

# ICvectorfields 0.0.1
* Accepted on CRAN (commit 6f85218)

# ICvectorfields 0.0.2
* Removed false fftw dependency (commit 487bd8de)
* Accepted by CRAN (commit 2d40b14)

# ICvectorfields 0.1.0
* Added PixelCt function which allows assessment of confidence placed in displacement estimates (commit a9cc8367)
* Added SubgridStats function for computing statistics on variables that may influence movement (commit 8450d7c7)
* Added SubgridMoransI function for computing Moran's I at sub-grid level (commit 4e499d96)
* Added DispStats function to compute shifted statistics at source or sink locations
* Added DispMoransI function to compute shifted Moran's I at source or sink locations
* Added PatternDetect function to classify vector field patterns as divergence or convergence
* Added RotationDetect function to classify vector field patterns as divergence or convergence
* Fixed errors caught by the majority of unit tests when terra was updated (commit 15ed0f27)

# ICvectorfields 0.1.1
* Updated for compatibility with the latest version of terra

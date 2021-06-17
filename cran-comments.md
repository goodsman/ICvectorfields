---
title: "cran-comments"
output: html_document
---

## Test environments
* Local windows 10 install, R 4.0.5

* Mac OSx using xcode12.5 (on travis-ci), R 4.1.0

* win-builder (devel and release), R 4.1.0

## R CMD check results

There were 0 ERRORs, 0 WARNINGs, and 1 NOTEs. 

There was 1 NOTE:

* Maintainer: 'Devin Goodsman <goodsman@ualberta.ca>'
  New submission

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Second resubmission 
This is a resubmission. I received a message indicating to shorten the title and add any methodological references to the description field of the DESCRIPTION file. To address these concerns I

* removed redundant section of title 

* hope to submit a manuscript describing the novel aspect of the algorithm referenced in the DESCRIPTION file      shortly in *The R Journal*. The Digital Image Correlation technique used in the package in general can not be       attributed to a single author or group of authors as it is a family of approaches and so I have not added  additional references. However, if this is still a concern for the reviewer, I will add a reference to a text book on the next revision (for example M. A. Sutton, J.-J. Orteu, H. W. Schreier, Book - Image Correlation for Shape, Motion and Deformation Measurements, Hardcover ISBN 978-0-387-78746-6)

## First resubmission 
This is a resubmission. The previous submission did not pass the incoming checks automatically. To mitigate the issues with failing the pre-submission tests I have:

* removed the + file LICENSE in the description that was inconsistent GPL (>= 3)

* Removed the for-loop from the exmaples of DispField and related functions so that they run faster

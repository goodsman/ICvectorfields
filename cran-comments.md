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

## Resubmission
This is a resubmission. The previous submission did not pass the incoming checks automatically. To mitigate the issues with failing the pre-submission tests I have:

* removed the + file LICENSE in the description that was inconsistent GPL (>= 3)

* Removed the for-loop from the exmaples of DispField and related functions so that they run faster

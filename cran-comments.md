---
title: "cran-comments"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Test environments
* Local windows 10 install R 4.0.5
*


## R CMD check results
There were 0 ERRORs, 1 WARNING and 0 NOTEs. 

There was 1 WARNING:

* checking data for ASCII and uncompressed saves ... OK
   WARNING
  'qpdf' is needed for checks on size reduction of PDFs
  
To address this WARNING I have used the following when checking
which did not resolve the issue.
```{r}
# devtools::check(build_args = c('--compact.vignettes=qpdf'))
```

When calling 

```{r}
# check_win_release()
```

I received this message, which may be related:

* Version contains large components (0.0.0.9000)

However, there were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:

* Ecoengine: this appears to be a failure related to config on 
  that machine. I couldn't reproduce it locally, and it doesn't 
  seem to be related to changes in httr (the same problem exists 
  with httr 0.4).

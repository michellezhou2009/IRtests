# IRtests

This R package is to implement the proposed information ratio (IR) test and in-and-out-of-sample pseudo (PIOS) likelihood ratio test in the manuscript titled "Information matrix equivalence in the presence of censoring: A goodness-of-fit test for semiparametric copula models with multivariate survival data".

## Instrall `IRtests`

To install this R package, you need to first install the `devtool` package via
```{r}
install.packages("devtools")
```
To install the `NCCIPW` package,
```{r}
library(devtools)
install_github("michellezhou2009/IRtests")
library(IRtests)
```

## Example

An example data `twins` is included in the package. This data is from the Australian NHMRC Twin Registry and available at https://genepi.qimr.edu.au/staff/davidD/Appendix/. The following gives the R code to test the goodness-of-fit of the Clayton copula for 567 monozygotic same-sex male twin pairs.
```{r}

```

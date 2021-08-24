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
load(app)
workdat = app[app$zyg==1,]
x1 = workdat$onset[workdat$id==1]
x2 = workdat$onset[workdat$id==2]
d1 =  workdat$app[workdat$id==1]
d2 =  workdat$app[workdat$id==2]
ncl = detectCores(logical = FALSE)
cl = makeCluster(ncl)
registerDoSNOW(cl)
out = IRtest_BiSurvCopula(x1=x1,x2=x2, d1=d1,d2=d2,copula.fam = "clayton",
                    control=list(yes.boot=TRUE,nboot=1000,same.cen=TRUE,seed1=20210823))
stopCluster(cl)
```


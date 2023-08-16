# IRtests

This R package is to implement the proposed information ratio (IR) test and in-and-out-of-sample pseudo (PIOS) likelihood ratio test in the manuscript titled "Information matrix equivalence in the presence of censoring: A goodness-of-fit test for semiparametric copula models with multivariate survival data".

## Instrall `IRtests`

To install this R package, you need to first install the `devtool` package via
```{r}
install.packages("devtools")
```
To install the `IRtests` package,
```{r}
library(devtools)
install_github("michellezhou2009/IRtests")
library(IRtests)
```

## Example

An example data `app` is included in the package. This data is from the Australian NHMRC Twin Registry and available at https://genepi.qimr.edu.au/staff/davidD/Appendix/. The following gives the R code to test the goodness-of-fit of the Clayton copula for 1231 monozygotic same-sex female twin pairs.

```{r}
data(app)
workdat = app[app$zyg == 3,]
x1 = workdat$onset[workdat$id==1]
x2 = workdat$onset[workdat$id==2]
d1 =  workdat$app[workdat$id==1]
d2 =  workdat$app[workdat$id==2]
cl = makeCluster(4)
registerDoSNOW(cl)
copula.fam.all = c("clayton", "frank", "gumbel", "gaussian")
res.all = lapply(copula.fam.all, function(copula.fam){
  out = IRtest_BiSurvCopula(x1 = x1, x2 = x2, d1 = d1, d2 = d2, 
                            copula.fam = copula.fam,
                            control=list(yes.boot = TRUE, nboot = 1000, 
                                         same.cen = TRUE, 
                                         seed1 = 20210823))
  data.frame(
    family = copula.fam, theta = out$theta_est, IR = out$IR, pval = out$pval
  )
}) 
do.call(rbind, res.all)
stopCluster(cl)
```



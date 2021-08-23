#' IR test for misspecification of copula functions in semi-parametric survival copula model for bivariate right-censored data
#'
#' @description
#'
#' @param x1 a vector, the first response
#' @param x2 a vector, the second response
#' @param d1 a vector of indicators whether each observation in \code{x1} is fully observed: \code{1} indicates the observation is fully observed, and \code{0} indicates the observation is censored
#' @param d2 a vector of indicators whether each observation in \code{x2} is fully observed: \code{1} indicates the observation is fully observed, and \code{0} indicates the observation is censored
#' @param copula.fam a character indicating which one of the following copula families: "clayton", "frank", "joe", "gumbel", and "normal"
#' @param control a list of the following components: \code{yes.boot}, \code{nboot}, \code{seed1}, and \code{same.cen}. \code{yes.boot} is a logical value indicating whether to implement the bootstrap procedure. \code{nboot} is the number of bootstrap samples.  \code{seed1} is the seed for generating the bootstrap samples. \code{same.cen} is a logical value indicating whether the censoring time is same for both event time.
#'
#' @import pbivnorm copula survival numDeriv foreach parallel doSNOW
#'
#' @export
#'
#' @return a list of the following components: \code{theta.est}, the PMLE of the copula parameter,  \code{IR}, IR test statistic, \code{IR.boot}, the bootstrapped resamples of the IR test statistics if \code{yes.boot=TRUE}, and \code{pval}, the one-sided and two-sided p-values if \code{yes.boot=TRUE}.
#'
#'

IRtest_BiSurvCopula <- function(x1,x2,d1,d2, copula.fam, control=list(yes.boot=TRUE, nboot=500, seed1=1234, same.cen=FALSE)){

  # checking d1/d2 has other values than 0/1
  # checking whether x1, x2, d1, d2 any missing data
  # checking dimensions
  # checking whether copula.fam is not included in the R package
  if (length(unique(c(length(x1),length(x2),length(d1),length(d2))))>1) stop("The length of \"x1\", \"x2\", \"d1\", and \"d2\" has to be same.")
  if (any(is.na(d1))) stop("Missing values in \"d1\".")
  if (any(is.na(d2))) stop("Missing values in \"d2\".")
  if (any(is.na(x1))) stop("Missing values in \"x1\".")
  if (any(is.na(x2))) stop("Missing values in \"x2\".")
  if (length(setdiff(unique(d1),c(0,1)))>0) stop("Values in \"d1\" can only be 0 or 1.")
  if (length(setdiff(unique(d2),c(0,1)))>0) stop("Values in \"d2\" can only be 0 or 1.")
  if(!(copula.fam%in% c("joe","gumbel","clayton","frank","normal"))){
    stop("'copula.fam' should be one of `joe`, 'gumbel','clayton','frank','normal'")
  }

  seed1 = control$seed1
  yes.boot = control$yes.boot
  nboot = control$nboot
  same.cen=control$same.cen

  nsamp = length(x1)
  update.out = update.data(x1,x2,d1,d2,type="surv")
  distfun.t1 = update.out$distfun.t1; distfun.t2 = update.out$distfun.t2

  out = IR.fun(u1=update.out$u1,u2=update.out$u2,d1,d2,copula.fam)
  theta_est = out$theta.est
  IR = out$IR

  if (yes.boot) {
    if (same.cen){
      distfun.c = survfit(Surv(pmax(x1,x2),1*(d1==0 | d2==0)) ~ 1,se.fit=F,type="fh2")
    } else {
      distfun.c1 = survfit(Surv(x1,1-d1) ~ 1,se.fit=F,type="fh2")
      distfun.c2 = survfit(Surv(x2,1-d2) ~ 1,se.fit=F,type="fh2")
    }
    set.seed(seed1)
    seeds_b = round(runif(nboot,10,100000))
    boot_out = foreach (b=1:nboot,
                        .packages=c("survival","copula","numDeriv","IRtests"),
                        .combine = c) %dopar% {
                          set.seed(seeds_b[b])
                          # generate the bivariate event times (T1,T2)
                          switch(copula.fam,
                                 "joe"={cc_b=joeCopula(theta_est)},
                                 "clayton"={cc_b=claytonCopula(theta_est)},
                                 "frank"={cc_b=frankCopula(theta_est)},
                                 "gumbel"={cc_b=gumbelCopula(theta_est)},
                                 "normal"={cc_b=normalCopula(theta_est)}
                          )
                          uu_b = rCopula(nsamp,cc_b)
                          t1_b = inv_surv(uu_b[,1],distfun.t1)
                          t2_b = inv_surv(uu_b[,2],distfun.t2)
                          # generate the censoring through the inverse of its empirical survival
                          if (same.cen) c1_b = c2_b = inv_surv(runif(nrow(uu_b)),distfun.c) else {
                            c1_b = inv_surv(runif(nrow(uu_b)),distfun.c1)
                            c2_b = inv_surv(runif(nrow(uu_b)),distfun.c2)
                          }
                          # generate bootstrap data
                          x1_b = pmin(t1_b,c1_b); x2_b = pmin(t2_b,c2_b)
                          d1_b = ifelse(t1_b<=c1_b,1,0); d2_b=ifelse(t2_b<=c2_b,1,0)
                          update.b = update.data(x1_b,x2_b,d1_b,d2_b)
                          out_sp = IR.fun(u1=update.b$u1,u2=update.b$u2,d1=d1_b,d2=d2_b,copula.fam)
                          out_sp$IR
                        }
    boot.u = boot_out - mean(boot_out,na.rm=T)
    if (!is.na(IR)) pval = c("onesided" = mean(boot.u>(IR-1),na.rm=T), "twosided" = mean(boot.u-abs(IR-1)>0,na.rm=T)+mean(boot.u+abs(IR-1)<0,na.rm=T)) else pval = c("onesided"=0,"twosided"=0)
  } else {boot_out = NA;pval=NA}
  return(list("theta_est"=theta_est,"IR"=IR,"IR.boot"=boot_out,"pval"=pval))
}

#' PIOS and IR tests of copula function for bivariate survival time
#'
#' @description
#'
#' @param data a data frame with four columns, in which the first two columns are censored event times, and the next two columns are censoring indicator variables
#' @param copula.fam a character indicating which one of the following copula families: "clayton", "frank", "gumbel", and "normal"
#' @param control, a list of the following components: \code{yes.PIOSexact}, \cdoe{yes.boot}, \code{nboot}, \code{seed1}. \code{yes.PIOSexact} is a logical value indicating whether to calculate the exact PIOS test statistic, and the default value is \code{FALSE}. \code{yes.boot} is a logical value indicating whether to implement the bootstrap procedure. \code{nboot} is the number of bootstrap samples.  \code{seed1} is the seed for generating the bootstrap samples.
#'
#' @import copula survival numDeriv foreach parallel doSNOW
#'
#' @export
#'
#' @return
#'
#'

GoFtests_BiSurv <- function(data, copula.fam, control=list(yes.PIOSexact = FALSE, yes.boot = TRUE, nboot=500, seed1=1234)){

  if(!(copula.fam%in% c("gumbel","clayton","frank","normal"))){
    stop("'copula.fam' should be one of 'gumbel','clayton','frank','normal'")
  }

  seed1 = control$seed1
  yes.boot = control$yes.boot
  nboot = control$nboot
  yes.PIOSexact = control$yes.PIOSexact

  nsamp = nrow(data)
  update.out = update.data(data)
  data = update.out$data
  surv.t1 = update.out$surv.t1; surv.t2 = update.out$surv.t2

  out = lr.test(data[,-c(1,2)],copula.fam,yes.PIOSexact=yes.PIOSexact)
  theta_est = out$theta.est
  PIOS.approx = out$PIOS.approx
  PIOS.exact = out$PIOS.exact
  IR = out$IR

  if (yes.boot) {
    surv.c = survival::survfit(survival::Surv(pmax(data[,1],data[,2]),1*(data[,3]==0 | data[,4]==0)) ~ 1)
    set.seed(seed1)
    seeds_b = round(runif(nboot,10,100000))
    library(parallel)
    library(doSNOW)
    library(foreach)
    ncl = detectCores(logical = FALSE)
    cl = makeCluster(ncl)
    registerDoSNOW(cl)

    boot_out = foreach (b=1:nboot,
                         .packages=c("survival","copula","numDeriv","IRtools"),
                         .combine = cbind) %dopar% {
      set.seed(seeds_b[b])
      switch(copula.fam,
             "clayton"={cc_b=copula::claytonCopula(theta_est)},
             "frank"={cc_b=copula::frankCopula(theta_est)},
             "gumbel"={cc_b=copula::gumbelCopula(theta_est)},
             "normal"={cc_b=copula::normalCopula(theta_est)}
      )
      uu_b = copula::rCopula(nsamp,cc_b)
      # generate the event times through the inverse of their respective empirical survival estimates
      t1_b = inv_surv(uu_b[,1],surv.t1)
      t2_b = inv_surv(uu_b[,2],surv.t2)
      # generate the censoring through the inverse of its empirical survival
      c_b =  inv_surv(runif(nrow(uu_b)),surv.c)
      # generate bootstrap data
      data_sp = data.frame(x1=pmin(t1_b,c_b),x2=pmin(t2_b,c_b),d1=ifelse(t1_b<=c_b,1,0),d2=ifelse(t2_b<=c_b,1,0))
      data_sp = update.data(data_sp)$data
      out_sp = lr.test(data_sp[,-c(1,2)],copula.fam,yes.PIOSexact=FALSE)
      c("sp_PIOS"=out_sp$PIOS.approx,"sp_IR"=out_sp$IR)
                         }
    stopCluster(cl)
  } else {boot_out = NA}
  return(list("theta_est"=theta_est,"PIOS"=PIOS.approx,"IR"=IR,"boot"=boot_out))
}

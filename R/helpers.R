update.data = function(x1,x2,d1,d2,type="surv"){
  if (type=="surv"){
    distfun.t1= survfit(Surv(x1,d1)~1,se.fit=F,type="fh2")
    u1 = summary(distfun.t1,sort(x1))$surv[rank(x1)]
    distfun.t2= survfit(Surv(x2,d2)~1,se.fit=F,type="fh2")
    u2 = summary(distfun.t2,sort(x2))$surv[rank(x2)]
    return(list(u1=u1,u2=u2,distfun.t1=distfun.t1,distfun.t2=distfun.t2))
  }
}

inv_surv <- function(uu,distfun){
  out = quantile(distfun,1-uu)
  out[is.na(out)] = max(distfun$time)
  out
}

inv_cdf <- function(uu,distfun){
  quantile(distfun,uu)
}

# Bivariate copula functions with derivative
CopulaFuns <- function(copula.family){
  switch(copula.family,
         "gumbel" = {
           gumbel = expression(exp(-((-log(u1))^(theta) + (-log(u2))^(theta))^(1/theta)))
           dg.u1 = function(u1,u2,theta){}; body(dg.u1) = D(gumbel, "u1") # first-order derivative of gumbel with respect to u1
           dg.u2 = function(u1,u2,theta){}; body(dg.u2) = D(gumbel, "u2") # first-order derivative of gumbel with respect to u2
           dg.u1u2 = function(u1,u2,theta){}; body(dg.u1u2) = D(D(gumbel, "u1"), "u2") # 2nd-order derivative of gumbel with respect to u1 and u2
           C.fun = function(u1,u2,theta){}; body(C.fun) = gumbel
           list(C=C.fun, c.u2=dg.u2, c.u1=dg.u1,c.u1u2=dg.u1u2)
         },
         "clayton" = {
           clayton = expression((u1^(-theta) + u2^(-theta) -1)^(-1/theta))
           dc.u1 = function(u1,u2,theta){}; body(dc.u1) = D(clayton, "u1")
           dc.u2 = function(u1,u2,theta){}; body(dc.u2) = D(clayton, "u2")
           dc.u1u2 = function(u1,u2,theta){}; body(dc.u1u2) = D(D(clayton, "u1"), "u2")
           C.fun = function(u1,u2,theta){}; body(C.fun) = clayton
           list(C=C.fun, c.u2=dc.u2, c.u1=dc.u1,c.u1u2=dc.u1u2)
         },
         "frank"= {
           frank = expression((-1/theta)*log(1 + (exp(-theta*u1) - 1)*(exp(-theta*u2) - 1)/(exp(-theta) - 1)))
           df.u1 = function(u1,u2,theta){}; body(df.u1) = D(frank, "u1")
           df.u2 = function(u1,u2,theta){}; body(df.u2) = D(frank, "u2")
           df.u1u2 = function(u1,u2,theta){}; body(df.u1u2) = D(D(frank, "u1"), "u2")
           C.fun = function(u1,u2,theta){}; body(C.fun) = frank
           list(C=C.fun, c.u2=df.u2, c.u1=df.u1,c.u1u2=df.u1u2)
         },
         "normal"={
           dn.u1 = function(u1,u2,theta){
             pnorm(qnorm(u2),mean=theta*qnorm(u1),sd=sqrt(1-theta^2))*dnorm(qnorm(u1))
           }
           dn.u2 = function(u1,u2,theta){
             pnorm(qnorm(u1),mean=theta*qnorm(u2),sd=sqrt(1-theta^2))*dnorm(qnorm(u2))
           }
           dn.u1u2 = function(u1,u2,theta){
             dnorm(qnorm(u1),mean=theta*qnorm(u2),sd=sqrt(1-theta^2))*dnorm(qnorm(u2))
           }
           C.fun = function(u1,u2,theta){
             pbivnorm::pbivnorm(x=qnorm(cbind(u1,u2)),rho=theta)
           }
           list(C=C.fun, c.u2=dn.u2, c.u1=dn.u1,c.u1u2=dn.u1u2)
         }
         )

  }

# log-likelihood calculation for one observation
nlik1 = function(theta,u1,u2,d1,d2,copula.fam){
  cop.fun <- CopulaFuns(copula.fam)
  code = 2*d1 + d2 + 1
  code.ls = sort(unique(code))
  index = unlist(lapply(code.ls,function(c){which(code==c)}))
  out.ls = unlist(sapply(code.ls,function(c){
    cop.fun[[c]](u1[code==c],u2[code==c],theta)
  }))
  out = code*0
  out[index] = out.ls
  #return(-sum(log(out),na.rm=T))
  return(-log(out))
}

# log likelihood calculation for more than one observations
nlik = function(theta, u1,u2,d1,d2, copula.fam){
  cop.fun <- CopulaFuns(copula.fam)
  code = 2*d1 + d2 + 1
  code.ls = sort(unique(code))
  out = unlist(sapply(code.ls,function(c){
    cop.fun[[c]](u1[code==c],u2[code==c],theta)
  }))
  return(-sum(log(out),na.rm=T))
}

#' Calculation of PIOS test statistic
#'
#' @description Given the data and copula family, calculate the PIOS test statistic.
#'
#' @param u1 a vector, the first pseudo-observations
#' @param u2 a vector, the second pseudo-observations
#' @param d1 a vector of indicators whether each observation in the first response is fully observed: \code{1} indicates the observation is fully observed, and \code{0} indicates the observation is censored
#' @param d2 a vector of indicators whether each observation in the second response is fully observed: \code{1} indicates the observation is fully observed, and \code{0} indicates the observation is censored
#' @param copula.fam a character indicating which one of the following copula families: "clayton", "frank", "gumbel", and "normal"
#' @param yes.exact a logical value indicating whether to calculate the exact test statistic; if \code{yes.exact=FALSE} (default value), the approximate test statistic is calculated.

#' @import pbivnorm copula numDeriv foreach parallel doSNOW
#'
#' @export
#'
#' @return a list of the following components: \code{theta.est}, the PMLE of the copula parameter, and \code{PIOS}, PIOS exact test statistic if \code{yes.exact=TRUE} or PIOS approximate test statistic if \code{yes.exact=FALSE} (default).
#'
#'
PIOS.fun = function(u1,u2,d1,d2,copula.fam,yes.exact=F){

  if (length(unique(c(length(u1),length(u2),length(d1),length(d2))))>1) stop("The length of \"x1\", \"x2\", \"d1\", and \"d2\" has to be same.")
  if (any(is.na(d1))) stop("Missing values in \"d1\".")
  if (any(is.na(d2))) stop("Missing values in \"d2\".")
  if (any(is.na(u1))) stop("Missing values in \"u1\".")
  if (any(is.na(u2))) stop("Missing values in \"u2\".")
  if (length(setdiff(unique(d1),c(0,1)))>0) stop("Values in \"d1\" can only be 0 or 1.")
  if (length(setdiff(unique(d2),c(0,1)))>0) stop("Values in \"d2\" can only be 0 or 1.")
  if(!(copula.fam%in% c("gumbel","clayton","frank","normal"))){
    stop("'copula.fam' should be one of 'gumbel','clayton','frank','normal'")
  }

  # obtain the initial value of theta using the relationship between copula parameter and kendall's tau
  switch(copula.fam,
         "clayton"={
           theta.ini = copClayton@iTau(cor(u1, u2, method = "kendall"))
         },
         "frank"={
           theta.ini = copFrank@iTau(cor(u1, u2, method = "kendall"))
         },
         "gumbel"={
           theta.ini = copGumbel@iTau(cor(u1, u2, method = "kendall"))
         },
         "normal" = {
           theta.ini = cor(u1, u2)
         }
  )

  # define boundary value of optimization
  low = c(1,rep(0,3)); up = c(rep(40,3),1);
  names(low) = names(up) = c("gumbel","clayton","frank","normal")

  pseudoMLE.out = optim(theta.ini, nlik, u1=u1,u2=u2,d1=d1,d2=d2, copula.fam = copula.fam, method = "Brent", lower = low[copula.fam], upper = up[copula.fam],hessian=T)
  theta.est = pseudoMLE.out$par

  # in-sample pseudo log-likelihood
  is.lik = pseudoMLE.out$value
  if (yes.exact){
    os.lik = foreach (k=1:length(u1),.packages=c("copula","IRtests"),.combine = c) %dopar% {
      os.mle = optim(theta.ini, nlik, copula.fam = copula.fam,u1=u1[-k],u2=u2[-k],d1=d1[-k],d2=d2[-k], method = "Brent", lower = low[copula.fam], upper = up[copula.fam])$par
      nlik1(theta=os.mle,u1=u1[k],u2=u2[k],d1=d1[k],d2=d2[k], copula.fam = copula.fam)
    }
  } else {
    ii = pseudoMLE.out$hessian
    if (!is.na(ii)) ii_inv = try(solve(ii),silent=T) else ii_inv=NA
    if (!is.character(ii_inv) & !is.na(ii_inv)) {
      ss = numDeriv::jacobian(nlik1,x=theta.est,u1=u1,u2=u2,d1=d1,d2=d2,copula.fam=copula.fam)
      # ss = do.call(rbind,lapply(1:length(u1),function(k){
      #   numDeriv::jacobian(nlik1,x=theta.est,u1=u1[k],u2=u2[k],d1=d1[k],d2=d2[k],copula.fam=copula.fam)}
      # ))
      os.mle = as.vector(theta.est + ii_inv%*%t(ss))
      os.lik = sapply(1:length(os.mle),function(i){
        if (!is.na(os.mle[i])) {
          nlik1(os.mle[i],u1=u1[i],u2=u2[i],d1=d1[i],d2=d2[i],copula.fam=copula.fam)
        } else return(NA)

      })
    } else {
    os.lik=NULL
    warning("The sensitivity matrix cannot be calucated.")
    }
  }
  return(list(theta.est = theta.est, PIOS=sum(os.lik,na.rm=T) - is.lik))

}

#' Calculation of IR test statistic
#'
#' @description Given the data and copula family, calculate the IR test statistic
#'
#' @param u1 a vector, the first pseudo-observations
#' @param u2 a vector, the second pseudo-observations
#' @param d1 a vector of indicators whether each observation in the first response is fully observed: \code{1} indicates the observation is fully observed, and \code{0} indicates the observation is censored
#' @param d2 a vector of indicators whether each observation in the second response is fully observed: \code{1} indicates the observation is fully observed, and \code{0} indicates the observation is censored
#' @param copula.fam a character indicating which one of the following copula families: "clayton", "frank", "gumbel", and "normal"

#' @import pbivnorm copula numDeriv
#'
#' @export
#'
#' @return a list of the following components: \code{theta.est}, the PMLE of the copula parameter, and \code{IR}, IR test statistic.
#'
#'
IR.fun = function(u1,u2,d1,d2, copula.fam){

  if (length(unique(c(length(u1),length(u2),length(d1),length(d2))))>1) stop("The length of \"x1\", \"x2\", \"d1\", and \"d2\" has to be same.")
  if (any(is.na(d1))) stop("Missing values in \"d1\".")
  if (any(is.na(d2))) stop("Missing values in \"d2\".")
  if (any(is.na(u1))) stop("Missing values in \"u1\".")
  if (any(is.na(u2))) stop("Missing values in \"u2\".")
  if (length(setdiff(unique(d1),c(0,1)))>0) stop("Values in \"d1\" can only be 0 or 1.")
  if (length(setdiff(unique(d2),c(0,1)))>0) stop("Values in \"d2\" can only be 0 or 1.")
  if(!(copula.fam%in% c("gumbel","clayton","frank","normal"))){
    stop("'copula.fam' should be one of 'gumbel','clayton','frank','normal'")
  }

  # obtain the initial value of theta using the relationship between copula parameter and kendall's tau
  switch(copula.fam,
         "clayton"={
           theta.ini = copClayton@iTau(cor(u1, u2, method = "kendall"))
         },
         "frank"={
           theta.ini = copFrank@iTau(cor(u1, u2, method = "kendall"))
         },
         "gumbel"={
           theta.ini = copGumbel@iTau(cor(u1, u2, method = "kendall"))
         },
         "normal" = {
           theta.ini = cor(u1, u2)
         }
  )

  # define boundary value of optimization
  low = c(1,rep(0,3)); up = c(rep(40,3),1);
  names(low) = names(up) = c("gumbel","clayton","frank","normal")

  pseudoMLE.out = optim(theta.ini, nlik, u1=u1,u2=u2,d1=d1,d2=d2, copula.fam = copula.fam, method = "Brent", lower = low[copula.fam], upper = up[copula.fam], hessian=T)
  theta.est = pseudoMLE.out$par
  ii = pseudoMLE.out$hessian
  if (!is.na(ii)) ii_inv = try(solve(ii),silent=T) else ii_inv=NA
  # ii = try(numDeriv::hessian(nlik,x=theta.est,u1=u1,u2=u2,d1=d1,d2=d2,copula.fam=copula.fam))
  # ii_inv = try(solve(ii),silent=T)

  if (!is.character(ii_inv) & !is.na(ii_inv)) {
    ss = numDeriv::jacobian(nlik1,x=theta.est,u1=u1,u2=u2,d1=d1,d2=d2,copula.fam=copula.fam)
    ss = ss[complete.cases(ss),,drop=F]
    # ss = do.call(rbind,lapply(1:length(u1),function(k){
    #   numDeriv::jacobian(nlik1,x=theta.est,u1=u1[k],u2=u2[k],d1=d1[k],d2=d2[k],copula.fam=copula.fam)}
    # ))
    V = t(ss)%*%ss; IR.stat = sum(diag(ii_inv%*%V))
  } else {
    IR.stat = NA
    warning("The sensitivity matrix cannot be calucated.")
  }
  return(list(theta.est = theta.est, IR=IR.stat))
}


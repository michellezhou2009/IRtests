inv_surv <- function(uu,surv.out){
  unlist(sapply(uu,function(u){
    index = which(surv.out$surv<u)[1]
    if(is.na(index)) max(surv.out$time) else
      if (index==1) min(surv.out$time) else surv.out$time[index-1]
  },simplify = T))
}

# Bivariate copula functions with derivative
CopulaFuns <- function(copula.family){
  switch(copula.family,
         "gumbel" = {
           gumbel = expression(exp(-((-log(s1))^(theta) + (-log(s2))^(theta))^(1/theta)))
           g1.der = function(s1,s2,theta){}; body(g1.der) = D(gumbel, "s1") # first-order derivative of gumbel with respect to s1
           g2.der = function(s1,s2,theta){}; body(g2.der) = D(gumbel, "s2") # first-order derivative of gumbel with respect to s2
           g3.der = function(s1,s2,theta){}; body(g3.der) = D(D(gumbel, "s1"), "s2") # 2nd-order derivative of gumbel with respect to s1 and s2
           g4.der = function(s1,s2,theta){}; body(g4.der) = gumbel
           list(C=g4.der, C.u2=g2.der, C.u1=g1.der,C.u1u2=g3.der)
         },
         "clayton" = {
           clayton = expression((s1^(-theta) + s2^(-theta) -1)^(-1/theta))
           c1.der = function(s1,s2,theta){}; body(c1.der) = D(clayton, "s1")
           c2.der = function(s1,s2,theta){}; body(c2.der) = D(clayton, "s2")
           c3.der = function(s1,s2,theta){}; body(c3.der) = D(D(clayton, "s1"), "s2")
           c4.der = function(s1,s2,theta){}; body(c4.der) = clayton
           list(C=c4.der, C.u2=c2.der, C.u1=c1.der,C.u1u2=c3.der)
         },
         "frank"= {
           frank = expression((-1/theta)*log(1 + (exp(-theta*s1) - 1)*(exp(-theta*s2) - 1)/(exp(-theta) - 1)))
           f1.der = function(s1,s2,theta){}; body(f1.der) = D(frank, "s1")
           f2.der = function(s1,s2,theta){}; body(f2.der) = D(frank, "s2")
           f3.der = function(s1,s2,theta){}; body(f3.der) = D(D(frank, "s1"), "s2")
           f4.der = function(s1,s2,theta){}; body(f4.der) = frank
           list(C=f4.der, C.u2=f2.der, C.u1=f1.der,C.u1u2=f3.der)
         },
         "normal"={
           n1.der = function(s1,s2,theta){
             pnorm(qnorm(s2),mean=theta*qnorm(s1),sd=sqrt(1-theta^2))*dnorm(qnorm(s1))
           }
           n2.der = function(s1,s2,theta){
             pnorm(qnorm(s1),mean=theta*qnorm(s2),sd=sqrt(1-theta^2))*dnorm(qnorm(s2))
           }
           n3.der = function(s1,s2,theta){
             dnorm(qnorm(s1),mean=theta*qnorm(s2),sd=sqrt(1-theta^2))*dnorm(qnorm(s2))
           }
           n4.der = function(s1,s2,theta){
             sigma = matrix(c(1,theta,theta,1),nrow=2,ncol=2)
             apply(cbind(s1,s2),1,function(s) mvtnorm::pmvnorm(upper=qnorm(s),lower=-Inf,sigma=sigma))
           }
           # n1.der = function(s1,s2,theta){
           #   sapply(1:length(s1),function(k){
           #     pnorm(qnorm(s2[k]),mean=theta*qnorm(s1[k]),sd=sqrt(1-theta^2))*dnorm(qnorm(s1[k]))
           #   })
           # }
           # n2.der = function(s1,s2,theta){
           #   sapply(1:length(s1),function(k){
           #     pnorm(qnorm(s1[k]),mean=theta*qnorm(s2[k]),sd=sqrt(1-theta^2))*dnorm(qnorm(s2[k]))
           #   })
           #   }
           # n3.der = function(s1,s2,theta){
           #   sapply(1:length(s1),function(k){
           #      dnorm(qnorm(s1[k]),mean=theta*qnorm(s2[k]),sd=sqrt(1-theta^2))*dnorm(qnorm(s2[k]))
           #   })
           #   }
           # n4.der = function(s1,s2,theta){
           #   sapply(1:length(s1),function(k){
           #     mvtnorm::pmvnorm(upper=qnorm(c(s1[k],s2[k])),lower=-Inf,sigma=matrix(c(1,theta,theta,1),nrow=2,ncol=2))
           #   })
           # }
           list(C=n4.der, C.u2=n2.der, C.u1=n1.der,C.u1u2=n3.der)
         }
         )

  }


nlik1 = function(theta, obs, copula.fam){
  cop.fun <- CopulaFuns(copula.fam)
  obs = as.numeric(obs)
  code = 2*obs[1] + obs[2] + 1
  out = cop.fun[[code]](obs[3],obs[4],theta)
  return(-log(out))
}

nlik = function(theta, obs, copula.fam){
  cop.fun <- CopulaFuns(copula.fam)

  code = 2*obs[,1] + obs[,2] + 1
  code.ls = sort(unique(code))
  out = unlist(sapply(code.ls,function(c){
    cop.fun[[c]](obs[code==c,3],obs[code==c,4],theta)
  }))
  return(-sum(log(out),na.rm=T))
}

#' PIOS and IR test statistic
#'
#' @description Given the data and copula family, calculate the PIOS and IR test.
#'
#' @param data a data frame or a matrix with four columns, in which the first two columns are the censoring indicator variables, and the last two columns are the estimated working variables.
#'
#' @param copula.fam a character indicating which one of the following copula families: "clayton", "frank", "gumbel", and "normal".
#'
#' @param yes.PIOSexact a logical value indicating whether to calculate the exact PIOS test statistic, and the default value is \code{FALSE}.
#'
#' @import copula numDeriv
#'
#' @export
#'
#' @return a list of the following components: \code{theta.est} (the PMLE of the copula parameter), \code{PIOS.exact} (exact PIOS test statistic), \code{PIOS.approx} (approximate PIOS test statistic), and \code{IR} (IR test statistic).
#'
#'
lr.test = function(data, copula.fam, yes.PIOSexact=F){

  if(!is.matrix(data) && !is.data.frame(data)){
    stop("'data' must be either matrix or data.frame")
  }
  if(!(copula.fam%in% c("gumbel","clayton","frank","normal"))){
    stop("'copula.fam' should be one of 'gumbel','clayton','frank','normal'")
  }

  # obtain the initial value of theta using the relationship between copula parameter and kendall's tau
  switch(copula.fam,
         "clayton"={
           theta.ini = copula::copClayton@iTau(cor(data[,3], data[,4], method = "kendall"))
         },
         "frank"={
           theta.ini = copula::copFrank@iTau(cor(data[,3], data[,4], method = "kendall"))
         },
         "gumbel"={
           theta.ini = copula::copGumbel@iTau(cor(data[,3], data[,4], method = "kendall"))
         },
         "normal" = {
           theta.ini = cor(data[,3], data[,4])
         }
         )

  # define boundary value of optimization
  low = c(1,rep(0,3)); up = c(rep(40,3),1);
  names(low) = names(up) = c("gumbel","clayton","frank","normal")

  pseudoMLE.out = optim(theta.ini, nlik, obs = data, copula.fam = copula.fam, method = "Brent", lower = low[copula.fam], upper = up[copula.fam])
  # in-sample pseudo log-likelihood
  is.lik = pseudoMLE.out$value
  theta.est = pseudoMLE.out$par

  ii = numDeriv::hessian(nlik,x=theta.est,obs=data,copula.fam=copula.fam)
  ii_inv = try(solve(ii),silent=T)

  if (!is.character(ii_inv) & !is.na(ii_inv)) {
    ss = sapply(1:nrow(data),function(k){
      numDeriv::jacobian(nlik1,x=theta.est,obs=data[k,],copula.fam=copula.fam)}
    )
    V = sum(ss^2); IR.stat = sum(diag(ii_inv%*%V))
    os.approx.all = as.vector(theta.est + ii_inv%*%ss)
    os.lik.approx = sapply(1:nrow(data),function(i){
      nlik1(os.approx.all[i],data[i,],copula.fam=copula.fam)
    })
    test.stat.approx=sum(os.lik.approx,na.rm = T) - is.lik
  } else {
    test.stat.approx = NA
    IR.stat = NA
    warning("The sensitivity matrix cannot be calucated.")
  }
  if (yes.PIOSexact){
    # out-of-sample pseudo log-likelihood
    library(parallel)
    library(doSNOW)
    library(foreach)
    ncl = detectCores(logical = FALSE)
    cl = makeCluster(ncl)
    registerDoSNOW(cl)
    os.lik = foreach (k=1:nrow(data),.packages=c("copula","IRtests"),.combine = c) %dopar% {
      os.mle = optim(theta.ini, nlik, copula.fam = copula.fam, obs = data[-k,,drop=F], method = "Brent", lower = low[copula.fam], upper = up[copula.fam])$par
      nlik1(theta = os.mle, obs = data[k,], copula.fam = copula.fam)
    }
    stopCluster(cl)
    test.stat.exact  = sum(os.lik,na.rm=T) - is.lik
  } else {
    test.stat.exact = NA
  }
  return(list(theta.est = theta.est, PIOS.exact=test.stat.exact, PIOS.approx=test.stat.approx, IR=IR.stat))
}

update.data = function(data){
  # calculate nelson estimate of survival function
  x1 = data[,1]; x2 = data[,2]; d1 = data[,3]; d2 = data[,4]

  surv.t1 = survival::survfit(survival::Surv(x1,d1)~1,se.fit=F,type="fh2")
  S.x1 = summary(surv.t1,sort(x1))$surv[rank(x1)]
  surv.t2 = survival::survfit(survival::Surv(x2,d2)~1,se.fit=F,type="fh2")
  S.x2 = summary(surv.t2,sort(x2))$surv[rank(x2)]

  return(list(data=cbind(data,S.x1,S.x2),surv.t1=surv.t1,surv.t2=surv.t2))
}


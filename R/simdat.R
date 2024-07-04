#' Add together two numbers
#'
#' @param x A number
#' @param y A number
#' @return The sum of \code{x} and \code{y}
#' @import refund
#' @import splines 
#' @import mgcv 
#' @import matrixcalc 
#' @import Matrix 
#' @import mvtnorm 
#' @import refund 
#' @import stats 
#' @import dplyr 
#' @import reshape2 
#' @import truncdist 
#' @import EnvStats
#' @import triangle
#' 
#' 
#' @export

# Simulate informative censoring data with covariate
# param: a list of setting values:
#  grid: grid of time point in [0,1]
#  n: number of subjects
#  rho: frailty parameter (0-1)
#  dist: distribution of frailty variable ('lognormal' or 'gamma')
#  C: censoring time (2: negatively correlatve with Zstar, 1: positively correlatve with Zstar.)
#  sigz: variance of fraily variable
#  sige: variance of white noise
#  pp: type of poisson process ('hpp' or 'nhpp')
#' @export
simple.datf = function(param){
  grid = param[['grid']]
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  pp = param[['pp']]
  mu = param[['mu']]
  f1 = param[['f1']]
  f2 = param[['f2']]
  f3 = param[['f3']]
  tau = param[['tau']]
  omg = param[['omg']]
  mix = param[['mix']]

  # parameters
  if (sigz == 1){
    shape=1;scale=1;
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  } else if (sigz == 0.5) {
    shape=2;scale = .5;
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  # generate frailty Z
  if(dist=="gamma"){
    Zstar = rgamma(n, shape=shape, scale=scale)
  }
  if(dist=="lognormal"){
    Zstar = rlnorm(n, lognm, logsd)
  }
  Z = (1-rho) + rho*Zstar

  # generate censoring time
  Cstar = Zstar
  if(C==1){ # pos. corr.
    Cstar[Cstar>tau] = tau
    Cstar = Cstar/tau
  }else if(C==2){ # neg. corr.
    Cstar[Cstar<tau] = tau
    Cstar = tau/Cstar
  }
  u = runif(n, 0.5 ,1)
  if(mix) {
    mix = rbinom(n, 1, omg)
    C = mix*Cstar + (1-mix)*u
  } else {
    C = omg*Cstar + (1-omg)*u
  }

  # generate baseline covariates
  a1 = rbinom(n,1,0.5)
  a2 = runif(n)
  A = cbind(a1, a2)

  # Counting process N(t)
  # Marginal proportion rate
  B = c(-0.5, 0.5)
  if (pp=="hpp"){
    m = rpois(n,10*C*Z* c(exp(A%*%B)))

    # idx = which(m>=1); C = C[idx]; Z = Z[idx]; m = m[idx]; n = length(m)
    # idx = which(m>=2); C = C[idx]; Z = Z[idx]; m = m[idx]; n = length(m)
    # m[m<2] = 2
    # m[m<1] = 1 
    m = m+2
    # m[m<1] = 1
    t = lapply(seq_along(m), function(i) sort(runif(m[i], 0, C[i])))
    # t = lapply(seq_along(m), function(i) c(0,sort(runif(m[i], 0, C[i])),C[i])); m=m+2
  }else if (pp == "nhpp"){ # scaled wang's simulation
    t <- lapply(1:n, function(i){
        lambda <- function(t) Z[i]*(6 + (10*t-6)^3/36) * c(exp(X[i,]%*%B))
        res_1 <- get_nhpp_realization(lambda, C[i])
    })
    m <- sapply(1:n, function(i) length(t[[i]]))
    idx = which(m>=2); C = C[idx]; Z = Z[idx]; t = t[idx]; m = m[idx]; n = length(m)
  }

  # Generate time-dependent functional covariate
  a3 = lapply(1:n, function(i) t[[i]] + rmvnorm(n=1, sigma=C.matern(t[[i]])))

  muf = function(i) Z[i]*mu(t[[i]]) + a1[i]*f1(t[[i]]) + a2[i]*f2(t[[i]]) + a3[[i]]*f3(t[[i]])

  sigx = lapply(1:n,function(i) sigxf(t[[i]])) # lapply(1:n,function(i) one.cor(t[[i]])) #
  mu = lapply(1:n,function(i) muf(i))
  X = lapply(1:n, function(i) rmvnorm(n=1, sigma=sigx[[i]]))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sqrt(sige)))
  Y = lapply(1:n, function(i){
    mu[[i]] + eps[[i]] + X[[i]]
  })

  ts = unlist(t)
  A3 = unlist(a3)
  ord = order(ts)
  mus = unlist(mu)
  ys = unlist(Y)
  Xs = unlist(X)
  Zs = Cs = id = A1 = A2 = c()
  for(i in seq_along(m)){id = c(id, rep(i, m[i]))}
  for(i in seq_along(m)){Cs = c(Cs, rep(C[i], m[i]))}
  for(i in seq_along(m)){Zs = c(Zs, rep(Z[i], m[i]))}
  for(i in seq_along(m)){A1 = c(A1, rep(a1[i], m[i]))}
  for(i in seq_along(m)){A2 = c(A2, rep(a2[i], m[i]))}
  data=data.frame(argvals=ts, subj=id, y=ys, A1=A1, A2=A2, A3=A3, C=Cs, Z=Zs, X=Xs, mu=mus)
  return(list(data=data, c.time=C, t=t, Z=Z))
}

#' @export
simple.true = function(param){
  grid = param[['grid']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  mu = param[['mu']]
  f1 = param[['f1']]
  f2 = param[['f2']]
  f3 = param[['f3']]
  tau = param[['tau']]
  omg = param[['omg']]
  ca = param[['ca']]
  cb = param[['cb']]

  # parameters
  # sigz=0.50: a=1.747273,b=1.229385
  # sigz=0.75: a=0.9051419,b=0.3007614
  if (sigz == 1){
    param$shape=1;param$scale=1;
    param$logmn = log(1^2 / sqrt(1 + 1^2))
    param$logsd = sqrt(log(1 + (1/ 1^2)))
  } else if (sigz == 0.5) {
    param$shape=2;param$scale = .5;
    param$logmn <- log(1^2 / sqrt(.5 + 1^2))
    param$logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  if(dist=="gamma" || dist=='tgamma'){
    param$beta = 1 / param$scale
    param$alpha = param$shape
    param$kappa = 1
  } else if (dist=='lognormal'){
    param$kappa = 1
  }

  if (C==1 || C==2){
    z.conditional = function(t){
      s_z_cgt(t, param)
    }
    z2.conditional = function(t){
      s_z2_cgt(t, param)
    }

    Cinv.conditional = function(t){
      s_cinv_cgt(t, param)
    }

    zCinv.conditional = function(t){
      s_zcinv_cgt(t, param)
    }
  }

  true.scale = function(t) sapply(t, function(s) z.conditional(s))
  subj.scale = function(t) sapply(t, function(s) zCinv.conditional(s)/Cinv.conditional(s))
  E_cov = exp(0.5)-exp(-0.5)
  obs.scale = function(t) sapply(t, function(s) (10*E_cov*z2.conditional(s)+2*zCinv.conditional(s))/(10*E_cov*z.conditional(s)+2*Cinv.conditional(s)))

  # plot(grid, zCinv.conditional(grid))
  # plot(grid, Cinv.conditional(grid))

  muf = mu(grid); f1=f1(grid); f2=f2(grid); f3=f3(grid)
  trus = true.scale(grid)
  subs = subj.scale(grid)
  if (any(is.na(subs))){
    subs = grid*0
  }
  obss = obs.scale(grid)
  truemu = muf * trus
  subjmu = muf * subs
  obsmu  = muf * obss

  # Covariance mu + X
  sigmuf = function(ts){
    mat = matrix(0,length(ts),length(ts))
    temp = sapply(ts, function(t) z2.conditional(t) - z.conditional(t)^2)
    for(i in seq_along(ts)){
      for(j in i:length(ts)){
        s = ts[i]
        t = ts[j]
        mat[i,j] = temp[j] * mu(s)*mu(t)
      }
    }
    mat = t(mat) + mat
    diag(mat) = diag(mat)/2
    return(mat)
  }


  truesigx  = sigxf(grid) #one.cor(grid) #
  truesigmu = sigmuf(grid)
  truecov   = truesigmu + truesigx

  return(list(truemu=truemu, truecov=truecov, subjmu=subjmu, obsmu=obsmu,
              f1=f1, f2=f2, f3=f3))
}

#' @export
simple.tru.mean = function(param){
  grid = param[['grid']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  mu = param[['mu']]
  f1 = param[['f1']]
  f2 = param[['f2']]
  f3 = param[['f3']]
  tau = param[['tau']]
  omg = param[['omg']]
  ca = param[['ca']]
  cb = param[['cb']]

  # parameters
  # sigz=0.50: a=1.747273,b=1.229385
  # sigz=0.75: a=0.9051419,b=0.3007614
  if (sigz == 1){
    param$shape=1;param$scale=1;
    param$logmn = log(1^2 / sqrt(1 + 1^2))
    param$logsd = sqrt(log(1 + (1/ 1^2)))
  } else if (sigz == 0.5) {
    param$shape=2;param$scale = .5;
    param$logmn <- log(1^2 / sqrt(.5 + 1^2))
    param$logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  if(dist=="gamma" || dist=='tgamma'){
    param$beta = 1 / param$scale
    param$alpha = param$shape
    param$kappa = 1
  } else if (dist=='lognormal'){
    param$kappa = 1
  }

  if (C==1 || C==2){
    z.conditional = function(t){
      s_z_cgt(t, param)
    }
    z2.conditional = function(t){
      s_z2_cgt(t, param)
    }

    Cinv.conditional = function(t){
      s_cinv_cgt(t, param)
    }

    zCinv.conditional = function(t){
      s_zcinv_cgt(t, param)
    }
  }

  true.scale = function(t) sapply(t, function(s) z.conditional(s))
  subj.scale = function(t) sapply(t, function(s) zCinv.conditional(s)/Cinv.conditional(s))
  E_cov = exp(0.5)-exp(-0.5)
  obs.scale = function(t) sapply(t, function(s) (10*E_cov*z2.conditional(s)+2*zCinv.conditional(s))/(10*E_cov*z.conditional(s)+2*Cinv.conditional(s)))

  # plot(grid, zCinv.conditional(grid))
  # plot(grid, Cinv.conditional(grid))

  muf = mu(grid); f1=f1(grid); f2=f2(grid); f3=f3(grid)
  trus = true.scale(grid)
  # subs = subj.scale(grid)
  # if (any(is.na(subs))){
  #   subs = grid*0
  # }
  # obss = obs.scale(grid)
  coef.funs = f1(grid) + f2(grid) + f3(grid)
  truemu = muf * trus + coef.funs
  # subjmu = muf * subs + coef.funs
  # obsmu  = muf * obss + coef.funs

  return(truemu)
}

#' @export
dyn.dat = function(dat){
  data = dat$data

  unique_id = unique(data$subj)
  mus = Xs = ys = ts = Zs = Cs = A1 = A2 = A3 = ids = c()
  id = 1600
  for (i in seq_along(unique_id)) {
    sel = which(data$subj == i)
    while (length(sel) >= 2) {
      id = id + 1
      ys = c(ys, data$y[sel])
      mus = c(mus, data$mu[sel])
      Xs = c(Xs, data$X[sel])
      ts = c(ts, data$argvals[sel])
      Zs = c(Zs, data$Z[sel])
      Cs = c(Cs, data$C[sel])
      A1 = c(A1, data$A1[sel])
      A2 = c(A2, data$A2[sel])
      A3 = c(A3, data$A3[sel])
      ids = c(ids, rep(id, length(sel)))

      sel = sel[-length(sel)]
    }
  }

  ys.na = ys
  for (i in sort(unique(ids))) {
    sel = max(which(ids == i))
    ys.na[sel] = NA
  }
  
  tru.dat = data.frame(argvals=ts, subj=ids, y=ys, A1=A1, A2=A2, A3=A3, C=Cs, Z=Zs)
  dyn.dat = data.frame(argvals=ts, subj=ids, y=ys.na, A1=A1, A2=A2, A3=A3, C=Cs, Z=Zs)
  sel = is.na(ys.na)
  return(list(tru.dat=tru.dat, dyn.dat=dyn.dat, mu=(mus+Xs)[sel], sel=sel))
}

# Simulate informative censoring data with covariate
# param: a list of setting values:
#  grid: grid of time point in [0,1]
#  n: number of subjects
#  rho: frailty parameter (0-1)
#  dist: distribution of frailty variable ('lognormal' or 'gamma')
#  C: censoring time (2: negatively correlatve with Zstar, 1: positively correlatve with Zstar.)
#  sigz: variance of fraily variable
#  sige: variance of white noise
#  pp: type of poisson process ('hpp' or 'nhpp')
#' @export
bayes.datf = function(param){
  grid = param[['grid']]
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  pp = param[['pp']]
  mu = param[['mu']]
  f1 = param[['f1']]
  f2 = param[['f2']]
  f3 = param[['f3']]
  tau = param[['tau']]
  omg = param[['omg']]
  ca = param[['ca']]
  cb = param[['cb']]


    # parameters
  if (sigz == 1){
    shape=1;scale=1;
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  } else if (sigz == 0.5) {
    shape=2;scale = .5;
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  # generate frailty Z
  if(dist=="gamma"){
    Zstar = rgamma(n, shape=shape, scale=scale)
  }
  if(dist=="lognormal"){
    Zstar = rlnorm(n, lognm, logsd)
  }
  Z = (1-rho) + rho*Zstar

  # generate censoring time
  if(C==1){
    Cstar = sapply(Zstar, function(z)rbeta(1, ca+tau*z, cb))
  }else if(C==2){
    Cstar = sapply(Zstar, function(z)rbeta(1, ca, cb+tau*z))
  } else if(C==3){ # neg. corr.
    Cstar = sapply(Zstar, function(z)rexp(1, ca*z))
    Cstar[Cstar>tau] = tau
    Cstar = Cstar/tau
  } else if(C==4){ # positive corr.
    Cstar = sapply(Zstar, function(z)rexp(1, ca*z))
    Cstar[Cstar<tau] = tau
    Cstar = tau/Cstar
  }
  C = Cstar


  # generate baseline covariates
  a1 = runif(n)
  a2 = rbinom(n,1,0.5)
  A = cbind(a1, a2)

  # Counting process N(t)
  # Marginal proportion rate
  B = c(0.5, -0.5)
  if (pp=="hpp"){
    m = rpois(n,10*C*Z* c(exp(A%*%B)))
    m = m+2
    t = lapply(seq_along(m), function(i) sort(runif(m[i], 0, C[i])))
  }else if (pp == "nhpp"){ # scaled wang's simulation
    t <- lapply(1:n, function(i){
        lambda <- function(t) Z[i]*(6 + (10*t-6)^3/36) * c(exp(X[i,]%*%B))
        res_1 <- get_nhpp_realization(lambda, C[i])
    })
    m <- sapply(1:n, function(i) length(t[[i]]))
    idx = which(m>=2); C = C[idx]; Z = Z[idx]; t = t[idx]; m = m[idx]; n = length(m)
  }

  # Generate time-dependent covariate
  a3 = lapply(m, function(n) rnorm(n, 0.5, 0.5))
  muf = function(i) Z[i]*mu(t[[i]]) + a1[i]*f1(t[[i]]) + a2[i]*f2(t[[i]]) + a3[[i]]*f3(t[[i]])

  sigx = lapply(1:n,function(i) sigxf(t[[i]])) # lapply(1:n,function(i) one.cor(t[[i]])) #
  mu = lapply(1:n,function(i) muf(i))
  X = lapply(1:n, function(i) rmvnorm(n=1, sigma=sigx[[i]]))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sqrt(sige)))
  Y = lapply(1:n, function(i){
    mu[[i]] + eps[[i]] + X[[i]]
  })

  ts = unlist(t)
  A3 = unlist(a3)
  ord = order(ts)
  mus = unlist(mu)
  ys = unlist(Y)
  Zs = Cs = id = A1 = A2 = c()
  for(i in seq_along(m)){id = c(id, rep(i, m[i]))}
  for(i in seq_along(m)){Cs = c(Cs, rep(C[i], m[i]))}
  for(i in seq_along(m)){Zs = c(Zs, rep(Z[i], m[i]))}
  for(i in seq_along(m)){A1 = c(A1, rep(a1[i], m[i]))}
  for(i in seq_along(m)){A2 = c(A2, rep(a2[i], m[i]))}
  data=data.frame(argvals=ts, subj=id, y=ys, A1=A1, A2=A2, A3=A3, C=Cs, Z=Zs)
  return(list(data=data, c.time=C, t=t, Z=Z))
}
#' @export
bayes.true = function(param){
  grid = param[['grid']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  mu = param[['mu']]
  f1 = param[['f1']]
  f2 = param[['f2']]
  f3 = param[['f3']]
  tau = param[['tau']]
  omg = param[['omg']]
  ca = param[['ca']]
  cb = param[['cb']]

  # parameters
  # sigz=0.50: a=1.747273,b=1.229385
  # sigz=0.75: a=0.9051419,b=0.3007614
  if (sigz == 1){
    param$shape=1;param$scale=1;
    param$logmn = log(1^2 / sqrt(1 + 1^2))
    param$logsd = sqrt(log(1 + (1/ 1^2)))
  } else if (sigz == 0.5) {
    param$shape=2;param$scale = .5;
    param$logmn <- log(1^2 / sqrt(.5 + 1^2))
    param$logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  if(dist=="gamma" || dist=='tgamma'){
    param$beta = 1 / param$scale
    param$alpha = param$shape
    param$kappa = 1
  } else if (dist=='lognormal'){
    param$kappa = 1
  }

  lambda = 0.1
  if (C==1 || C==2 || C==3 || C==4){
    z.conditional = function(t){
      b_z_cgt(t, param)
    }
    z2.conditional = function(t){
      b_z2_cgt(t, param)
    }

    Cinv.conditional = function(t){
      b_cinv_cgt(t, param)
    }

    zCinv.conditional = function(t){
      b_zcinv_cgt(t, param)
    }
  }
  
  true.scale = function(t) sapply(t, function(s) z.conditional(s))
  subj.scale = function(t) sapply(t, function(s) zCinv.conditional(s)/Cinv.conditional(s))
  E_cov = exp(0.5)-exp(-0.5)
  obs.scale = function(t) sapply(t, function(s) (10*E_cov*z2.conditional(s)+2*zCinv.conditional(s))/(10*E_cov*z.conditional(s)+2*Cinv.conditional(s)))

  # plot(grid, zCinv.conditional(grid))
  # plot(grid, Cinv.conditional(grid))

  muf = mu(grid); f1=f1(grid); f2=f2(grid); f3=f3(grid)
  trus = true.scale(grid)
  subs = subj.scale(grid)
  if (any(is.na(subs))){
    subs = grid*0
  }
  obss = obs.scale(grid)
  truemu = muf * trus
  subjmu = muf * subs
  obsmu  = muf * obss

  # Covariance mu + X
  sigmuf = function(ts){
    mat = matrix(0,length(ts),length(ts))
    temp = sapply(ts, function(t) z2.conditional(t) - z.conditional(t)^2)
    for(i in seq_along(ts)){
      for(j in i:length(ts)){
        s = ts[i]
        t = ts[j]
        mat[i,j] = temp[j] * mu(s)*mu(t)
      }
    }
    mat = t(mat) + mat
    diag(mat) = diag(mat)/2
    return(mat)
  }

  
  truesigx  = sigxf(grid) #one.cor(grid) #
  truesigmu = sigmuf(grid)
  truecov   = truesigmu + truesigx

  return(list(truemu=truemu, truecov=truecov, subjmu=subjmu, obsmu=obsmu,
              f1=f1, f2=f2, f3=f3))
}

#' @export
baseline = function(fitted) predict(fitted,data.frame(argvals=x, A1=0, A2=0, A3=0))
#' @export
get_coef = function(fitted) coef(fitted)[2:3]
#' @export
rise_diff = function(v1, v2) mean((v1-v2)^2)/mean(v1^2)


# Correlation between C and Z
cor.cz = function(param){
  grid = seq(0.025, 0.975, 0.010)
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  pp = param[['pp']]
  if(sigz==1){
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  }
  if(sigz==.5){
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }
  if (sigz > 0){
    scale = sigz; shape = 1/scale
  }

  if(sigz==0){
    Zstar = rep(1,n)
  }  else {
    if(dist=="gamma"){
      beta = 1 / scale
      alpha = shape
      Zstar = rgamma(n,shape, scale=scale)
      incomplete.gamma = function(t, alpha, beta) gamma(alpha)*pgamma(beta*t, alpha, 1, lower = T)
      z.conditional = function(t){
        temp = zstar(t)
        if(temp==Inf) return(1)
        if(temp==0) return(0)
        lower.gamma = incomplete.gamma(temp, alpha, beta)
        temp = temp*beta
        alpha/beta - temp^alpha*exp(-temp) / (lower.gamma*beta)
      }
    }
    if(dist=="lognormal") Zstar = rlnorm(n,logmn,logsd)
    if(dist=="uniform"){
      Zstar = runif(n,0,2)
      z.conditional = function(t){
        temp = zstar(t)
        temp/2
      }
    }
  }
  Z = (1-rho) + rho*Zstar


  if(C==0){
    C = rep(1, n)
  }else if(C==1){
    C = runif(n)
  }else if(C==2){
    u1 = runif(n)
    u2 = runif(n,.2,1)
    C = ifelse(u1 < .2, 1, u2)
  }else if(C==3){
    C = sapply(seq(n), function(i) rtrunc(1, 'exp', rate = 1 / Z[i], a = 0.1, b = 1))
  }else if(C==4){
    C = 1 / log(Z + exp(1) - (1-rho))
    zstar = function(t) (exp(1/t) - exp(1))/rho
  }else if(C==5){
    C = 1 / (Z + 1 - (1-rho))
    zstar = function(t) (1/t - 1 )/rho
    subj.scale = function(t){
      temp = zstar(t)
      if(temp==Inf) return(1)
      if(temp==0) return(0)
      temp1 = incomplete.gamma(temp, alpha, beta)
      temp2 = incomplete.gamma(temp, alpha+1, beta)
      temp3 = incomplete.gamma(temp, alpha+2, beta)

      up = (1-rho)*beta^2*temp1 + (2-rho)*rho*beta*temp2 + rho^2*temp3
      down = beta^2*temp1 + rho*beta*temp2
      up/down
    }
    obs.scale = function(t){
      temp = zstar(t)
      if(temp==Inf) return(1)
      if(temp==0) return(0)
      temp1 = incomplete.gamma(temp, alpha, beta)
      temp2 = incomplete.gamma(temp, alpha+1, beta)
      temp3 = incomplete.gamma(temp, alpha+2, beta)

      up = (1-rho)^2*beta^2*temp1 + 2*(1-rho)*rho*beta*temp2 + rho^2*temp3
      down = (1-rho)*beta^2*temp1 + rho*beta*temp2
      up/down
    }
  }else if(C==6){
    a = 10
    C = 1 / (a*Z + 1 - a*(1-rho))
    zstar = function(t) (1/t - 1 )/(a*rho)
  }else if(C==10){
    lambda = 0.5
    C = 1 / sapply(seq(n), function(i) 1+rexp(1,lambda*Zstar[i]))
    z.conditional = function(t){
      temp1 = 1/t - 1
      temp2 = beta / (beta + lambda*temp1)
      shape*scale* (1-temp2^(shape + 1)) / (1-temp2^(shape))
    }

    c.conditional = function(t){
      temp1 = 1/t - 1
      temp2 = beta / (beta + lambda*temp1)
      if (alpha == 1) {
        temp3 = log(lambda*temp1/beta+1)/lambda
      } else {
        temp3 = - 1/(alpha-1)*((lambda*temp1+beta)^{1-alpha}-beta^{1-alpha})/lambda
      }
      beta^alpha/(1-temp2^alpha) * (-temp1*(lambda*temp1+beta)^{-alpha} + temp3)
    }

    zc.conditional = function(t){
      temp1 = 1/t - 1
      temp2 = beta / (beta + lambda*temp1)
      1/lambda - ((alpha*temp1)/(beta+lambda*temp1) + 1/lambda) * temp2^alpha
    }

    subj.scale = function(t){
      z=z.conditional(t); c=c.conditional(t); zc = zc.conditional(t)
      (rho*zc + (1-rho)*c + rho*z + (1-rho)) / (1+c)
    }
  }
  true.scale = function(t){
    z=z.conditional(t)
    rho*z + (1-rho)
  }

  c(cor(C, Z), var(C), var(Z))
}

# Simulate the Observation and Censoring Only
datf.oc = function(param){
  grid = seq(0.025, 0.975, 0.025)
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  pp = param[['pp']]

  if(sigz==1){
    shape=1;scale=1;
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  }
  if(sigz==.5){
    shape=2;scale = .5;
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  if(sigz==0){
    Zstar = rep(1,n)
  }  else {
    if(dist=="gamma") Zstar = rgamma(n,shape, scale=scale)
    if(dist=="lognormal") Zstar = rlnorm(n,logmn,logsd)
    if(dist=="uniform") Zstar = runif(n,0.5,1.5)
  }
  Z = (1-rho) + rho*Zstar


  if(C==0){
    C = rep(1, n)
  }else if(C==1){
    C = runif(n)
  }else if(C==2){
    u1 = runif(n)
    u2 = runif(n,.2,1)
    C = ifelse(u1 < .2, 1, u2)
  }else if(C==3){
    C = sapply(seq(n), function(i) rtrunc(1, 'exp', rate = 1/Z[i], a = 0.1, b = 1))
  }

  if (pp=="hpp"){
    m = rpois(n,10*Z*C)
    m[m<2] = 2
    t = lapply(seq_along(m), function(i) sort(runif(m[i], 0, C[i])))
  }else if (pp == "nhpp"){ # scaled wang's simulation
    t <- lapply(1:n, function(i){
        lambda <- function(t) Z[i]*(6 + (6-6*t)^3/36)
        res_1 <- get_nhpp_realization(lambda, C[i])
    })
    m <- sapply(1:n, function(i) length(t[[i]]))
    #C = C[m>=2]; Z = Z[m>=2]; t = t[m>=2]; m = m[m>=2]; n = length(m)
  }

  return(list(t = t, m = m, c.time=C))
}

# Simulate for bias for subj
datf.bias.subj = function(param){
  grid = seq(0.025, 0.975, 0.025)
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  pp = param[['pp']]
  if(sigz==1){
    shape=1;scale=1;
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  }
  if(sigz==.5){
    shape=2;scale = .5;
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  if(sigz==0){
    Zstar = rep(1,n)
  }  else {
    if(dist=="gamma") Zstar = rgamma(n,shape, scale=scale)
    if(dist=="lognormal") Zstar = rlnorm(n,logmn,logsd)
  }
  Z = (1-rho) + rho*Zstar

  if (pp=="hpp"){
    CRF = function(t) 10*t
  }else if (pp == "nhpp"){ # scaled wang's simulation
    CRF = function(t) (6*t + (10*t-6)^4/1440 - 1296/1440)
  }


  if(C==0){
    C = rep(1, n)
  }else if(C==1){
    C = runif(n)
  }else if(C==2){
    u1 = runif(n)
    u2 = runif(n,.2,1)
    C = ifelse(u1 < .2, 1, u2)
  }else if(C==3){
    C = sapply(seq(n), function(i){
      temp = rtrunc(100, 'exp', rate = 0.1*Z[i], a = 0.1, b = 1)
      mean(1/CRF(temp))
    } )
  }else if(C==4){
    C = 1 / log(exp(1) + Z)
  }else if(C==5){
    C = sapply(Z, function(u) min(u,1))
  }else if(C==6){
    C = sapply(Z, function(u) min(runif(1,0,u), 1))
  }else if(C==7){
    C = 1/sapply(Z, function(u) max(runif(1,u,10), 1))
  }else if(C==8){
    C = sapply(Z, function(u) min(1/u,1))
  }else if(C==9){
    C = 1 / sapply(seq(n), function(i) max(1,rexp(1,1/Z[i])))
  }else if(C==10){
    C = 1 / sapply(seq(n), function(i) 1+rexp(1,1/Z[i]))
  }else if(C==11){
    C = 1 / sapply(seq(n), function(i) rtrunc(1,'exp',rate=1/Z[i],a=1.2, b = Inf))
  }

  
  mean(C*Z) / mean(C)
}

datf.bias.obs = function(param){
  grid = seq(0.025, 0.975, 0.025)
  n = param[['n']]
  rho = param[['rho']]
  dist = param[['dist']]
  C = param[['C']]
  sigz = param[['sigz']]
  sige = param[['sige']]
  pp = param[['pp']]
  if(sigz==1){
    shape=1;scale=1;
    logmn = log(1^2 / sqrt(1 + 1^2))
    logsd = sqrt(log(1 + (1/ 1^2)))
  }
  if(sigz==.5){
    shape=2;scale = .5;
    logmn <- log(1^2 / sqrt(.5 + 1^2))
    logsd <- sqrt(log(1 + (.5 / 1^2)))
  }

  if(sigz==0){
    Zstar = rep(1,n)
  }  else {
    if(dist=="gamma") Zstar = rgamma(n,shape, scale=scale)
    if(dist=="lognormal") Zstar = rlnorm(n,logmn,logsd)
  }
  Z = (1-rho) + rho*Zstar

  if (pp=="hpp"){
    CRF = function(t) 10*t
  }else if (pp == "nhpp"){ # scaled wang's simulation
    CRF = function(t) (6*t + (10*t-6)^4/1440 - 1296/1440)
  }


  if(C==0){
    C = rep(1, n)
  }else if(C==1){
    C = runif(n)
  }else if(C==2){
    u1 = runif(n)
    u2 = runif(n,.2,1)
    C = ifelse(u1 < .2, 1, u2)
  }else if(C==3){
    C = sapply(seq(n), function(i){
      temp = rtrunc(100, 'exp', rate = 0.1*Z[i], a = 0.1, b = 1)
      mean(1/CRF(temp))*mean(CRF(temp))
    } )
  }else if(C==4){
    C = 1 / log(exp(1) + Z)
  }else if(C==5){
    C = sapply(Z, function(u) min(u,1))
  }else if(C==6){
    C = sapply(Z, function(u) min(runif(1,0,u), 1))
  }else if(C==7){
    C = 1/sapply(Z, function(u) max(runif(1,u,10), 1))
  }else if(C==8){
    C = sapply(Z, function(u) min(1/u,1))
  }else if(C==9){
    C = 1 / sapply(seq(n), function(i) max(1,rexp(1,1/Z[i])))
  }else if(C==10){
    C = 1 / sapply(seq(n), function(i) 1+rexp(1,1/Z[i]))
  }else if(C==11){
    C = 1 / sapply(seq(n), function(i) rtrunc(1,'exp',rate=1/Z[i],a=1.2, b = Inf))
  }

  
  mean(C*Z^2) / mean(C*Z)
}

get_nhpp_realization <- function(lambda, t_max=10){
    #t_max <- 10
    t <- 0
    Lambda <- function(tupper) integrate(f = lambda, lower = 0, upper = tupper)$value
    Ft <- function(x) Lambda(x)/ Lambda(t_max)
    Ft_inv <- function(u){
        a <- 0
        b <- t_max+0.1
        eps <- 1e-6
        while(abs(a-b)>eps){
            if(Ft((a + b)/2)<=u) a <- (a+b)/2
            if(Ft((a + b)/2)>u) b <- (a+b)/2
        }
        return(0.5*(a+b))
    }
    n <- rpois(1, Lambda(t_max))
    if (n<2) n = 2 # make sure at least 2 observations per subject
    X <- sapply(1:n, function(z) Ft_inv(runif(1)))
    return(sort(X))
}

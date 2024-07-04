# Modified from fcr package

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


createPhi <- function(obj, data, argvals, nPhi = NULL){

    if(is.null(nPhi)){
        nPhi <- ncol(obj$eigenfunctions)
    }

    ut <- obj$argvals.new
    for(i in 1:nPhi){
        cur    <- vapply(data[[argvals]], function(x) obj$eigenfunctions[ut==x,i], numeric(1))
        data[[paste("phi",i,sep="")]] <- cur

        cur_sp <- obj$sigma2/obj$eigenvalues[i]
        assign(paste("sp",i,sep=""),cur_sp , envir = parent.frame())
    }

    data
}


createFormula <- function(obj, formula = NULL, sp = FALSE, nPhi = NULL){
    stopifnot(!is.null(formula))
    if(is.null(nPhi)){
        nPhi <- ncol(obj$eigenfunctions)
    }

    if(!sp){
        phiForm <- paste("s(g, by = phi",1:nPhi,", bs='re') ", collapse="+",sep="")
    }
    if(sp){
        phiForm <- paste("s(g, by = phi",1:nPhi,", bs='re', sp = sp",1:nPhi,") ", collapse="+",sep="")
    }
    model <- as.character(formula)
    paste(model[2], model[1], model[3], "+", phiForm)
}

plot.infcr <- function(x, plot.covariance = FALSE, ...){
        stopifnot(class(x) == "infcr")

        if(plot.covariance) {
                oldpar <- par()$ask
                par(ask=TRUE)

                tnew <- x$face.object$argvals.new
                inx <- which(tnew %in% seq(min(tnew), max(tnew), len = 100))

                image.plot(tnew[inx],tnew[inx], x$face.object$Cor.new[inx,inx],
                           xlab="Functional Domain", ylab="Functional Domain", main = "Correlation Function")

                plot(tnew[inx], x$face.object$Chat.raw.diag.new[inx], xlab="Functional Domain",
                     ylab="", main = "Variance Function", type="l")

                matplot(x$face.object$eigenfunctions[inx,], xlab="Functional Domain", ylab="",
                     main="Eigenfunctions of the Covariance Function", type='l')
                evals <- x$face.object$eigenvalues
                evals <- sprintf("%5.3f",evals)
                evals[evals == "0.000"] <- "<0.000"
                legend("top", legend = paste("eval", 1:length(evals), " = ", evals ,sep=""),
                       col=1:length(evals), lty = 1:length(evals),
                       bty='n',ncol=3)

                par(ask=oldpar)
        }

        if(!plot.covariance) {
                plot(x$fit, ...)
        }
}


predict.face.sparse <- function(object,newdata,...){
  
  ## check inputs
  if(!inherits(object, "face.sparse")) stop("'fit' has to be a face.sparse object")
  check.data(newdata,type="predict") 

  subj.pred = newdata$subj
  subj_unique.pred = unique(subj.pred)
  y.pred = newdata$y
  
  Theta = object$Theta
  knots = object$knots
  npc = object$npc
  p = object$p
  
  center = object$center
  #const.var.error = object$const.var.error
  fit_mean = object$fit_mean
  #fit_var_error = object$fit_var_error
  sigma2 = object$sigma2
  
  mu.pred <- rep(0,nrow(newdata))
  var.error.pred <- rep(max(sigma2,0.000001),nrow(newdata))
  cov.pred <-  matrix(0,length(mu.pred),length(mu.pred))
    
  if(center) {mu.pred <- predict.pspline.face(fit_mean,newdata$argvals) }
  #if(!const.var.error){
  #  var.error.pred <- predict.pspline(fit_var_error,newdata$argvals)
  #  var.error.pred <-  sapply(var.error.pred,function(x) max(x,0))### added Sept 8, 2015 by Luo
  #}

  
  Bnew = spline.des(knots=knots, x=object$argvals.new, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  scores = list(subj=subj_unique.pred,
                scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta))
                )
  
  for(i in 1:length(subj_unique.pred)){
    sel.pred = which(subj.pred==subj_unique.pred[i])
    pred.points <- newdata$argvals[sel.pred]
    mu.predi <- mu.pred[sel.pred]
    var.error.predi <- var.error.pred[sel.pred]
    
    y.predi = y.pred[sel.pred] - mu.predi
    sel.pred.obs = which(!is.na(y.predi))
    obs.points <- pred.points[sel.pred.obs]
     
    if(!is.null(obs.points)){
      var <- mean(var.error.predi[sel.pred.obs])
      if(var==0&length(sel.pred.obs) < npc)
        stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
               cannot be estimated.")
      B3i.pred = spline.des(knots=knots, x=pred.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
      B3i = spline.des(knots=knots, x=obs.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
      Chati = tcrossprod(B3i%*%Theta,B3i)
      
      if(length(sel.pred.obs)==1) Ri = var.error.predi[sel.pred.obs]
      if(length(sel.pred.obs)>1) Ri = diag(var.error.predi[sel.pred.obs])
      Vi.inv = as.matrix(solve(Chati + Ri))
      Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta,B3i.pred))
      Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
      ui =tcrossprod(Theta,B3i)%*%Vi.inv %*%y.predi[sel.pred.obs]
      scores$u[i,] = as.vector(ui)
      y.pred[sel.pred] = as.numeric(Hi%*%y.predi[sel.pred.obs]) + mu.predi
      temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
      if(length(sel.pred.obs) >1){
        cov.pred[sel.pred,sel.pred] = Vi.pred - temp%*%Vi.inv%*%t(temp)
      }
      if(length(sel.pred.obs) ==1){
        cov.pred[sel.pred,sel.pred] = Vi.pred - Vi.inv[1,1]*temp%*%t(temp)
      }
      ## predict scores
      if(object$calculate.scores==TRUE){ 
        temp = matrix(t(object$eigenfunctions),nrow=npc)%*%(as.matrix(Bnew)%*%ui)/sum(object$eigenfunctions[,1]^2)
        temp = as.matrix(temp)
        scores$scores[i,1:npc] = temp[,1]
      }
    }
  }
  
  B = spline.des(knots=knots, x=newdata$argvals, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  Chat.pred = as.matrix(tcrossprod(B%*%Matrix(Theta),B))
  
  
  return(list(object=object,newdata=newdata,
              y.pred = y.pred, mu.pred=mu.pred,
              var.error.pred=var.error.pred,
              rand_eff = scores, cov.pred = cov.pred, se.pred = sqrt(diag(cov.pred)),
              Chat.pred = Chat.pred,Chat.diag.pred = diag(Chat.pred)))  
}

predict.infcr <- function(object, newdata, ...){
  ## error checking
  stopifnot(class(object) == "infcr")
  check.data(newdata,type="predict")


  # predict mean function
  mu.pred = predict(object$fit, newdata, ...)

  # predict subject deviations
  rs = newdata$y - mu.pred
  resdat = data.frame(y=rs, subj=newdata$subj, argvals=newdata$argvals)
  face.pred = predict.face.sparse(object$face.object, resdat)
  
  sel.pred = which(is.na(newdata$y))
  y.pred = as.vector(face.pred$y.pred + mu.pred)
  y.pred.na = y.pred[sel.pred]
  
  return(list(y.pred = y.pred, y.pred.na=y.pred[sel.pred], sel.pred=sel.pred,
              face.pred = face.pred, mu.pred = mu.pred))
}

# Drop the diagonal element, don't estimate error variance; use 0.0000011 as default error variance estimate
# Add aSUBJ weight for informative censor
informative.weights <- function(data,
                                mean.weight = c("aSUBJ", "SUBJ","OBS"), cov.weight = c("aSUBJ", "SUBJ","OBS"),
                                err.weight = c("aSUBJ", "SUBJ","OBS"),
                                crf_smoothed = F){
  
  #########################
  ####step 0: read in data
  #########################
  check.data(data)
  
  t <- data$argvals
  subj <- data$subj
  
  if(is.null(data$C)){ 
    c.time <- rep(1, length(unique(data$subj)))
  } else {
    c.time <- sapply(unique(subj), function(i) data$C[which(subj==i)[1]])}
  t.list <- lapply(unique(subj), function(i) t[which(subj==i)])
  nobs <- table(subj)
  nobs = nobs[which(nobs > 0)]
  
  #########################
  ####step 0.5: initialize weights for mean function and cov function
  #########################
  # calculate mean weights
  if(mean.weight == "aSUBJ"){
    crf = fit.crf(list(c.time=c.time, t=t.list))
    if (crf_smoothed){
      crf = mgcv::gam(crf ~ s(c.time, bs="ps", m=c(2,2), k=10), sp=NULL)$fitted.values;
    }
    # Normalize
    crf <- crf / sum(crf) 
    # Divide by m
    mean.weight <- crf / nobs
  } else if (mean.weight == "SUBJ"){
    mean.weight <- 1 / nobs / length(nobs) 
  } else if (mean.weight == "OBS"){
    mean.weight = rep(1, length(nobs))
  }
  
  # calculate cov weights
  if(cov.weight == "aSUBJ"){
    crf = fit.crf(list(c.time=c.time, m=nobs, t=t.list))
    if (crf_smoothed){
      crf = mgcv::gam(crf ~ s(c.time, bs="ps", m=c(2,2), k=10), sp=NULL)$fitted.values;
    }
    crf2 = crf^2
    # Normalize 
    crf2 <- crf2 / sum(crf2)
    # Divide by m  
    cov.weight <- crf2 / (nobs^2-nobs)
  } else if (cov.weight == "SUBJ"){
    cov.weight <- 1 / (nobs^2-nobs) / length(nobs) 
  } else if (cov.weight == "OBS"){
    cov.weight = rep(1, length(nobs)) / sum(nobs^2-nobs)

  }

  # calculate err variance weights
  if(err.weight == "aSUBJ"){
    err.weight = mean.weight
  } else {
    crf = fit.crf(list(c.time=c.time, m=nobs, t=t.list))
    if (crf_smoothed){
      crf = mgcv::gam(crf ~ s(c.time, bs="ps", m=c(2,2), k=10), sp=NULL)$fitted.values
    }

    # in case crf is too small
    crf[which(crf <= 0)] = 1e-2

    if (err.weight == "SUBJ"){
      err.weight <- 1 / (crf * nobs)
    } else if (err.weight == "OBS"){
      err.weight = (nobs-1) / crf
    }

    # normalize
    err.weight = err.weight / sum(err.weight*nobs)
  }

  # Calculate consistent mean.weight for mean function
  W.mean <- rep(1, length(subj))
  J = length(subj)
  subj_unique = unique(subj)
  for(i in seq_along(subj_unique)){
    seq <- which(subj == subj_unique[i])
    W.mean[seq] <- mean.weight[i]
  }

  # Calculate consistent err.weight for error variance
  W.err <- rep(1, length(subj))
  J = length(subj)
  subj_unique = unique(subj)
  for(i in seq_along(subj_unique)){
    seq <- (1:J)[subj == subj_unique[i]]
    W.err[seq] <- err.weight[i]
  }

  res <- list(mean.weight=mean.weight, cov.weight=cov.weight, err.weight=err.weight, W.mean=W.mean, W.err=W.err) 
  
  class(res) <- "informative.weights"
  return(res)
}


# Drop the diagonal element, don't estimate error variance; use 0.0000011 as default error variance estimate
# Add aSUBJ weight for informative censor
#' @export
infcr <- function(formula, argvals='argvals', subj='subj', argvals.new = NULL, data = NULL, niter=1, 
                sp = FALSE, nPhi = NULL, use_bam = FALSE, discrete = FALSE,
                 mean.weight = c("aSUBJ", "SUBJ","OBS"), cov.weight = c("aSUBJ", "SUBJ","OBS"),
                                err.weight = c("aSUBJ", "SUBJ","OBS"), crf_smoothed = F,
                face.args = list(knots = 12, lower = -3, pve=0.95),
                ...){
        stopifnot(class (data) == "data.frame")
        stopifnot(class(subj) == "character")
        stopifnot(is.list(face.args))
        stopifnot(!any(c("data","newdata","center","argvals.new") %in% names(face.args)))
        stopifnot(all(subj %in% names(data), argvals %in% names(data)))
        if(any(grepl("^phi[0-9]+|^sp[0-9]+",colnames(data)))){
                stop("column names `sp[0-9]+` and `phi[0-9]+ are reserved`")
        }
        if("g" %in% colnames(data)){
                stop("column name `g` is reserved`")
        }

        ## check to see if new predictions are in the range of the data used to fit the model
        if(!is.null(argvals.new)){
           if(max(argvals.new) > max(data[[argvals]]) | min(argvals.new) < min(data[[argvals]])) {
                   warning("Range of arvals.new exceeds the range of the funcitonal domain in the data.")
           }
        }

        model <- as.character(formula)
        stopifnot(model[1] == "~" & length(model) == 3)
        outcome <- model[[2]]

        if(!outcome %in% names(data)) stop("Outcome variable must be named and included in data.")

        #########################
        ####step 0.5: informative weighting
        #########################
        
        inf_w = informative.weights(data, mean.weight, cov.weight, err.weight, crf_smoothed)
        mean_w = inf_w$W.mean
        cov_w = inf_w$cov.weight
        err_w = inf_w$W.err
        
        #########################
        ####step 1: Mean estimating
        #########################
        data$mean_w = mean_w
        fit <- mgcv::gam(as.formula(formula), data = data, weight = mean_w, ...)
        curEst_fx <- fit$fitted.values

        if(length(curEst_fx) != nrow(data)) stop("Method not implemented to handle missing data in model fitting. Please remove missing data and refit.")

        

        data$g  <- factor(data[[subj]])
        if(is.null(argvals.new)){
                ut <- sort(unique(c(data[[argvals]],
                                    seq(min(data[[argvals]]), max(data[[argvals]]), len =100))))
        } else {
                #ut <- sort(unique(c(argvals.new, data[[argvals]],
                #                    seq(min(data[[argvals]]), max(data[[argvals]]), len =100))))
                ut = argvals.new
        }

        #########################
        ####step 2: covariance estimating
        #########################
        run_time <- c()
        Cest1 = NULL
        for(n in seq_len(niter)){
                start_time <- proc.time()

                resid1  <- curEst_fx - data[[outcome]]
                datCest <- data.frame("argvals" = as.vector(data[[argvals]]),
                                      "subj" = as.vector(data[[subj]]),
                                      "y" = as.vector(resid1))

                Cest1  <- do.call("infda.face",
                                  c(list(data = datCest, newdata = datCest,
                                         mean.weight=mean.weight, cov.weight=cov_w, err.weight=err_w,
                                         center = FALSE, argvals.new = ut),
                                    face.args))
                
                # if(is.null(nPhi)){
                #         nPhi <- length(Cest1$eigenvalues)
                # } else {
                #         nPhi <- min(length(Cest1$eigenvalues), nPhi)
                # }
                # message(paste("Number of Eigenfunctions used in estimation =", nPhi))
                # data <- createPhi(Cest1, data = data, argvals = argvals, nPhi = nPhi)
                # if(use_bam){
                #         fit <- mgcv::bam(formula = as.formula(createFormula(Cest1, formula = formula, sp = sp, nPhi = nPhi)),
                #                    data=data, weight = mean_w, discrete = discrete, ...)
                # } else if (!use_bam){
                #         fit <- mgcv::gam(formula = as.formula(createFormula(Cest1, formula = formula, sp = sp, nPhi = nPhi)),
                #                    data=data, weight = mean_w, ...)
                # }

                run_time[[n]] <- proc.time() - start_time

                # if(niter > 1){
                #         coef_names <- names(coef(fit))
                #         phi_names  <- unique(regmatches(coef_names, regexpr("phi[0-9]+", coef_names)))
                #         curEst_fx  <- predict(fit, exclude = paste("s(g):", phi_names, sep=""))
                #         nPhi <- NULL
                #         rm_cols <- which(grepl("^phi[0-9]+",colnames(data)))
                #         data[rm_cols] <- NULL
                # }
                # rm(list=ls()[which(ls() %in% paste("sp",0:10000,sep=""))])

        }

        ret <- list("fit" = fit,
              "face.object" = Cest1,
              "runtime" = run_time,
              "argvals" = argvals,
              "sp" = sp)
        class(ret) <- "infcr"
        return(ret)
}

#' @export
infda.face <- function(data, newdata = NULL,
                      mean.weight = c("aSUBJ","SUBJ","OBS"), cov.weight = c("aSUBJ","SUBJ","OBS"), err.weight = c("aSUBJ", "SUBJ","OBS"), crf_smoothed = F,
                      center=TRUE,argvals.new=NULL, use_mgcv = TRUE, est_cov = TRUE,
                      knots=7, knots.option="equally-spaced",
                      p=3,m=2,lambda=NULL,lambda_mean=NULL,
                      search.length=21,
                      lower=-10,upper=10,
                      calculate.scores=FALSE,pve=0.99){
  
  #########################
  ####step 0: read in data
  #########################
  check.data(data)
  if(!is.null(newdata)){ check.data(newdata,type="predict")}
  
  y <- data$y
  t <- data$argvals
  subj <- data$subj
  tnew <- argvals.new
  if(is.null(tnew)) tnew <- seq(min(t),max(t),length=100)
  fit_mean <- NULL
  knots.initial <- knots

  #########################
  ####step 0.5: informative weighting
  #########################
  if (is.character(mean.weight) && is.character(cov.weight) && is.character(err.weight)) {
    weights = informative.weights(data, mean.weight, cov.weight, err.weight, crf_smoothed)
    mean.weight = weights$W.mean
    cov.weight = weights$cov.weight
    err.weight = weights$W.err
  }
  

  #########################
  ####step 1: demean
  #########################
  r <- y
  mu.new <- rep(0,length(tnew))
  if(center){
    # use mgcv::gam to estimate mean function
    if (use_mgcv){
      fit_mean <- mgcv::gam(y ~ s(t, bs="ps", m=c(2,2), k=10), weights = mean.weight, sp=lambda_mean)
      mu.new = predict(fit_mean,data.frame(t = tnew))
    } else {
      fit_mean <- pspline(data, argvals.new=tnew, weight = mean.weight, knots=knots.initial, lambda=lambda_mean,
                          knots.option=knots.option, lower=-20, upper=-3, search.length = 100)
      mu.new <- fit_mean$mu.new
    }
    
    if (!est_cov){
      return(list(mu.new =mu.new))
    }

    r <- y - fit_mean$fitted.values
  }
  
  #########################
  ####step 2:raw estimates
  #########################  
  raw <- raw.construct(data.frame("argvals" = t, "subj" = subj, "y" = as.vector(r)), include.diag=F)
  C <- raw$C
  st <- raw$st
  N <- raw$st
  N2 <- raw$N2
  W <- raw$W # retrieve weights
  n0 <- raw$n0
  W = lapply(seq_along(W), function(i) W[[i]]*cov.weight[i])
  
  # delta <- Matrix((st[,1]==st[,2]) * 1) # sparse
  
  #########################
  ####step 3: smooth
  #########################
  knots <- construct.knots(t,knots,knots.option,p)
  
  List <- pspline.setting(st[,1],knots=knots,p,m,type="simple",knots.option=knots.option)
  B1 <- List$B
  B1 <- Matrix(B1)
  DtD <- List$P
  
  B2 = spline.des(knots=knots, x=st[,2], ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  c = dim(B1)[2]
  c2 = c*(c+1)/2
  B = Matrix(t(KhatriRao(Matrix(t(B2)),Matrix(t(B1)))))
  G = Matrix(duplication.matrix(c))
  
  BtWB = matrix(0,nrow=c^2,ncol=c^2)
  # Wdelta = c()
  WC = c()
  for(i in 1:n0){
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    B3 = Matrix(matrix(B[seq,],nrow=length(seq)))
    W3 = W[[i]] # don't form a large W
    BtWB = BtWB + crossprod(B3, W3%*%B3)
    # Wdelta <- c(Wdelta,as.matrix(W3 %*% delta[seq]))
    WC <- c(WC,as.matrix(W3 %*% C[seq]))
  }
  
  GtBtWBG = crossprod(G,BtWB%*%G)
  
  BG = B%*%G # sparse
  # detWde <- crossprod(delta,Wdelta) # detWde = sum(delta)
  # GtBtWdelta <- crossprod(BG,Wdelta)
  # XtWX <- rbind(cbind(GtBtWBG,GtBtWdelta), cbind(t(GtBtWdelta),detWde))
  XtWX <- GtBtWBG
  
  eSig = eigen(XtWX,symmetric=TRUE)
  V = eSig$vectors
  E = eSig$values
  E = E + 0.000001*max(E)
  Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)
  
  P = crossprod(G,Matrix(suppressMessages(kronecker(diag(c),DtD))))%*%G
  
  Q = P #bdiag(P,0) 
  tUQU = crossprod(Sigi_sqrt,(Q%*%Sigi_sqrt))
  Esig = eigen(tUQU,symmetric=TRUE)
  
  U = Esig$vectors
  s = Esig$values
  A0 <- Sigi_sqrt%*%U
  X <- cbind(BG)
  A = as.matrix(X%*%A0) # F=XA dense
  
  AtA = crossprod(A) # diff
  f = crossprod(A,C) # diff
  ftilde = crossprod(A,WC) # diff
  
  c2 <- c2 #+ 1
  g <- rep(0, c2)
  G1 <- matrix(0,c2,c2)
  mat_list <- list()
  
  
  for(i in 1:n0){
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    Ai = matrix(A[seq,],nrow=length(seq))
    #AitAi = crossprod(Ai) #t(Ai)%*%Ai
    Wi = W[[i]]
    
    fi = crossprod(Ai,C[seq]) # t(Fi)Ci
    Ji = crossprod(Ai,Wi%*%C[seq])
    Li = crossprod(Ai,Wi%*%Ai)
    AitAi = Li
    if(i==1) {AtA = 0*AtA +  AitAi}
    else{AtA = AtA + AitAi}
    
    g = g + Ji*Ji
    G1 = G1 + AitAi*(Ji%*%t(ftilde))
    
    LList <- list()
    LList[[1]] = AitAi
    LList[[2]] = Li
    mat_list[[i]] = LList
    
  }
  
  Lambda <- seq(lower,upper,length=search.length)
  Gcv <- 0*Lambda
  gcv <- function(x){
    lambda <- exp(x)
    d <- 1/(1+lambda*s)
    ftilde_d <- ftilde*d
    cv0 <- -2*sum(ftilde_d*ftilde)
    cv1 <-  sum(ftilde_d*(AtA%*%ftilde_d))
    cv2 <-  2*sum(d*g)
    cv3 <-  -4*sum(d*(G1%*%d))
    cv4 <- sum(unlist(sapply(mat_list,function(x){
      a <- x[[1]]%*%ftilde_d
      b <- x[[2]]%*%ftilde_d
      2*sum(a*b*d)
    })))
    cv <- cv0 + cv1 + cv2 + cv3 + cv4
    return(cv)
  }
  if(is.null(lambda)){
    Lambda <- seq(lower,upper,length=search.length)
    Length <- length(Lambda)
    Gcv <- rep(0,Length)
    for(i in 1:Length) 
      Gcv[i] <- gcv(Lambda[i])
    i0 <- which.min(Gcv)
    lambda <- exp(Lambda[i0])
  }
  
  alpha <- matrix.multiply(A0,1/(1+lambda*s))%*%ftilde
  Theta <- G %*% alpha[1:c2] #c2-1
  Theta <- matrix(Theta,c,c)         # parameter estimated (sym)


  # estimating the homogeneous error variance
  B_old = spline.des(knots=knots, x=data$argvals, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  C.diag.est = rowSums((B_old%*%Matrix(Theta)) * B_old)
  sigma2 <- sum((as.vector(r)^2 - C.diag.est) * err.weight) / sum(err.weight)
  tryCatch(
    expr = {
      if(sigma2 <= 0.000001) {                                               
        warning("error variance cannot be non-positive, reset to 1e-6!")    
        sigma2 <- 0.000001             
      }
    },
    error = function(e) {
      print(e)
      print(mean.weight)
      # print(c(param$omg, param$rho, param$sigz))
      print(sigma2)
      print(data$argvals)
      print(data$subj)
      # print(data$C)
      save(data, file = paste(mean.weight,param$omg, param$rho, param$sigz,'.RData'))
      print(format(round(err.weight, 3), nsmall = 3))
      sigma2 <- 0.000001
    },
    warning = function(w) {
      print(w)
    }
  )
  
  
  Eigen <- eigen(Theta,symmetric=TRUE)
  Eigen$values[Eigen$values<0] <- 0
  npc <- sum(Eigen$values>0) #which.max(cumsum(Eigen$values)/sum(Eigen$values)>pve)[1]
  if(npc >1){
    Theta <- matrix.multiply(Eigen$vectors[,1:npc],Eigen$values[1:npc])%*%t(Eigen$vectors[,1:npc])
    Theta_half <- matrix.multiply(Eigen$vectors[,1:npc],sqrt(Eigen$values[1:npc]))
  }
  if(npc==1){
    Theta <- Eigen$values[1]*suppressMessages(kronecker(Eigen$vectors[,1],t(Eigen$vectors[,1])))
    Theta_half <- sqrt(Eigen$values[1])*Eigen$vectors[,1]
  }
  Eigen <- eigen(Theta,symmetric=TRUE)
  
  #########################
  ####step 4: calculate estimated covariance function
  #########################
  Bnew = spline.des(knots=knots, x=tnew, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  Gmat <- crossprod(Bnew) / nrow(Bnew)
  eig_G <- eigen(Gmat, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
  G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)
  
  Chat.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta),Bnew)) 
  Chat.diag.new = as.vector(diag(Chat.new))  
  Cor.new = diag(1/sqrt(Chat.diag.new))%*%Chat.new%*%diag(1/sqrt(Chat.diag.new))
  
  Eigen.new = eigen(as.matrix(G_half%*%Matrix(Theta)%*%G_half),symmetric=TRUE)
  # Eigen.new = eigen(Chat.new,symmetric=TRUE)
  
  npc = which.max(cumsum(Eigen$values)/sum(Eigen$values)>pve)[1] #determine number of PCs
  # 
  eigenfunctions = matrix(Bnew%*%G_invhalf%*%Eigen.new$vectors[,1:min(npc,length(tnew))],ncol=min(npc,length(tnew)))
  eigenvalues = Eigen.new$values[1:min(npc,length(tnew))]
  # eigenfunctions = eigenfunctions*sqrt(length(tnew))/sqrt(max(tnew)-min(tnew))
  # eigenvalues = eigenvalues/length(tnew)*(max(tnew)-min(tnew))


  #########################
  ####step 5: calculate variance
  #########################
  var.error.hat <- rep(sigma2,length(t))
  var.error.new <- rep(sigma2,length(tnew))
  # 
  # 
  # 
  Chat.raw.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta),Bnew)) + diag(var.error.new) 
  Chat.raw.diag.new = as.vector(diag(Chat.raw.new)) 
  Cor.raw.new = diag(1/sqrt(Chat.raw.diag.new))%*%Chat.raw.new%*%diag(1/sqrt(Chat.raw.diag.new))
  
  #########################
  ####step 6: prediction
  #########################
  if(!is.null(newdata)){
    
    mu.pred <- rep(0,length(newdata$argvals))
    var.error.pred <- rep(sigma2,length(newdata$argvals))
    if(center){
      mu.pred <- predict.pspline.face(fit_mean,newdata$argvals)
    }
    
    subj.pred = newdata$subj
    subj_unique.pred = unique(subj.pred)
    y.pred = newdata$y
    Chat.diag.pred = 0*y.pred
    se.pred = 0*y.pred
    
    scores = list(subj=subj_unique.pred,
                  scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                  u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta))
    )
    
    for(i in 1:length(subj_unique.pred)){
      sel.pred = which(subj.pred==subj_unique.pred[i])
      lengthi = length(sel.pred)
      
      pred.points <- newdata$argvals[sel.pred]
      mu.predi <- mu.pred[sel.pred]
      var.error.predi <- var.error.pred[sel.pred]
      
      y.predi = y.pred[sel.pred] - mu.predi
      sel.pred.obs = which(!is.na(y.predi))
      obs.points <- pred.points[sel.pred.obs]
      if(!is.null(obs.points)){
        var <- mean(var.error.predi[sel.pred.obs])
        if(var==0&length(sel.pred.obs) < npc)
          stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
               cannot be estimated.")
        B3i.pred = spline.des(knots=knots, x=pred.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
        B3i = spline.des(knots=knots, x=obs.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
        Chati = tcrossprod(B3i%*%Theta,B3i)
        Chat.diag.pred[sel.pred] = diag(Chati)
        if(length(sel.pred.obs)==1) Ri = var.error.predi[sel.pred.obs]
        if(length(sel.pred.obs)>1) Ri = diag(var.error.predi[sel.pred.obs])
        Vi.inv = as.matrix(solve(Chati + Ri))
        Vi.pred = tcrossprod(B3i.pred%*%Theta,B3i.pred)
        Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
        ui =tcrossprod(Theta,B3i)%*%Vi.inv %*%y.predi[sel.pred.obs]
        scores$u[i,] = as.vector(ui)
        y.pred[sel.pred] = as.numeric(Hi%*%y.predi[sel.pred.obs]) + mu.predi
        temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
        if(length(sel.pred.obs) >1){
          se.pred[sel.pred] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))
        }
        if(length(sel.pred.obs) ==1){
          se.pred[sel.pred] = sqrt(Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp))
        }
        
        ## predict scores
        if(calculate.scores==TRUE){ 
          temp = matrix(t(eigenfunctions),nrow=npc)%*%(as.matrix(Bnew)%*%ui)/sum(eigenfunctions[,1]^2)
          temp = as.matrix(temp)
          scores$scores[i,1:npc] = temp[,1]
        }
      }
    }
  }## if(is.null(newdata))
  if(is.null(newdata)){
    y.pred=NULL
    mu.pred = NULL
    var.error.pred = NULL
    Chat.diag.pred = NULL
    se.pred = NULL
    scores=NULL
    
  }
  
  res <- list(newdata=newdata, mu.new = mu.new, Chat.new=Chat.new,
              W = W, Theta=Theta, sigma2=sigma2, 
              argvals.new=tnew, lambda=lambda, Gcv=Gcv, 
              y.pred = y.pred,var.error.new = var.error.new,
              Cor.new = Cor.new, eigenfunctions = eigenfunctions, eigenvalues = eigenvalues,
              Cor.raw.new = Cor.raw.new, Chat.raw.diag.new = Chat.raw.diag.new,
              scores = scores, calculate.scores=calculate.scores,
              mu.hat = fit_mean$fitted.values,var.error.hat = var.error.hat,
              mu.pred = mu.pred, var.error.pred = var.error.pred, Chat.diag.pred = Chat.diag.pred,
              se.pred = se.pred,
              fit_mean = fit_mean,
              Lambda=Lambda,knots=knots,knots.option=knots.option,s=s,npc=npc, p = p, m=m,
              center=center,pve=pve, r = r, DtD = DtD,
              U = Eigen.new$vectors[,1:npc],G_invhalf = G_invhalf)
  
  class(res) <- "face.sparse"
  return(res)
}


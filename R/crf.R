# Nonparametric estimator for CRF (cumulative rate function)
# data: a list of setting values:
#  t: time point
#  nobs: individual obs
#  c.time: censored time point

fit.crf <- function(sim_dat, newt=NULL, onlyShape=T){
    # info from sim_dat
    c.time <- sim_dat$c.time
    time <- sim_dat$t
    un.time <- unlist(time)
    # precalculation
    time.grid <- sort(unique(un.time))
    time.count <- table(un.time)

    c.temp <- c.time
    t.temp <- time
    n.cnt <- c()
    for(i in seq_along(time.grid)){
        #index for subject whose censor time is smaller than s_l
        idx <- which(time.grid[i] > c.temp)
        if(length(idx)>0){
            c.temp <- c.temp[-idx]
            t.temp <- t.temp[-idx]
        }
        temp <- unlist(t.temp)
        n.cnt <- c(n.cnt, length(temp[temp <= time.grid[i]]))
    }
    prod.vtr <- 1-time.count/n.cnt
    shape <- function(t) prod(prod.vtr[which(time.grid > t)])

    c.shape <- sapply(c.time, shape)
    # Lambda(C_0)
    
    # Lambda(t)
    if (onlyShape){
        crf <- function(t) shape(t)
        crf.fitted = c.shape
    } else {
        nobs <- sapply(time, length)
        crf.censor <- mean(nobs/ c.shape)
        crf.fitted = c.shape*crf.censor
        crf <- function(t) crf.censor *shape(t)
    }

    if (is.null(newt)){
        return(crf.fitted)
    }else{
        return(sapply(newt, crf))
    }
}


get_nhpp_realization <- function(lambda, t_max=10){
    #t_max <- 10
    t <- 0
    Lambda <- function(tupper) integrate(f = lambda, lower = 0, upper = tupper)$value
    Ft <- function(x) Lambda(x)/ Lambda(t_max)
    Ft_inv <- function(u){
        a <- 0
        b <- t_max+2
        eps <- 1e-6
        while(abs(a-b)>eps){
            if(Ft((a + b)/2)<=u) a <- (a+b)/2
            if(Ft((a + b)/2)>u) b <- (a+b)/2
        }
        return(0.5*(a+b))
    }
    n <- rpois(1, Lambda(t_max))
    if (n<2) n = 2
    X <- sapply(1:n, function(z) Ft_inv(runif(1)))
    return(sort(X))
}

#Simulation of Non-Homogeneous Poisson Processes
datf.wang = function(param){
  grid = seq(0.025, 0.975, 0.025)
  n = param[['n']]

  Z <- runif(n,0.5,1.5)

  C <- sapply(1:n, function(i) rtrunc(1, 'exp', rate = 0.1*Z[i], a = 1, b = 10))
  ord <- order(C)
  Z <- Z[ord]
  C <- C[ord]

  t <- lapply(1:n, function(i){
        lambda <- function(t) Z[i]*(0.6 + (t-6)^3/360)
        res_1 <- get_nhpp_realization(lambda, C[i])
        res_1[res_1<=C[i]] # may modify the nhpp to have censor
    })
  m <- sapply(1:n, function(i) length(t[[i]])) # or discard the zeros
  idx <- which(m==0)
  #m = m[-idx]; Z = Z[-idx]; C = C[-idx]
  return(list(t = t, m = m, c.time=C))
}


fit.wang <- function(param){
    a <- datf.wang(param)
    # info from generation
    c.time <- a$c.time
    m <- a$m
    time <- a$t
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
    crf.censor <- mean(m/ c.shape)
    # Lambda(t)
    crf <- function(t) crf.censor *shape(t)
    return(sapply(seq(0,10,0.1), shape))
}


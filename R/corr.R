matern <- function(x,y=NULL,nu=1,rho=1)
{
    if(is.null(y)) y <- x
    
    G <- expand.grid(x,x)
    
    S <- apply(G,1,function(z){
        delta <- abs(z[1]-z[2])/rho
        if(delta == 0)
            1
        else
            (sqrt(2*nu)*delta)^nu * besselK(sqrt(2*nu)*delta,nu=nu) / (2^(nu-1)*gamma(nu))
    })
    
    C <- matrix(S,nrow=length(x),ncol=length(y),byrow=F)
    return(C)
}

C.matern.sig = function(t)  rep(1, length(t))* sqrt(0.6) #*sqrt(t) * exp(-(t-0.5)^2/10+1))
C.matern.stf = function(s,t) C.matern.sig(s) * matern(s,t,nu=0.5,rho=1)* C.matern.sig(t)
C.matern = function(ts){
  sig = C.matern.sig(ts)
  mat = matern(ts,nu=1,rho=1)
  mat = tcrossprod(sig) * mat
  #mat = mat*0   # no sigma_x
  return(mat)
}

sigxf = function(ts){
  mat = matrix(NA,length(ts),length(ts))
  for(i in 1:length(ts)){
    for(j in 1:length(ts)){
      s = ts[i]
      t = ts[j]
      mat[i,j] = (0.5)*exp(s^2) * (0.5)*exp(t^2) * (exp(-(s-t)^2))
    }
  }
  mat
}

# nonperiodic and nonsmooth covariance in Lin
f_rho = function(s,t) exp(-(s-t)^2)
f_sigmoid = function(t) 1/(1+exp(-t))
f_v = function(t){
  sqrt(f_sigmoid(pi*(t-0.5)))
  #sqrt(1/2*(sin(pi*t - pi/2)+1))
}

C.2.coef = function(ts) {
  n=length(ts)
  diff = ts[2]-ts[1]
  temp = matrix(1, nrow=n, ncol=n)
  for (i in seq(n-1)) temp[abs(row(temp) - col(temp)) == i] <- exp(-(i*diff)^2)
  temp
}

C.2= function(ts) {
  2*tcrossprod(sapply(ts, f_v)) * C.2.coef(ts)
}




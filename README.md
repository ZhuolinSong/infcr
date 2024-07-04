# INFormative Functional Concurrent Regression
Faces with object-specific (OBS), subject-specific (SUBJ), and adaptive subject-specific (aSUBJ) weighting.

## Install
devtools::install_github("ZhuolinSong/infcr")

## Demonstration
library(infcr)

set.seed(1)

### Choose mean/covariate functions and generate a functional dataset

x = seq(0.025, 0.975,,100)

mu = function(t) -((t-1)^2+1/3)

f1 = function(t) 1

f2 = function(t) -1

f3 = function(t) 0

param <- list(dist = 'gamma', omg = 1, tau = 2, sigz = 0.5, mu=mu, f1=f1, f2=f2, f3=f3, n = 300, rho = 0.5, C = 4, pp='hpp', grid = x, ca=0.5, cb=2, sige = 0.5)

tru = bayes.true(param)

dat <- bayes.datf(param)

### Estimate semi-parametric model using the aSUBJ weights

fit <- infcr(formula = y ~ s(argvals, k=10, bs="ps") + A1 + A2,
                    argvals = "argvals", subj="subj", argvals.new = x, niter=1,
                    data=dat$data, cov.weight = 'aSUBJ', mean.weight = 'aSUBJ', err.weight = 'aSUBJ', 
                    face.args = list(knots=7, pve=0.99))
        
mu = baseline(fit$fit)

beta = get_coef(fit$fit)

cov = fit$face.object$Chat.new

sigma2 = fit$face.object$sigma2


## References
Song, Z., Xiao, L., Lu, W.(2018). Functional data analysis for longitudinal data with informative observation and censoring times. manuscript

Xiao, L., Li, C., Checkley, W., & Crainiceanu, C. (2018). [Fast covariance estimation for sparse functional data](https://link.springer.com/article/10.1007/s11222-017-9744-8). Statistics and computing, 28, 511-522.

Leroux, A., Xiao, L., Crainiceanu, C., & Checkley, W. (2018). [Dynamic prediction in functional concurrent regression with an application to child growth](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.7582). Statistics in medicine, 37(8), 1376-1388.

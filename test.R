devtools::install_github("ZhuolinSong/infcr", force=T)
devtools::install(force=T)
library(infcr)

x = seq(0.025, 0.975,,100)
mu = function(t) -((t-1)^2+1/3)
f1 = function(t) 1
f2 = function(t) -1
f3 = function(t) 0
param <- list(dist = 'gamma', omg = 1, tau = 2, sigz = 0.5, mu=mu, f1=f1, f2=f2, f3=f3, n = 300, rho = 0.5, C = 4, pp='hpp', grid = x, ca=0.5, cb=2, sige = 0.5)
tru = bayes.true(param)
dat <- bayes.datf(param)
fit <- infcr(formula = y ~ s(argvals, k=10, bs="ps") + A1 + A2,
                    argvals = "argvals", subj="subj", argvals.new = x, niter=1,
                    data=dat$data, cov.weight = 'aSUBJ', mean.weight = 'aSUBJ', err.weight = 'aSUBJ', 
                    face.args = list(knots=7, pve=0.99))
        
mu = baseline(fit$fit)
beta = get_coef(fit$fit)
cov = fit$face.object$Chat.new
sigma2 = fit$face.object$sigma2

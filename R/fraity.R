# ICC distribution
s_z_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    sigz = param[['sigz']]
    mix = param[['mix']]

    if(mix) {
        mix = omg
        omg = 1
    } else {
        mix = 1
    }

    if (dist=='gamma'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param[['logmn']]
        logsd = param[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }

    if (omg==0){ #independent
        return(rep(1,length(ts)))
    }

    if(C==1){
        l_lim = function(ts) {sapply(ts, function(t) max((t-(1-omg))/(omg/tau), 0))}
        u_lim = Inf
        f_p = function(zs, t) {sapply(zs, function(z) {(1 - max((t-omg*min(z,tau)/tau)/(1-omg), 0))*f_den(z)})}
        f_z = function(zs, t) {sapply(zs, function(z) {(rho*z+(1-rho))*(1 - max((t-omg*min(z,tau)/tau)/(1-omg), 0))*f_den(z)})}
        num = sapply(ts, function(t) {integrate(f_z, l_lim(t), u_lim, t=t)$value})
        den = sapply(ts, function(t) {integrate(f_p, l_lim(t), u_lim, t=t)$value})
        (mix*num + (1-mix)*(1-ts)) / (mix*den + (1-mix)*(1-ts))
        #return(kappa*sapply(ts, function(t) {integrate(f_z, l_lim(t), u_lim, t=t)$value / integrate(f_p, l_lim(t), u_lim, t=t)$value}))
    } else if (C==2){
        l_lim = 0
        u_lim = function(ts) {sapply(ts, function(t) if(t>1-omg) tau*omg/(t+omg-1) else Inf)}
        f_p = function(zs, t) {sapply(zs, function(z) {(1 - max((t-omg*tau/max(z,tau))/(1-omg), 0))*f_den(z)})}
        f_z = function(zs, t) {sapply(zs, function(z) {(rho*z+(1-rho))*(1 - max((t-omg*tau/max(z,tau))/(1-omg), 0))*f_den(z)})}
        num = sapply(ts, function(t) {integrate(f_z, l_lim, u_lim(t), t=t)$value})
        den = sapply(ts, function(t) {integrate(f_p, l_lim, u_lim(t), t=t)$value})
        (mix*num + (1-mix)*(1-ts)) / (mix*den + (1-mix)*(1-ts))
        
        #return(kappa*sapply(ts, function(t) {integrate(f_z, l_lim, u_lim(t), t=t)$value / integrate(f_p, l_lim, u_lim(t), t=t)$value}))
    }

}
s_z2_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    sigz = param[['sigz']]
    mix = param[['mix']]

    if(mix) {
        mix = omg
        omg = 1
    } else {
        mix = 1
    }

    if (dist=='gamma'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param[['logmn']]
        logsd = param[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }

    if (omg==0){ #independent
        return(rep(1+rho^2*sigz,length(ts)))
    }
    if(C==1){
        l_lim = function(ts) {sapply(ts, function(t) max((t-(1-omg))/(omg/tau), 0))}
        u_lim = Inf
        f_p = function(zs, t) {sapply(zs, function(z) {(1 - max((t-omg*min(z,tau)/tau)/(1-omg), 0))*f_den(z)})}
        f_z = function(zs, t) {sapply(zs, function(z) {(rho*z+(1-rho))^2*(1 - max((t-omg*min(z,tau)/tau)/(1-omg), 0))*f_den(z)})}
        num = sapply(ts, function(t) {integrate(f_z, l_lim(t), u_lim, t=t)$value})
        den = sapply(ts, function(t) {integrate(f_p, l_lim(t), u_lim, t=t)$value})
        (mix*num + (1-mix)*(1-ts)*(1+rho^2*sigz)) / (mix*den + (1-mix)*(1-ts))
    } else if (C==2){
        l_lim = 0
        u_lim = function(ts) {sapply(ts, function(t) if(t>1-omg) tau*omg/(t+omg-1) else Inf)}
        f_p = function(zs, t) {sapply(zs, function(z) {(1 - max((t-omg*tau/max(z,tau))/(1-omg), 0))*f_den(z)})}
        f_z = function(zs, t) {sapply(zs, function(z) {(rho*z+(1-rho))^2 * (1 - max((t-omg*tau/max(z,tau))/(1-omg), 0))*f_den(z)})}
        num = sapply(ts, function(t) {integrate(f_z, l_lim, u_lim(t), t=t)$value})
        den = sapply(ts, function(t) {integrate(f_p, l_lim, u_lim(t), t=t)$value})
        (mix*num + (1-mix)*(1-ts)*(1+rho^2*sigz)) / (mix*den + (1-mix)*(1-ts))
    }
}
s_cinv_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    sigz = param[['sigz']]
    mix = param[['mix']]

    if(mix) {
        mix = omg
        omg = 1
    } else {
        mix = 1
    }

    if (dist=='gamma'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param[['logmn']]
        logsd = param[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }

    if(C==1){
        l_lim = function(ts) {sapply(ts, function(t) max((t-(1-omg))/(omg/tau), 0))}
        u_lim = Inf
        if(omg != 1) cinv = function(z, t)  sapply(z, function(z) (log(omg*min(z,tau)/tau+(1-omg)) - log(omg*min(z,tau)/tau+(1-omg)*max((t-omg*min(z,tau)/tau)/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) cinv = function(z, t)  sapply(z, function(z) (tau/min(z,tau))*f_den(z))
        return(mix*sapply(ts, function(t) integrate(cinv, l_lim(t), u_lim, t=t)$value)-(1-mix)*log(ts))
        #sapply(ts, function(t) integrate(cinv, l_lim(t), u_lim, t=t)$value)
    } else if (C==2){
        l_lim = 0
        u_lim = function(ts) {sapply(ts, function(t) if(t>1-omg)tau*omg/(t+omg-1) else Inf)}
        if(omg != 1) cinv = function(z, t)  sapply(z, function(z) (log(omg*tau/max(z,tau)+(1-omg)) - log(omg*tau/max(z,tau)+(1-omg)*max((t-omg*tau/max(z,tau))/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) cinv = function(z, t)  sapply(z, function(z) (max(z,tau)/tau)*f_den(z))
        return(mix*sapply(ts, function(t) integrate(cinv, l_lim, u_lim(t), t=t)$value)-(1-mix)*log(ts))
        #sapply(ts, function(t) integrate(cinv, l_lim, u_lim(t), t=t)$value)
    }

}
s_zcinv_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    sigz = param[['sigz']]
    mix = param[['mix']]

    if(mix) {
        mix = omg
        omg = 1
    } else {
        mix = 1
    }

    if (dist=='gamma'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param[['logmn']]
        logsd = param[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }

    if(C==1){
        l_lim = function(ts) {sapply(ts, function(t) max((t-(1-omg))/(omg/tau), 0))}
        u_lim = Inf
        if(omg != 1) zcinv = function(z, t)  sapply(z, function(z) (rho*z*kappa+(1-rho)) * (log(omg*min(z,tau)/tau+(1-omg)) - log(omg*min(z,tau)/tau+(1-omg)*max((t-omg*min(z,tau)/tau)/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) zcinv = function(z, t)  sapply(z, function(z) (rho*z*kappa+(1-rho)) * (tau/min(z,tau))*f_den(z))
        mix*sapply(ts, function(t) integrate(zcinv, l_lim(t), u_lim, t=t)$value)-(1-mix)*log(ts)
    } else if (C==2){
        l_lim = 0
        u_lim = function(ts) {sapply(ts, function(t) if(t>1-omg)tau*omg/(t+omg-1) else Inf)}
        if(omg != 1) zcinv = function(z, t)  sapply(z, function(z) (rho*z*kappa+(1-rho)) * (log(omg*tau/max(z,tau)+(1-omg)) - log(omg*tau/max(z,tau)+(1-omg)*max((t-omg*tau/max(z,tau))/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) zcinv = function(z, t)  sapply(z, function(z) (rho*z*kappa+(1-rho)) * (max(z,tau)/tau)*f_den(z))
        mix*sapply(ts, function(t) integrate(zcinv, l_lim, u_lim(t), t=t)$value)-(1-mix)*log(ts)
    }

}


# Bayes GEN
b_z_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b, lower=F) / pgamma(tau*c, a, b, lower=F)) + (1-rho)
        }
    } else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        Ez_given_c = function(c){
            rho * extrunc('lnorm', 0, tau*c, meanlog=mn, sdlog=sd) + (1-rho)
        }
    }


    if(C==1){
        f_c = function(cs, a, b, ca, cb, tau) {
            tmp = function(z, c) dgamma(z, a, b) * dbeta(c, ca+z*tau, cb)
            sapply(cs, function(c) integrate(tmp, 0, Inf, c=c)$value)
        }

        f_joint_cgt = function(z, t, a, b, ca, cb, tau) {
            tmp = function(c) dgamma(z, a, b) * dbeta(c, ca+z*tau, cb)
            if(t==1) return(0)
            integrate(tmp, t, 1)$value /integrate(f_c, t, 1, a=a,b=b,ca=ca,cb=cb,tau=tau)$value
        }

        z_cgt = function(zs, t, a, b, ca, cb, tau) {
            zs*sapply(zs, f_joint_cgt, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)
        }
        
        Ez_cgt = function(ts, a, b, ca, cb, tau) {
            sapply(ts, function(t) integrate(z_cgt,0,Inf, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)$value)
        }

        return(Ez_cgt(ts, a, b, ca, cb, tau))
    } else if (C==2){
        f_c = function(cs, a, b, ca, cb, tau) {
            tmp = function(z, c) dgamma(z, a, b) * dbeta(c, ca, cb+z*tau)
            sapply(cs, function(c) integrate(tmp, 0, Inf, c=c)$value)
        }
        f_joint_cgt = function(z, t, a, b, ca, cb, tau) {
            tmp = function(c) dgamma(z, a, b) * dbeta(c, ca, cb+z*tau)
            if(t==1) return(0)
            integrate(tmp, t, 1)$value /integrate(f_c, t, 1, a=a,b=b,ca=ca,cb=cb,tau=tau)$value
        }
        z_cgt = function(zs, t, a, b, ca, cb, tau) {
            zs*sapply(zs, f_joint_cgt, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)
        }
        Ez_cgt = function(ts, a, b, ca, cb, tau) {
            sapply(ts, function(t) integrate(z_cgt,0,Inf, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)$value)
        }
        return(Ez_cgt(ts, a, b, ca, cb, tau))
    } else if (C==3){
        Ez_cgt = function(ts, a, b, ca, tau) {
            if(omg==1) a/(b+ca*ts*tau)
            else a*omg/(1-omg)/ca/tau*(log(omg*b+ca*tau*ts)-log(omg*b+ca*tau*(ts-1+omg)))
        }
        return(rho*sapply(ts, Ez_cgt, a=a,b=b,ca=ca,tau=tau) + (1-rho))
    } else if (C==4){
        Ez_cgt = function(ts, a, b, ca, tau) {
            k=tau/ts
            a*((b+ca*k)^a/b - b^a/(b+ca*k)) / ((b+ca*k)^a-b^a)
        }
        return(rho*sapply(ts, Ez_cgt, a=a,b=b,ca=ca,tau=tau) + (1-rho))
    }
}
b_z2_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b, lower=F) / pgamma(tau*c, a, b, lower=F)) + (1-rho)
        }
    } else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        Ez_given_c = function(c){
            rho * extrunc('lnorm', 0, tau*c, meanlog=mn, sdlog=sd) + (1-rho)
        }
    }


    if(C==1){
        f_c = function(cs, a, b, ca, cb, tau) {
            tmp = function(z, c) dgamma(z, a, b) * dbeta(c, ca+z*tau, cb)
            sapply(cs, function(c) integrate(tmp, 0, Inf, c=c)$value)
        }

        f_joint_cgt = function(z, t, a, b, ca, cb, tau) {
            tmp = function(c) dgamma(z, a, b) * dbeta(c, ca+z*tau, cb)
            if(t==1) return(0)
            integrate(tmp, t, 1)$value /integrate(f_c, t, 1, a=a,b=b,ca=ca,cb=cb,tau=tau)$value
        }

        z2_cgt = function(zs, t, a, b, ca, cb, tau) {
            zs^2*sapply(zs, f_joint_cgt, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)
        }
        
        Ez2_cgt = function(ts, a, b, ca, cb, tau) {
            sapply(ts, function(t) integrate(z2_cgt,0,Inf, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)$value)
        }

        return(Ez2_cgt(ts, a, b, ca, cb, tau) )
    } else if (C==2){
        f_c = function(cs, a, b, ca, cb, tau) {
            tmp = function(z, c) dgamma(z, a, b) * dbeta(c, ca, cb+z*tau)
            sapply(cs, function(c) integrate(tmp, 0, Inf, c=c)$value)
        }
        f_joint_cgt = function(z, t, a, b, ca, cb, tau) {
            tmp = function(c) dgamma(z, a, b) * dbeta(c, ca, cb+z*tau)
            if(t==1) return(0)
            integrate(tmp, t, 1)$value /integrate(f_c, t, 1, a=a,b=b,ca=ca,cb=cb,tau=tau)$value
        }
        z2_cgt = function(zs, t, a, b, ca, cb, tau) {
            zs^2*sapply(zs, f_joint_cgt, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)
        }
        Ez2_cgt = function(ts, a, b, ca, cb, tau) {
            sapply(ts, function(t) integrate(z2_cgt,0,Inf, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)$value)
        }
        return(Ez2_cgt(ts, a, b, ca, cb, tau) )
    } else if (C==3){
        Ez_cgt = function(ts, a, b, ca, tau) {
            if(omg==1) a/(b+ca*ts*tau)
            else a*omg/(1-omg)/ca/tau*(log(omg*b+ca*tau*ts)-log(omg*b+ca*tau*(ts-1+omg)))
        }
        Ez2_cgt = function(ts, a, b, ca, tau) {
            if (omg==1) (a+1)*a / (b+ca*ts*tau)^2
            else a*(a+1)/(ca*tau*(1-omg))*(1/(b+ca*tau*(ts-1+omg)/omg) - 1/(b+ca*tau*ts/omg))
        }
        return(rho^2*Ez2_cgt(ts,a,b,ca,tau)+ 2*rho*(1-rho)*Ez_cgt(ts,a,b,ca,tau) + (1-rho)^2)
    } else if (C==4){
        Ez_cgt = function(ts, a, b, ca, tau) {
            k=tau/ts
            a*((b+ca*k)^a/b - b^a/(b+ca*k)) / ((b+ca*k)^a-b^a)
        }
        Ez2_cgt = function(ts, a, b, ca, tau) {
            k=tau/ts
            (a+1)*a*((b+ca*k)^a/b^2 - b^a/(b+ca*k)^2) / ((b+ca*k)^a-b^a)
        }
        return(rho^2*Ez2_cgt(ts,a,b,ca,tau)+ 2*rho*(1-rho)*Ez_cgt(ts,a,b,ca,tau) + (1-rho)^2)
    }
}
b_cinv_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b, lower=F) / pgamma(tau*c, a, b, lower=F)) + (1-rho)
        }
    } else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        Ez_given_c = function(c){
            rho * extrunc('lnorm', 0, tau*c, meanlog=mn, sdlog=sd) + (1-rho)
        }
    }


    if(C==1){
        f_c = function(cs, a, b, ca, cb, tau) {
            tmp = function(z, c) dgamma(z, a, b) * dbeta(c, ca+z*tau, cb)
            sapply(cs, function(c) integrate(tmp, 0, Inf, c=c)$value)
        }
        
        Einvc = function(ts, a, b, ca, cb, tau) {
            tmp = function(cs) 1/cs*f_c(cs, a,b,ca,cb,tau)
            sapply(ts, function(t) integrate(tmp, t, 1)$value)
        }

        return(Einvc(ts, a, b, ca, cb, tau) )
    } else if (C==2){

        f_c = function(cs, a, b, ca, cb, tau) {
            tmp = function(z, c) dgamma(z, a, b) * dbeta(c, ca, cb+z*tau)
            sapply(cs, function(c) integrate(tmp, 0, Inf, c=c)$value)
        }

        Einvc = function(ts, a, b, ca, cb, tau) {
            tmp = function(cs) 1/cs*f_c(cs, a,b,ca,cb,tau)
            sapply(ts, function(t) integrate(tmp, t, 1)$value)
        }

        return(Einvc(ts, a, b, ca, cb, tau) )
    } else if (C==3){
        f_joint = function(z, c, a, b, ca) {
            dgamma(z, a, b) * dexp(c, ca*z)
        }
        Ecinv_given_z = function(z, t, a, b, ca, tau) {
            tmp = function(cs) sapply(cs, function(c) tau/min(c,tau)*f_joint(z,c, a,b,ca))
            integrate(tmp, t*tau, Inf)$value
        }
        Ecinv = function(ts, a,b,ca,tau){
            tmp = function(zs, t) sapply(zs, Ecinv_given_z, t=t,a=a,b=b,ca=ca,tau=tau)
            sapply(ts, function(t) integrate(tmp, 0, Inf, t=t)$value)
        }

        return(Ecinv(ts, a, b, ca, tau))
    } else if (C==4){
        f_joint = function(z, c, a, b, ca) {
            dgamma(z, a, b) * dexp(c, ca*z)
        }
        Ecinv_given_z = function(z, t, a, b, ca, tau) {
            tmp = function(cs) sapply(cs, function(c) max(c,tau)/tau*f_joint(z,c, a,b,ca))
            integrate(tmp, 0, tau/t)$value
        }
        Ecinv = function(ts, a,b,ca,tau){
            tmp = function(zs, t) sapply(zs, Ecinv_given_z, t=t,a=a,b=b,ca=ca,tau=tau)
            sapply(ts, function(t) integrate(tmp, 0, Inf, t=t)$value)
        }

        return(Ecinv(ts, a, b, ca, tau))
    }
}
b_zcinv_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b, lower=F) / pgamma(tau*c, a, b, lower=F)) + (1-rho)
        }
    } else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        Ez_given_c = function(c){
            rho * extrunc('lnorm', 0, tau*c, meanlog=mn, sdlog=sd) + (1-rho)
        }
    }


    if(C==1){
        f_joint = function(z, c, a, b, ca, cb, tau) {
            dgamma(z, a, b) * dbeta(c, ca+z*tau, cb)
        }
        Ezcinv_given_z = function(z, t, a, b, ca, cb, tau) {
            tmp = function(c) z/c*f_joint(z, c, a,b,ca,cb,tau)
            integrate(tmp, t, 1)$value
        }
        Ezcinv = function(ts, a,b,ca,cb,tau){
            tmp = function(zs, t) sapply(zs, Ezcinv_given_z, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)
            sapply(ts, function(t) integrate(tmp, 0, Inf, t=t)$value)
        }
        return(Ezcinv(ts, a, b, ca, cb, tau) )
    } else if (C==2){
        f_joint = function(z, c, a, b, ca, cb, tau) {
            dgamma(z, a, b) * dbeta(c, ca, cb+z*tau)
        }
        Ezcinv_given_z = function(z, t, a, b, ca, cb, tau) {
            tmp = function(c) z/c*f_joint(z, c, a,b,ca,cb,tau)
            integrate(tmp, t, 1)$value
        }
        Ezcinv = function(ts, a,b,ca,cb,tau){
            tmp = function(zs, t) sapply(zs, Ezcinv_given_z, t=t,a=a,b=b,ca=ca,cb=cb,tau=tau)
            sapply(ts, function(t) integrate(tmp, 0, Inf, t=t)$value)
        }
        return(Ezcinv(ts, a, b, ca, cb, tau) )
    } else if (C==3){
        f_joint = function(z, c, a, b, ca) {
            dgamma(z, a, b) * dexp(c, ca*z)
        }
        Ezcinv_given_z = function(z, t, a, b, ca, tau) {
            tmp = function(cs) sapply(cs, function(c) (rho*z+(1-rho)) * tau/min(c,tau)*f_joint(z,c, a,b,ca))
            integrate(tmp, t*tau, Inf)$value
        }
        Ezcinv = function(ts, a,b,ca,tau){
            tmp = function(zs, t) sapply(zs, Ezcinv_given_z, t=t,a=a,b=b,ca=ca,tau=tau)
            sapply(ts, function(t) integrate(tmp, 0, Inf, t=t)$value)
        }
        return(Ezcinv(ts, a, b, ca, tau) )
    }  else if (C==4){
        f_joint = function(z, c, a, b, ca) {
            dgamma(z, a, b) * dexp(c, ca*z)
        }
        Ezcinv_given_z = function(z, t, a, b, ca, tau) {
            tmp = function(cs) sapply(cs, function(c) (rho*z+(1-rho)) * max(c,tau)/tau*f_joint(z,c, a,b,ca))
            integrate(tmp, 0, tau/t)$value
        }
        Ezcinv = function(ts, a,b,ca,tau){
            tmp = function(zs, t) sapply(zs, Ezcinv_given_z, t=t,a=a,b=b,ca=ca,tau=tau)
            sapply(ts, function(t) integrate(tmp, 0, Inf, t=t)$value)
        }
        return(Ezcinv(ts, a, b, ca, tau) )
    }
}


# NEW GEN
z_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b) / pgamma(tau*c, a, b)) + (1-rho)
        }   
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b, lower=F) / pgamma(tau*c, a, b, lower=F)) + (1-rho)
        }
    } else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        Ez_given_c = function(c){
            rho * extrunc('lnorm', 0, tau*c, meanlog=mn, sdlog=sd) + (1-rho)
        }
    }else if (dist=='beta'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param[['a']]
        b = param[['b']]
        c = param[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }


    if(C==1){
        l_limit = function(tt) {sapply(tt, function(t) max((t-(1-omg))/omg, 0))}
        u_limit = 1
        p_cgt = function(cs, t) {sapply(cs, function(c) {(1 - max((t-omg*c)/(1-omg), 0))*dbeta(c, ca, cb)})}
        z_cgt = function(cs, t) {sapply(cs, function(c) {Ez_given_c(c)*(1 - max((t-omg*c)/(1-omg), 0))*dbeta(c, ca, cb)})}
        return(sapply(ts, function(t) {integrate(z_cgt, l_limit(t), u_limit, t=t)$value / integrate(p_cgt, l_limit(t), u_limit, t=t)$value}))
    } else if (C==2){
        l_limit = 0
        u_limit = function(tt) {sapply(tt, function(t) min((1-t)/omg, 1))}
        p_cgt = function(cs, t) {sapply(cs, function(c) {(1 - max((t-omg*(1-c))/(1-omg), 0))*dbeta(c, ca, cb)})}
        z_cgt = function(cs, t) {sapply(cs, function(c) {Ez_given_c(c)*(1 - max((t-omg*(1-c))/(1-omg), 0))*dbeta(c, ca, cb)})}
        return(sapply(ts, function(t) {integrate(z_cgt, l_limit, u_limit(t), t=t)$value / integrate(p_cgt, l_limit, u_limit(t), t=t)$value}))
    }
}
z2_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']] 
        Ez2_given_c = function(c){
            a/b*(rho^2*(a+1)/b*pgamma(tau*c, a+2, b) + 2*rho*(1-rho)*pgamma(tau*c, a+1, b))/pgamma(tau*c, a, b) + (1-rho)^2
        }
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez2_given_c = function(c){
            a/b*(rho^2*(a+1)/b*pgamma(tau*c, a+2, b, lower=F) + 2*rho*(1-rho)*pgamma(tau*c, a+1, b, lower=F))/pgamma(tau*c, a, b, lower=F) + (1-rho)^2
        }
    }else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        z2_given_c = function(z) {
            (rho*z+(1-rho))^2 * dlnorm(z,mn,sd) 
        }
        Ez2_given_c = function(c){
            integrate(z2_given_c, 0, tau*c)$value / plnorm(tau*c,mn,sd)
        }
    }else if (dist=='beta'){
        a = param[['alpha']]
        b = param[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param[['a']]
        b = param[['b']]
        c = param[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }


    if(C==1){
        l_limit = function(tt) {sapply(tt, function(t) max((t-(1-omg))/omg, 0))}
        u_limit = 1
        p_cgt = function(cs, t) {sapply(cs, function(c) {(1 - max((t-omg*c)/(1-omg), 0))*dbeta(c, ca, cb)})}
        z2_cgt = function(cs, t) {sapply(cs, function(c) {Ez2_given_c(c)*(1 - max((t-omg*c)/(1-omg), 0))*dbeta(c, ca, cb)})}
        return(sapply(ts, function(t) {integrate(z2_cgt, l_limit(t), u_limit, t=t)$value / integrate(p_cgt, l_limit(t), u_limit, t=t)$value}))
    } else if (C==2){
        l_limit = 0
        u_limit = function(tt) {sapply(tt, function(t) min((1-t)/omg, 1))}
        p_cgt = function(cs, t) {sapply(cs, function(c) {(1 - max((t-omg*(1-c))/(1-omg), 0))*dbeta(c, ca, cb)})}
        z2_cgt = function(cs, t) {sapply(cs, function(c) {Ez2_given_c(c)*(1 - max((t-omg*(1-c))/(1-omg), 0))*dbeta(c, ca, cb)})}
        return(sapply(ts, function(t) {integrate(z2_cgt, l_limit, u_limit(t), t=t)$value / integrate(p_cgt, l_limit, u_limit(t), t=t)$value}))
    }
}
cinv_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
    }else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
    }else if (dist=='beta'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param[['a']]
        b = param[['b']]
        c = param[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }

    if(C==1){
        l_limit = function(tt) {sapply(tt, function(t) max((t-(1-omg))/omg, 0))}
        u_limit = 1
        if(omg != 1) cinv_cgt = function(cs, t)  sapply(cs, function(c) (log(omg*c+(1-omg)) - log(omg*c+(1-omg)*max((t-omg*c)/(1-omg), 0)))*dbeta(c, ca, cb))/(1-omg)
        if(omg == 1) cinv_cgt = function(cs, t)  sapply(cs, function(c) (1/c)*dbeta(c, ca, cb))
        return(sapply(ts, function(t) {integrate(cinv_cgt, l_limit(t), u_limit, t=t)$value}))
    } else if (C==2){
        l_limit = 0
        u_limit = function(tt) {sapply(tt, function(t) min((1-t)/omg, 1))}
        if(omg != 1) cinv_cgt = function(cs, t) sapply(cs, function(c) (log(omg*(1-c)+(1-omg)) - log(omg*(1-c)+(1-omg)*max((t-omg*(1-c))/(1-omg), 0)))*dbeta(c, ca, cb))/(1-omg)
        if(omg == 1) cinv_cgt = function(cs, t) sapply(cs, function(c) (1/(1-c))*dbeta(c, ca, cb))
        return(sapply(ts, function(t) {integrate(cinv_cgt, l_limit, u_limit(t), t=t)$value}))
    }
}
zcinv_cgt = function(ts, param){
    rho = param[['rho']]
    dist = param[['dist']]
    C = param[['C']]
    tau = param[['tau']]
    kappa = param[['kappa']]
    omg = param[['omg']]
    ca = param[['ca']]
    cb = param[['cb']]

    if (dist=='gamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b) / pgamma(tau*c, a, b)) + (1-rho)
        }
    }else if (dist=='tgamma'){
        a = param[['alpha']]
        b = param[['beta']]
        Ez_given_c = function(c){
            rho* (a/b * pgamma(tau*c, a+1, b, lower=F) / pgamma(tau*c, a, b, lower=F)) + (1-rho)
        }
    } else if (dist=='lognormal'){
        mn = param[['logmn']]
        sd = param[['logsd']]
        Ez_given_c = function(c){
            rho * extrunc('lnorm', 0, tau*c, meanlog=mn, sdlog=sd) + (1-rho)
        }
    }else if (dist=='beta'){
        alpha = param[['alpha']]
        beta = param[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param[['a']]
        b = param[['b']]
        c = param[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }

    if(C==1){
        l_limit = function(tt) {sapply(tt, function(t) max((t-(1-omg))/omg, 0))}
        u_limit = 1
        if(omg != 1) zcinv_cgt = function(cs, t) sapply(cs, function(c) Ez_given_c(c)*(log(omg*c+(1-omg)) - log(omg*c+(1-omg)*max((t-omg*c)/(1-omg), 0)))*dbeta(c, ca, cb))/(1-omg)
        if(omg == 1) zcinv_cgt = function(cs, t) sapply(cs, function(c) Ez_given_c(c)*(1/c)*dbeta(c, ca, cb))
        return(sapply(ts, function(t) {integrate(zcinv_cgt, l_limit(t), u_limit, t=t)$value}))
    } else if (C==2){
        l_limit = 0
        u_limit = function(tt) {sapply(tt, function(t) min((1-t)/omg, 1))}
        if(omg != 1) zcinv_cgt = function(cs, t) sapply(cs, function(c) Ez_given_c(c)*(log(omg*(1-c)+(1-omg)) - log(omg*(1-c)+(1-omg)*max((t-omg*(1-c))/(1-omg), 0)))*dbeta(c, ca, cb))/(1-omg)
        if(omg == 1) zcinv_cgt = function(cs, t) sapply(cs, function(c) Ez_given_c(c)*(1/(1-c))*dbeta(c, ca, cb))
        return(sapply(ts, function(t) {integrate(zcinv_cgt, l_limit, u_limit(t), t=t)$value}))
    }
}



# ICC distribution
# Z conditional expectation
z_cond = function(ts, param1, param2){
    rho = param1[['rho']]
    dist = param1[['dist']]
    C = param1[['C']]
    tau = param1[['tau']]
    kappa = param1[['kappa']]
    omg = param1[['omg']]

    if (dist=='gamma'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param2[['logmn']]
        logsd = param2[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }else if (dist=='beta'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param2[['a']]
        b = param2[['b']]
        c = param2[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }

    if (omg==0){ #independent
        return(rep(1,length(ts)))
    }

    if(C==1){
        l_temp = function(tt) {sapply(tt, function(t) max((t-(1-omg))/(omg/tau), 0))}
        u_temp = tau
        f_ptemp = function(zs, tt) {sapply(zs, function(z) {(1 - max((tt-omg*z/tau)/(1-omg), 0))*f_den(z)})}
        f_ztemp = function(zs, tt) {sapply(zs, function(z) {(rho*z+(1-rho))*(1 - max((tt-omg*z/tau)/(1-omg), 0))*f_den(z)})}
        return(kappa*sapply(ts, function(t) {integrate(f_ztemp, l_temp(t), u_temp, tt=t)$value / integrate(f_ptemp, l_temp(t), u_temp, tt=t)$value}))
    } else if (C==2){
        l_temp = 0
        u_temp = function(tt) {sapply(tt, function(t) min((1-t)/(omg/tau), tau))}
        f_ptemp = function(zs, tt) {sapply(zs, function(z) {(1 - max((tt-omg*(1-z/tau))/(1-omg), 0))*f_den(z)})}
        f_ztemp = function(zs, tt) {sapply(zs, function(z) {(rho*z+(1-rho))*(1 - max((tt-omg*(1-z/tau))/(1-omg), 0))*f_den(z)})}
        return(kappa*sapply(ts, function(t) {integrate(f_ztemp, l_temp, u_temp(t), tt=t)$value / integrate(f_ptemp, l_temp, u_temp(t), tt=t)$value}))
    }else if (C==4){
        return(1)
    }else if (C==5){
        cinv = function(z) (1+z)*f_den(z)
        temp = function(z) (1/z - 1)
        sapply(ts, function(x) integrate(cinv, lower = 0, upper = temp(x))$value / f_p(temp(x)))
    }else if (C==6){
        cinv = function(z) (1+z)/z*f_den(z)
        temp = function(z) (1/(1-z) - 1)
        sapply(t, function(x) integrate(cinv, lower = temp(x), upper = Inf)$value / f_p(temp(x), lower=F))
    }

}
# Z^2 conditional expectation
z2_cond = function(ts, param1, param2){
    rho = param1[['rho']]
    dist = param1[['dist']]
    C = param1[['C']]
    tau = param1[['tau']]
    kappa = param1[['kappa']]
    omg = param1[['omg']]
    sigz = param1[['sigz']]

    if (dist=='gamma'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param2[['logmn']]
        logsd = param2[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }else if (dist=='beta'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param2[['a']]
        b = param2[['b']]
        c = param2[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }

    if (omg==0){ #independent
        return(rep(1+rho^2*sigz,length(ts)))
    }

    if(C==1){
        l_temp = function(tt) {sapply(tt, function(t) max((t-(1-omg))/(omg/tau), 0))}
        u_temp = tau
        f_ptemp = function(zs, tt) {sapply(zs, function(z) {(1 - max((tt-omg*z/tau)/(1-omg), 0))*f_den(z)})}
        f_ztemp = function(zs, tt) {sapply(zs, function(z) {(rho*z+(1-rho))^2*(1 - max((tt-omg*z/tau)/(1-omg), 0))*f_den(z)})}
        return(kappa^2*sapply(ts, function(t) {integrate(f_ztemp, l_temp(t), u_temp, tt=t)$value / integrate(f_ptemp, l_temp(t), u_temp, tt=t)$value}))
    } else if (C==2){
        l_temp = 0
        u_temp = function(tt) {sapply(tt, function(t) min((1-t)/(omg/tau), tau))}
        f_ptemp = function(zs, tt) {sapply(zs, function(z) {(1 - max((tt-omg*(1-z/tau))/(1-omg), 0))*f_den(z)})}
        f_ztemp = function(zs, tt) {sapply(zs, function(z) {(rho*z+(1-rho))^2*(1 - max((tt-omg*(1-z/tau))/(1-omg), 0))*f_den(z)})}
        return(kappa^2*sapply(ts, function(t) {integrate(f_ztemp, l_temp, u_temp(t), tt=t)$value / integrate(f_ptemp, l_temp, u_temp(t), tt=t)$value}))
    }else if (C==4){
        return(1)
    }else if (C==5){
        cinv = function(z) (1+z)*f_den(z)
        temp = function(z) (1/z - 1)
        sapply(t, function(x) integrate(cinv, lower = 0, upper = temp(x))$value / f_p(temp(x)))
    }else if (C==6){
        cinv = function(z) (1+z)/z*f_den(z)
        temp = function(z) (1/(1-z) - 1)
        sapply(t, function(x) integrate(cinv, lower = temp(x), upper = Inf)$value / f_p(temp(x), lower=F))
    }

}
# C^-1 conditional expectation
cinv_cond = function(t, param1, param2){
    rho = param1[['rho']]
    dist = param1[['dist']]
    C = param1[['C']]
    tau = param1[['tau']]
    omg = param1[['omg']]

    if (dist=='gamma'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param2[['logmn']]
        logsd = param2[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }else if (dist=='beta'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param2[['a']]
        b = param2[['b']]
        c = param2[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }

    if(C==1){
        l_temp = function(z) sapply(z, function(x) max((x-(1-omg))/omg*tau, 0))
        u_temp = tau
        if(omg != 1) cinv = function(z, t)  sapply(z, function(z) (log(omg*z/tau+(1-omg)) - log(omg*z/tau+(1-omg))*max((t-omg*z/tau)/(1-omg), 0))*f_den(z))/(1-omg)
        if(omg == 1) cinv = function(z, t)  sapply(z, function(z) (tau/z)*f_den(z))
        sapply(t, function(x) integrate(cinv, l_temp(x), u_temp, t=x)$value)
    }else if (C==2){
        l_temp = 0
        u_temp = function(z) sapply(z, function(x) min((1-x)/omg*tau, tau))
        if(omg != 1) cinv = function(z) sapply(z, function(z) (log(omg*(1-z/tau)+(1-omg)) - log(omg*(1-z/tau)+(1-omg)*max((t-omg*(1-z/tau))/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) cinv = function(z) sapply(z, function(z) (1/(1-z/tau))*f_den(z))
        sapply(t, function(x) integrate(cinv, l_temp, u_temp(x))$value)
    }else if (C==4){
        return(1)
    }else if (C==5){
        cinv = function(z) (1+z)*f_den(z)
        temp = function(z) (1/z - 1)
        sapply(t, function(x) integrate(cinv, lower = 0, upper = temp(x))$value / f_p(temp(x)))
    }else if (C==6){
        cinv = function(z) (1+z)/z*f_den(z)
        temp = function(z) (1/(1-z) - 1)
        sapply(t, function(x) integrate(cinv, lower = temp(x), upper = Inf)$value / f_p(temp(x), lower=F))
    }

}
# ZC^-1 conditional expectation
zcinv_cond = function(t, param1, param2){
    rho = param1[['rho']]
    dist = param1[['dist']]
    C = param1[['C']]
    tau = param1[['tau']]
    kappa = param1[['kappa']]
    omg = param1[['omg']]

    if (dist=='gamma'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dgamma(z, shape=alpha, rate=beta)
        f_p = function(z, lower = T) pgamma(z, shape=alpha, rate=beta, lower=lower)
    }else if (dist=='lognormal'){
        logmn = param2[['logmn']]
        logsd = param2[['logsd']]
        f_den = function(z) dlnorm(z, meanlog=logmn, sdlog=logsd)
        f_p = function(z, lower = T) plnorm(z, meanlog=logmn, sdlog=logsd, lower=lower)
    }else if (dist=='beta'){
        alpha = param2[['alpha']]
        beta = param2[['beta']]
        f_den = function(z) dbeta(z, alpha, beta)
        f_p = function(z, lower = T) pbeta(z, alpha, beta, lower=lower)
    }else if (dist=='triangular'){
        a = param2[['a']]
        b = param2[['b']]
        c = param2[['c']]
        f_den = function(z) dtriangle(z, a=a, b=b, c=c)
        f_p = function(z, lower = T) {
            if (lower) return(ptriangle(z, a=a, b=b, c=c))
            else return(1-ptriangle(z, a=a, b=b, c=c))
        }
    }

    if(C==1){
        l_temp = function(z) sapply(z, function(x) max((x-(1-omg))/omg*tau, 0))
        u_temp = tau
        if(omg != 1) cinv = function(z) sapply(z, function(z) (rho*z*kappa + (1-rho))*(log(omg*z/tau+(1-omg))-log(omg*z/tau+(1-omg)*max((t-omg*z/tau)/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) cinv = function(z) sapply(z, function(z) (rho*z*kappa + (1-rho))*(tau/z)*f_den(z))
        sapply(t, function(x) integrate(cinv, l_temp(x), u_temp)$value)
    }else if (C==2){
        l_temp = 0
        u_temp = function(z) sapply(z, function(x) min((1-x)/omg*tau, tau))
        if(omg != 1) cinv = function(z) sapply(z, function(z) (rho*z*kappa + (1-rho))*(log(omg*(1-z/tau)+(1-omg)) - log(omg*(1-z/tau)+(1-omg)*max((t-omg*(1-z/tau))/(1-omg), 0)))*f_den(z))/(1-omg)
        if(omg == 1) cinv = function(z) sapply(z, function(z) (rho*z*kappa + (1-rho))*(1/(1-z/tau))*f_den(z))
        sapply(t, function(x) integrate(cinv, l_temp, u_temp(x))$value)
    }else if (C==4){
        return(1)
    }else if (C==5){
        cinv = function(z) (rho*z+(1-rho))*(1+z)*f_den(z)
        temp = function(z) (1/z - 1)
        sapply(t, function(x) integrate(cinv, lower = 0, upper = temp(x))$value / f_p(temp(x)))
    }else if (C==6){
        cinv = function(z) (rho*z+(1-rho))*(1+z)/z*f_den(z)
        temp = function(z) (1/(1-z) - 1)
        sapply(t, function(x) integrate(cinv, lower = temp(x), upper = Inf)$value / f_p(temp(x), lower=F))
    }

}

# Z expectation for lognormal distribution
Zlog = function(z, logmn, logsd) z*dlnorm(z, meanlog=logmn, sdlog=logsd)
zlogcond = function(t, logmn, logsd, tau=Inf, lower.tail = T){
    t[t>9999] = Inf
    if (lower.tail){
        sapply(t, function(x) integrate(Zlog, lower = 0, upper = x, logmn = logmn, logsd = logsd)$value / plnorm(x, meanlog=logmn, sdlog=logsd))
    } else{
        sapply(t, function(x) integrate(Zlog, lower = x, upper = tau, logmn = logmn, logsd = logsd)$value / (plnorm(tau, meanlog=logmn, sdlog=logsd)-plnorm(x, meanlog=logmn, sdlog=logsd)))
}}
# Z^2 expectation for lognormal distribution
Z2log = function(z, logmn, logsd, tau=Inf) z^2*dlnorm(z, meanlog=logmn, sdlog=logsd)
z2logcond = function(t, logmn, logsd, tau=Inf, lower.tail = T){
    t[t>9999] = Inf
    if (lower.tail){
        sapply(t, function(x) integrate(Z2log, lower = 0, upper = x, logmn = logmn, logsd = logsd)$value / plnorm(x, meanlog=logmn, sdlog=logsd))
    } else{
        sapply(t, function(x) integrate(Z2log, lower = x, upper = tau, logmn = logmn, logsd = logsd)$value / (plnorm(tau, meanlog=logmn, sdlog=logsd)-plnorm(x, meanlog=logmn, sdlog=logsd)))
}}


# C conditional probability
dCZ = function(z, c = 0.5, shape = 1, rate = 1, lambda = 0.1) {
    sapply(z, function(z) truncdist::dtrunc(c, spec = 'exp', b = 10, rate = lambda*z) * dgamma(z, shape = shape, rate = rate))
}
dC = function(c, shape=1, rate=1, lambda=0.1){
    sapply(c, function(c) integrate(dCZ, lower = 0, upper = Inf, c = c, shape = shape, rate = rate, lambda = lambda)$value)
}
pC = function(t, shape=1, rate=1, lambda=0.1){
    sapply(t, function(t) integrate(dC, lower = 10*t, upper = 10, shape = shape, rate = rate, lambda = lambda)$value)
}

# Z conditional expectation
exCZ = function(z, c = 0.5, shape = 1, rate = 1, lambda = 0.1) {
    sapply(z, function(z) z*truncdist::dtrunc(c, spec = 'exp', b = 10, rate = lambda*z) * dgamma(z, shape = shape, rate = rate) )
}
exZ = function(c, shape=1, rate=1, lambda=0.1){
    sapply(c, function(c) integrate(exCZ, lower = 0, upper = Inf, c = c, shape = shape, rate = rate, lambda = lambda)$value)
}
exZpC = function(t, shape=1, rate=1, lambda=0.1){
    sapply(t, function(t) integrate(exZ, lower = 10*t, upper = 10, shape = shape, rate = rate, lambda = lambda)$value)
}
exZcondpC = function(x, shape=1, rate=1, lambda=0.1){
    sapply(x, function(x) exZpC(x, shape = shape, rate = rate, lambda = lambda)/pC(x, shape = shape, rate = rate, lambda = lambda))
}

# Z^2 conditional expectation
exCZ2 = function(z, c = 0.5, shape = 1, rate = 1, lambda = 0.1) {
    sapply(z, function(z) z^2*truncdist::dtrunc(c, spec = 'exp', b = 10, rate = lambda*z) * dgamma(z, shape = shape, rate = rate) )
}
exZ2 = function(c, shape=1, rate=1, lambda=0.1){
    sapply(c, function(c) integrate(exCZ2, lower = 0, upper = Inf, c = c, shape = shape, rate = rate, lambda = lambda)$value)
}
exZ2pC = function(t, shape=1, rate=1, lambda=0.1){
    sapply(t, function(t) integrate(exZ2, lower = 10*t, upper = 10, shape = shape, rate = rate, lambda = lambda)$value)
}
exZ2condpC = function(x, shape=1, rate=1, lambda=0.1){
    sapply(x, function(x) exZ2pC(x, shape = shape, rate = rate, lambda = lambda)/pC(x, shape = shape, rate = rate, lambda = lambda))
}

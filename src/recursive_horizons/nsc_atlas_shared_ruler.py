"""Full-covariance shared BAO ruler with one-dimensional SN posterior propagation."""
from __future__ import annotations

import math
import time

import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar,brentq
from scipy.special import logsumexp,ndtr,log_ndtr

from .nsc_atlas_parent_supernova import flrw_distance,parent_fractions,SupernovaShapeLikelihood

C_OVER_100=2997.92458  # c[km/s]/(100 km/s/Mpc); r_d h is then in Mpc.


def fractions(model,theta):
    if model=='parent':
        d=parent_fractions(math.exp(-float(theta)))
        return dict(Omega_m=d['matter'],Omega_Lambda=d['vacuum'],Omega_k=d['curvature'])
    if model=='flat_LCDM' and 0<=theta<=1:
        return dict(Omega_m=float(theta),Omega_Lambda=1-float(theta),Omega_k=0.)
    raise ValueError('parent or flat_LCDM in its admitted shape interval required')


def distance(model,theta,z,hel):
    f=fractions(model,theta)
    return flrw_distance(z,hel,matter=f['Omega_m'],curvature=f['Omega_k'])


def bao_shape(model,theta,redshifts,labels):
    z=np.asarray(redshifts,float)
    f=fractions(model,theta)
    DM=distance(model,theta,z,z)/(1+z)
    DH=1/np.sqrt(f['Omega_m']*(1+z)**3+f['Omega_k']*(1+z)**2+f['Omega_Lambda'])
    DV=np.cbrt(z*DM*DM*DH)
    values={'DM_over_rs':DM,'DH_over_rs':DH,'DV_over_rs':DV}
    if any(k not in values for k in labels):raise ValueError('unrecognized BAO observable')
    return np.array([values[k][i] for i,k in enumerate(labels)])


class SharedRuler:
    def __init__(self,mean,covariance):
        self.mean=np.asarray(mean,float)
        C=np.asarray(covariance,float)
        if (C.shape!=(len(self.mean),len(self.mean)) or not np.isfinite(C).all()
                or not np.isfinite(self.mean).all() or np.max(abs(C-C.T))>1e-12):
            raise ValueError('finite aligned BAO mean and symmetric covariance required')
        self.factor=cho_factor(C,lower=True)
        self.precision_mean=cho_solve(self.factor,self.mean)
        self.c=float(self.mean@self.precision_mean)

    def coefficients(self,shape):
        t=np.asarray(shape,float)
        if t.shape!=self.mean.shape or np.any(t<=0):raise ValueError('positive aligned BAO shape required')
        a=float(t@cho_solve(self.factor,t));b=float(t@self.precision_mean)
        if a<=0 or b<=0:raise ValueError('positive physical scale not resolved')
        return a,b,self.c

    def profile(self,shape):
        a,b,c=self.coefficients(shape)
        beta=b/a
        # Evaluate the residual directly to avoid cancellation in c-b^2/a.
        residual=self.mean-beta*np.asarray(shape,float)
        chi=float(residual@cho_solve(self.factor,residual))
        return dict(beta=beta,rd_h_Mpc=C_OVER_100/beta,chi_squared=chi,a=a,b=b,c=c)


def marginal_score(beta,log_weights,a,b,c):
    return float(-2*logsumexp(log_weights-.5*(a*beta*beta-2*b*beta+c)))


def posterior_calculation(model,z_hd,z_hel,magnitude,covariance,bao_z,bao_labels,
                          bao_mean,bao_covariance,*,order=48):
    started=time.process_time()
    if model not in ('parent','flat_LCDM') or order not in (48,96):raise ValueError('declared model/order required')
    z=np.asarray(z_hd,float);hel=np.asarray(z_hel,float);bz=np.asarray(bao_z,float)
    sn=SupernovaShapeLikelihood(magnitude,covariance);bao=SharedRuler(bao_mean,bao_covariance)
    lo_sn=max(0.,math.log((1+max(z))/3)) if model=='parent' else 0.
    lo_bao=max(lo_sn,math.log((1+max(bz))/3)) if model=='parent' else 0.
    hi=math.log(20) if model=='parent' else 1.
    cache={};calls=0
    def sn_chi(theta):
        nonlocal calls
        key=float(theta)
        if key not in cache:
            calls+=1
            if calls>3500 or time.process_time()-started>180:raise RuntimeError('posterior worker compute allowance exhausted')
            cache[key]=sn.evaluate(distance(model,key,z,hel))['chi_squared']
        return cache[key]
    fit=minimize_scalar(sn_chi,bounds=(lo_sn,hi),method='bounded',options=dict(xatol=2e-11))
    if not fit.success:raise ValueError('SN shape maximum unresolved')
    theta0=float(fit.x);chi0=float(fit.fun)
    step=1e-4
    curvature=(sn_chi(theta0+step)-2*chi0+sn_chi(theta0-step))/step**2
    if curvature<=0:raise ValueError('positive likelihood curvature required for quadrature placement')
    sigma=math.sqrt(2/curvature)
    # Integrate the entire stated prior; split only to resolve narrow peaks.
    breaks=sorted(set([lo_sn,hi,lo_bao,*[max(lo_sn,min(hi,theta0+k*sigma)) for k in (-8,-4,0,4,8)]]))
    gx,gw=np.polynomial.legendre.leggauss(order)
    nodes=[];qweights=[]
    for left,right in zip(breaks[:-1],breaks[1:]):
        if right>left:
            nodes.extend((left+right)/2+(right-left)*gx/2)
            qweights.extend((right-left)*gw/2)
    nodes=np.asarray(nodes);qweights=np.asarray(qweights)
    sn_values=np.array([sn_chi(x) for x in nodes])
    logw=np.log(qweights)-.5*(sn_values-chi0)
    lognorm=float(logsumexp(logw));w_all=np.exp(logw-lognorm)
    allowed=nodes>=lo_bao
    retained=float(w_all[allowed].sum());excluded=float(w_all[~allowed].sum())
    nodes_joint=nodes[allowed];sn_joint=sn_values[allowed]
    logw_sn=logw[allowed]-logsumexp(logw[allowed])
    shapes=np.array([bao_shape(model,x,bz,bao_labels) for x in nodes_joint])
    coefficients=np.array([bao.coefficients(t) for t in shapes])
    a,b,c=coefficients.T
    profiles=np.array([bao.profile(t)['chi_squared'] for t in shapes])
    def joint_objective(theta):
        return sn_chi(theta)+bao.profile(bao_shape(model,theta,bz,bao_labels))['chi_squared']
    joint_fit=minimize_scalar(joint_objective,bounds=(lo_bao,hi),method='bounded',options=dict(xatol=2e-11))
    if not joint_fit.success:raise ValueError('joint shape optimum unresolved')
    best=float(joint_fit.x);scale=bao.profile(bao_shape(model,best,bz,bao_labels))
    max_beta=float(np.max(b/a+12/np.sqrt(a)))
    marg=minimize_scalar(lambda beta:marginal_score(beta,logw_sn,a,b,c),bounds=(0,max_beta),method='bounded',
                         options=dict(xatol=1e-10))
    if not marg.success:raise ValueError('posterior predictive ruler optimum unresolved')
    # Integrate beta>0 analytically with a common flat beta prior.
    log_joint=logw_sn-.5*profiles-.5*np.log(a)+log_ndtr(b/np.sqrt(a))
    wjoint=np.exp(log_joint-logsumexp(log_joint))
    beta_means=b/a;beta_sd=1/np.sqrt(a);lower_cdf=ndtr(-beta_means/beta_sd)
    positive_mass=ndtr(beta_means/beta_sd)
    def beta_cdf(value):
        return float(wjoint@((ndtr((value-beta_means)/beta_sd)-lower_cdf)/positive_mass))
    beta_quantiles=[brentq(lambda value:beta_cdf(value)-p,0,max_beta,xtol=2e-11) for p in (.16,.5,.84)]
    rd_quantiles=[C_OVER_100/v for v in beta_quantiles[::-1]]
    def moments(values,weights):
        values=np.asarray(values,float)
        mean=float(weights@values)
        return dict(mean=mean,standard_deviation=math.sqrt(max(0.,float(weights@(values-mean)**2))))
    f_all=[fractions(model,x) for x in nodes]
    f_joint=[fractions(model,x) for x in nodes_joint]
    sn_moments={key:moments([f[key] for f in f_all],w_all) for key in ('Omega_m','Omega_Lambda','Omega_k')}
    joint_moments={key:moments([f[key] for f in f_joint],wjoint) for key in ('Omega_m','Omega_Lambda','Omega_k')}
    f_best=fractions(model,best)
    omegaL=joint_moments['Omega_Lambda']
    score_comparison=dict(predicted=3*f_best['Omega_Lambda'],comparison=2.054,comparison_sigma=.022,
        posterior_mean=3*omegaL['mean'],posterior_sigma=3*omegaL['standard_deviation'],
        standardized_difference_including_posterior=(3*omegaL['mean']-2.054)/math.hypot(3*omegaL['standard_deviation'],.022))
    return dict(model=model,order=order,
        SN_calibration=dict(best_theta=theta0,chi_squared=chi0,shape_parameters=1,magnitude_intercepts=1,
                            prior=('uniform ln(Ao) on [0,ln(20)], conditioned on proved SN emitter domain' if model=='parent' else
                                   'uniform Omega_m on [0,1]'),moments=sn_moments),
        posterior_propagation=dict(nodes=len(nodes),entire_prior_integrated=True,prior_bounds=[lo_sn,hi],
            BAO_admissible_lower=lo_bao,SN_posterior_mass_outside_BAO_domain=excluded,
            predictive_log_integral_normalizer=lognorm,
            marginal_BAO_deviance=float(marg.fun),marginal_best_rd_h_Mpc=C_OVER_100/float(marg.x),
            joint_prior='SN shape prior times a common flat beta prior on beta>0; no model Bayes factor',
            joint_rd_h_quantiles_16_50_84_Mpc=rd_quantiles,joint_fraction_moments=joint_moments),
        joint_best=dict(theta=best,**f_best,rd_h_Mpc=scale['rd_h_Mpc'],
            chi_squared_SN=sn_chi(best),chi_squared_BAO=scale['chi_squared'],
            chi_squared_total=sn_chi(best)+scale['chi_squared'],fitted_parameters=3,
            parameter_meaning=['one shape','one SN intercept','one common BAO ruler'],
            observer_scale=math.exp(best) if model=='parent' else None),
        three_Omega_Lambda_comparison=score_comparison,
        numerical=dict(SN_likelihood_evaluations=calls,process_CPU_seconds=time.process_time()-started,
                       new_field_evolutions=0))

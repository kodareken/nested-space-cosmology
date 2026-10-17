"""One fixed-cap/material-wall stage with a common exterior scalar value.

The exterior callback supplies first jets on that same value. Its bulk
constraints belong to the caller. This adapter takes no time step and
retains strong/weak scalar and temporal Einstein residuals explicitly.
"""
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np

from .nsc_atlas_cap_fixed_reference import fixed_cap_reference
from .nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from .nsc_atlas_polar_cap import _scalar
from .nsc_atlas_background_cap import _budget


def fixed_cap_hybrid_reference(basis,chi,canonical_density,Z,*,alpha,exterior_wall,
                               C_wall,eta_plus,eta_plus_x,budget,rtol=2e-13,atol=2e-15,
                               background=None,background_time=None,radial_step_fraction=None,
                               background_panel_sides=False,compiled_coefficients=False):
    """Derive the three cap boundary inputs from one exterior/hybrid event.

    exterior_wall receives immutable {r,c} and returns r_V,c_V,log_C_V.
    C_wall is an evolved metric value; eta_plus and its x derivative define
    the existing fixed-wall-coordinate gauge. No physical lapse rate or
    scalar waveform can be prescribed through this interface.
    """
    budget=_budget(budget);before=budget.evaluations
    if not callable(exterior_wall): raise ValueError('exterior wall must be callable')
    C=_scalar(C_wall,'wall C',positive=True)
    boost=_scalar(eta_plus,'wall boost');boost_x=_scalar(eta_plus_x,'wall boost derivative')
    connection={}

    def couple(packet):
        R,f=packet['radius'],packet['f']
        q=np.sqrt(packet['alpha']/6)*f
        u=f*(np.sinh(q)/q if q else 1.)
        c=float(R*u)
        if not np.isfinite(c): raise ValueError('nonfinite common scalar value')
        common=MappingProxyType(dict(r=R,c=c))
        raw=exterior_wall(common)
        names={'r_V','c_V','log_C_V'}
        if not isinstance(raw,Mapping) or set(raw)!=names:
            raise ValueError('exterior wall must return the three outgoing first jets')
        outgoing={key:_scalar(raw[key],'exterior '+key) for key in names}
        state=dict(common,eta_plus=boost)
        cap=dict(A=packet['A'],delta=packet['delta'],w_minus=(packet['p']-packet['h'])/R)
        first=hybrid_wall_boundary(state,outgoing,cap,dict(C_wall=C,log_C_s=0.),alpha=packet['alpha'])
        xs=float(first['Vdot_plus'])
        Lx=2*boost_x-2*float(first['k_plus'])/xs+2*outgoing['log_C_V']
        hybrid=hybrid_wall_boundary(state,outgoing,cap,dict(C_wall=C,log_C_s=Lx*xs),alpha=packet['alpha'])
        rates=dict(r=float(hybrid['r_s'])/xs,L=Lx,T=float(hybrid['T_s'])/xs,
                   s=1/xs,c=float(hybrid['c_s'])/xs)
        connection.update(common=dict(common),outgoing=outgoing,hybrid=hybrid,
                          wall_rates=rates,rho=packet['rho'])
        budget.check(512)
        return dict(T_x=rates['T'],radius_rate=rates['r'],incoming_w_plus=float(hybrid['w_plus']))

    cap=fixed_cap_reference(basis,chi,canonical_density,Z,alpha=alpha,wall_coupler=couple,
                            budget=budget,rtol=rtol,atol=atol,
                            background=background,background_time=background_time,
                            radial_step_fraction=radial_step_fraction,
                            background_panel_sides=background_panel_sides,
                            compiled_coefficients=compiled_coefficients)
    hybrid=connection['hybrid'];rates=connection['wall_rates'];wall=cap['wall']
    R=cap['radius'];tx,rx=rates['T'],rates['r'];mu=float(np.exp(cap['lapse']))
    cap_f_x=(tx*wall['velocity']+rx*wall['h_trace'])/R
    target_f_x=float(hybrid['f_s'])*rates['s']
    Ef=cap_f_x-target_f_x
    dp,dh=wall['momentum_discrepancy'],wall['normal_discrepancy']
    factored=((rx-tx*wall['B'])*dp+(rx+tx*wall['B'])*dh)/(2*R)
    rootD=float(np.sqrt(hybrid['D']))
    cap_c_x=(connection['common']['c']/R)*rx+R*rootD*cap_f_x
    temporal=(cap['mass_rate']-cap['moving_density'])/tx
    factor=wall['A']**2*mu/2
    actual=factor*wall['p_trace']*wall['h_trace']
    matched=factor*R**2*float(hybrid['Pi'])*float(hybrid['Phi'])
    velocity_channel=tx*wall['velocity_discrepancy']/R
    incoming_channel=(rx+tx*wall['B'])*wall['incoming_discrepancy']/2
    residual=dict(scalar_tangent=Ef,scalar_tangent_factored=float(factored),
        scalar_velocity_channel=float(velocity_channel),
        scalar_incoming_channel=float(incoming_channel),
        scalar_channel_identity_defect=float(Ef-velocity_channel-incoming_channel),
        canonical_wall_flux=R*Ef,conformal_wall_flux=cap_c_x-rates['c'],
        conformal_wall_flux_factored=R*rootD*factored,
        temporal_mass_actual=temporal-actual,temporal_mass_matched=temporal-matched,
        weak_flux_actual=cap['physical_weak_flux']-actual,
        weak_flux_matched=cap['physical_weak_flux']-matched,
        physical_temporal_mass_actual=cap['alpha']*(temporal-actual),
        physical_temporal_mass_matched=cap['alpha']*(temporal-matched),
        strong_minus_weak_outgoing=(dp-dh)/R)
    physical_mass_rate=1.5*R*R*rx+cap['alpha']*cap['mass_rate']
    physical_density=1.5*R*R+cap['alpha']*connection['rho']
    if not all(np.isfinite(v) for v in (*rates.values(),*residual.values(),
                                       physical_mass_rate,physical_density)):
        raise ValueError('nonfinite fixed-cap hybrid rate or residual')
    budget.check(144)
    return dict(cap=cap,common=connection['common'],outgoing=connection['outgoing'],
        hybrid=hybrid,wall_rates=rates,residual=residual,cap_f_x=cap_f_x,cap_c_x=cap_c_x,
        physical_mass_rate=physical_mass_rate,physical_mass_density=physical_density,
        charged_evaluations=budget.evaluations-before,
        scope='one shared-value finite-dimensional stage; no bulk exterior constraint or evolution certificate')

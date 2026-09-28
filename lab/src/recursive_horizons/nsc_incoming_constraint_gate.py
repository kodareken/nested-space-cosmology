"""Retained-state incoming momentum bound for the compact response family.

Apply the owned Killing-power normalization and current-normalized scattering
to the fixed NSC occupation law. No stress quadrature, mode solve or Einstein
equation is derived here. The certificate concerns the frozen incoming data
of SuppliedKSHarmonicMetric, not every transmitting spherical geometry.
"""
from fractions import Fraction
from math import exp, factorial, pi


def retained_incoming_flux_bound(ledger, channels, preparation):
    if sorted(c['index'] for c in channels) != list(range(33)):
        raise ValueError('the authenticated retained 33-group inventory is required')
    q, omega, kappa = ledger['magnetic_flux'], ledger['Omega'], preparation['surface_gravity']
    if q != preparation['magnetic_flux'] or omega != preparation['omega']:
        raise ValueError('source and locked ledger disagree')
    if q != 4 or not 3.97 < omega < 3.98 or not .238 < kappa < .239 or ledger['A'] <= 0:
        raise ValueError('parameters outside the certified retained branch')
    lll = channels[0]
    if lll['compact_mass'] != 0 or lll['angular_eigenvalue'] != 0 or lll['degeneracy'] != abs(q) or lll['copy_count'] != 1:
        raise ValueError('LLL multiplicity or mass changed')
    beta = 2*pi/kappa
    power_lll = abs(q)*kappa*kappa*(1-omega*omega)/(48*pi)
    massive = {}; nonpositive_massless = []
    for c in channels[1:]:
        m = c['compact_mass']; count = int(c['copy_count'])*int(c['degeneracy'])
        if m == 0:
            nonpositive_massless.append(c['index']); continue
        if m < 1.5 or count <= 0:
            raise ValueError('mass-gap or multiplicity outside the declared bound')
        key = str(c['compact_level'])
        row = massive.setdefault(key, {'mass': m, 'multiplicity': 0, 'channel_indices': []})
        if row['mass'] != m: raise ValueError('compact mass not common within level')
        row['multiplicity'] += count; row['channel_indices'].append(c['index'])
    total = sum(row['multiplicity'] for row in massive.values())
    if total != 1024: raise ValueError('retained massive multiplicity changed')
    for row in massive.values():
        m = row['mass']
        row['positive_power_upper_bound'] = row['multiplicity']/pi*exp(-beta*m)*(m/beta+1/beta**2)
    power_upper = power_lll+sum(row['positive_power_upper_bound'] for row in massive.values())
    a2, r2 = 3*pi/2-4, 2.
    # A separate rational witness proves the sign without a near-zero float
    # decision. These are outward arithmetic bounds, not parameter changes.
    F = Fraction
    pi_lower, pi_upper = F(314159,100000), F(314160,100000)
    k_lower, k_upper, omega_lower = F(238,1000), F(239,1000), F(397,100)
    beta_lower, mass_lower = F(26), F(3,2)
    e_lower = F(27,10)
    if not (pi_lower < pi < pi_upper and 2*pi_lower/k_upper > beta_lower
            and sum(F(1,factorial(i)) for i in range(6)) > e_lower
            and e_lower**39 > 10**16):
        raise ArithmeticError('elementary outward-bound witness failed')
    lll_magnitude_lower = abs(q)*k_lower**2*(omega_lower**2-1)/(48*pi_upper)
    # exp(-beta*m) < exp(-39) < 10^-16, and the polynomial factor
    # is bounded by the integral from m_lower at beta_lower.
    positive_upper = F(total)/pi_lower*F(1,10**16)*(mass_lower/beta_lower+1/beta_lower**2)
    margin = lll_magnitude_lower-positive_upper
    if margin <= 0: raise ArithmeticError('incoming flux sign not decided')
    t01_lower = margin/(8*pi_upper*(3*pi_upper/2-4))
    normalized_miss = t01_lower/(2*F(str(ledger['A'])))
    return {
        'beta_H': beta, 'LLL_Killing_power': power_lll, 'massive_emission_bounds': massive,
        'nonpositive_massless_angular_channels': nonpositive_massless,
        'Killing_power_upper_bound': power_upper,
        'incoming_geometry': {'rho': 1., 'r_squared': r2, 'a_squared': a2},
        'T01_lower_bound_from_sharp_estimates': -power_upper/(4*pi*r2*a2),
        'rational_witness': {
            'pi_interval': [str(pi_lower),str(pi_upper)], 'kappa_interval': [str(k_lower),str(k_upper)],
            'Omega_lower': str(omega_lower), 'mass_lower': str(mass_lower),
            'beta_lower': str(beta_lower), 'exp_minus39_upper': '1/10000000000000000',
            'LLL_magnitude_lower': str(lll_magnitude_lower), 'massive_positive_upper': str(positive_upper),
            'strict_negative_power_margin': str(margin), 'T01_strict_lower': str(t01_lower),
            'normalized_shift_miss_lower': str(normalized_miss),
            'strict_sign': True, 'T01_lower_decimal': float(t01_lower),
            'normalized_shift_miss_lower_decimal': float(normalized_miss)},
    }

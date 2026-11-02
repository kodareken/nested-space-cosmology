"""Generate the 4-panel Non-Local Multiverse to Local Observer Atlas figure and summary JSON."""
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from recursive_horizons.nsc_atlas_offset_generator import (
    nonlocal_to_local_atlas_readout,
    unified_least_resistance_cycle,
    nested_zero_crossing_map,
    buckling_to_goldilocks_calculation,
)

ROOT = Path(__file__).resolve().parents[1]


def main():
    atlas = nonlocal_to_local_atlas_readout()
    pipeline = buckling_to_goldilocks_calculation()
    summary_path = ROOT / 'results/summaries/nsc-nonlocal-to-local-atlas-v1.json'
    record = {
        'atlas': atlas,
        'pipeline': pipeline,
    }
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    plt.rcParams.update({
        'font.size': 10,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)

    # Panel A: Circle-Pendulum Duality & Smooth 0-Crossing vs Unsigned |R| Cusp
    ax = axes[0, 0]
    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    r_signed = np.cos(theta)
    r_unsigned = np.abs(r_signed)
    det_f = r_signed ** 3
    ax.plot(theta / np.pi, r_signed, color='#2357A6', lw=2.2, label=r'Signed radius $R(\tau)=R_0\cos(\omega\tau)$ (smooth)')
    ax.plot(theta / np.pi, r_unsigned, color='#C2410C', lw=1.8, ls='--', label=r'Unsigned size $r(\tau)=|R(\tau)|$ (cusp at $0$)')
    ax.plot(theta / np.pi, det_f, color='#21766C', lw=1.5, ls=':', label=r'3D orientation $\det F=(R/R_0)^3$')
    ax.axhline(0.0, color='#526779', lw=0.8)
    ax.scatter([0.5, 1.5], [0.0, 0.0], color='#C2410C', s=42, zorder=5)
    ax.annotate('Smooth $0$ crossing\n($\\Delta\\dot{R}=0$, $\\Delta\\dot{r}=+2R_0\\omega$ cusp)',
                xy=(0.5, 0.0), xytext=(0.04, -0.82),
                arrowprops=dict(arrowstyle='->', color='#122033', lw=0.9),
                fontsize=8.2, color='#122033')
    ax.set_title('A  Circle-Pendulum Duality & Smooth $0$-Crossing\n$(-)\\times(-)=(+)$ across two half-turns', fontsize=10.5, fontweight='bold')
    ax.set_xlabel(r'Offset phase $\omega\tau/\pi$ (turns $\times 2$)')
    ax.set_ylabel('Normalized coordinate / Jacobian')
    ax.set_ylim(-1.22, 1.48)
    ax.legend(fontsize=7.8, loc='upper right', frameon=True, facecolor='white', framealpha=0.9, edgecolor='none')
    ax.grid(alpha=0.22)

    # Panel B: Unified 5-Variable Least-Resistance Cycle
    ax = axes[0, 1]
    chis = np.linspace(0.0, 3.5, 250)
    rows = [unified_least_resistance_cycle(float(c)) for c in chis]
    v_plus = [r['plus_infall']['parcel_volume'] for r in rows]
    f_drop = [r['derived_cooker_takeover']['minus_volume_fraction'] for r in rows]
    aspect = [min(r['plus_infall']['filament_aspect_ratio'], 12.0) for r in rows]
    chi_spag = 0.5 * math.log(2.0)
    ax.plot(chis, v_plus, color='#2357A6', lw=2.0, label=r'Infalling $+$ parcel volume $V_+^*(\chi)\to 0$')
    ax.plot(chis, f_drop, color='#21766C', lw=2.2, label=r'Expanded $-$ sphere share $V_-^*/(V_+^*+V_-^*)\to 1$')
    ax.plot(chis, np.array(aspect) / 10.0, color='#C27517', lw=1.8, ls='-.', label=r'$+$ filament aspect ratio $\times 0.1$')
    ax.axvline(chi_spag, color='#8C4264', lw=1.0, ls=':')
    ax.annotate(f'Spaghettification onset\n$\\chi={chi_spag:.4f}$ ($\\theta^*=1/2$)',
                xy=(chi_spag, 0.5), xytext=(0.45, 0.86),
                arrowprops=dict(arrowstyle='->', color='#8C4264', lw=0.9),
                fontsize=8.2, color='#8C4264',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.88))
    ax.set_title('B  Single 5-Variable Cost Minimum ($E=Q_{\\mathrm{stretch}}+PV$)\nParcel shrinkage, 1D filament & 3D sphere without fluid rules', fontsize=10.5, fontweight='bold')
    ax.set_xlabel(r'Inherited log-length compression $\chi$')
    ax.set_ylabel('Derived geometric fraction / scaled ratio')
    ax.set_ylim(0.0, 1.08)
    ax.legend(fontsize=7.8, loc='center right', frameon=True, facecolor='white', framealpha=0.9, edgecolor='none')
    ax.grid(alpha=0.22)

    # Panel C: Two-Valve Steady-State Fixed Point chi_* & 95.02% Dark Sector
    ax = axes[1, 0]
    chis_c = np.linspace(0.2, 4.2, 240)
    c_rows = [nested_zero_crossing_map(1, float(c))['nested_chain'][0] for c in chis_c]
    om_de = [100.0 * r['dark_energy_like_fraction'] for r in c_rows]
    om_dm = [100.0 * r['dark_matter_like_fraction'] for r in c_rows]
    om_vis = [100.0 * r['visible_filament_like_fraction'] for r in c_rows]
    om_dark = [100.0 * r['dark_sector_combined_fraction'] for r in c_rows]
    fp = atlas['two_valve_fixed_point']
    chi_star = fp['chi_star']
    ax.plot(chis_c, om_dark, color='#0A1730', lw=2.0, ls='--', label='Combined Dark Sector ($95.02\\%$ at $\\chi_*$)')
    ax.plot(chis_c, om_de, color='#2357A6', lw=2.2, label='Dark Energy share $\\Omega_{\\mathrm{DE}}$ ($71.14\\%$)')
    ax.plot(chis_c, om_dm, color='#21766C', lw=2.0, label='Dark Matter share $\\Omega_{\\mathrm{DM}}$ ($23.88\\%$)')
    ax.plot(chis_c, om_vis, color='#C27517', lw=2.0, label='Visible Filaments $\\Omega_{\\mathrm{vis}}$ ($e^{-3}=4.98\\%$)')
    ax.axvline(chi_star, color='#8C4264', lw=1.2, ls=':')
    ax.scatter([chi_star]*3, [100*fp['dark_energy_fraction'], 100*fp['dark_matter_fraction'], 100*fp['visible_filament_fraction']],
               color='#8C4264', s=32, zorder=5)
    ax.annotate(f'Two-valve fixed point\n$\\chi_*={chi_star:.4f}$\n($71.14\\% / 23.88\\% / 4.98\\%$)',
                xy=(chi_star, 71.14), xytext=(2.05, 42.0),
                arrowprops=dict(arrowstyle='->', color='#8C4264', lw=0.9),
                fontsize=8.2, color='#8C4264',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.88))
    ax.set_title('C  Two-Valve Steady State ($\\Omega_{\\mathrm{vis}}(\\chi_*)=e^{-3}$)\n3:1 Virial lock & $1-e^{-3}=95.02\\%$ exponential horizon', fontsize=10.5, fontweight='bold')
    ax.set_xlabel(r'Inherited compression $\chi$')
    ax.set_ylabel('Local observer budget share (%)')
    ax.set_ylim(0.0, 108.0)
    ax.legend(fontsize=7.6, loc='upper left', frameon=True, facecolor='white', framealpha=0.9, edgecolor='none')
    ax.grid(alpha=0.22)

    # Panel D: Discrete 3D Volumetric Harmonic Asymmetry Ladder eta_K
    ax = axes[1, 1]
    ks = np.arange(1, 9)
    etas = [math.tanh(0.25 * math.pi * math.exp(-3.0 * k)) for k in ks]
    ax.semilogy(ks, etas, 'o-', color='#2357A6', lw=2.0, markersize=6,
                label=r'$\eta_K = \tanh\!\left(\frac{\pi}{4}e^{-3K}\right)$')
    ax.axhspan(5.8e-10, 6.5e-10, color='#21766C', alpha=0.25, label='Observed BBN band $[5.8, 6.5]\\times 10^{-10}$')
    ax.axhline(6.1e-10, color='#21766C', lw=1.1, ls='--')
    eta_7 = etas[6]
    ax.scatter([7], [eta_7], color='#C2410C', s=55, zorder=6)
    ax.annotate(f'$K=7$ ($21$ e-folds)\n$\\eta_7 = {eta_7:.3e}$',
                xy=(7, eta_7), xytext=(4.0, 3.0e-9),
                arrowprops=dict(arrowstyle='->', color='#C2410C', lw=1.0),
                fontsize=8.5, color='#C2410C', fontweight='bold')
    ax.set_title('D  Offset Time-Decay Asymmetry Ladder $\\eta_K$\n$K=7$ harmonic inside the $[\\lambda_{\\min},\\lambda_{\\max}]$ band', fontsize=10.5, fontweight='bold')
    ax.set_xlabel(r'3D volumetric e-fold harmonic stage $K$ ($3K$ e-folds)')
    ax.set_ylabel(r'Matter-antimatter residual $\eta_K$')
    ax.set_xticks(ks)
    ax.legend(fontsize=8, loc='upper right', frameon=False)
    ax.grid(alpha=0.22, which='both')

    fig.suptitle('Nested-Space Cosmology: Non-Local Multiverse Geometry to Local Observer Readout (Axioms A0–A4)',
                 fontsize=12.5, fontweight='bold')
    out_fig = ROOT / 'paper/figures/nsc-nonlocal-to-local-atlas.png'
    fig.savefig(out_fig, dpi=180, metadata={'Software': 'NSC non-local to local atlas'})
    plt.close(fig)
    print(f'Wrote {summary_path.relative_to(ROOT)} and {out_fig.relative_to(ROOT)}')


if __name__ == '__main__':
    main()

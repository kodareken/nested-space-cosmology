"""Plot compact supernova and AP comparisons; no likelihood or field replay."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def main():
    read=lambda name:json.loads((ROOT/'results/summaries'/name).read_text())
    old=read('nsc-atlas-parent-supernova-v1.json')['comparisons']['all_noncalibrators_z_ge_0_025']['models']
    new=read('nsc-atlas-parent-supernova-v2.json')['comparisons']['all_noncalibrators_z_ge_0_025']['models']
    bao=read('nsc-atlas-parent-bao-validation-v1.json')
    fig,axs=plt.subplots(1,2,figsize=(11.8,4.35),layout='constrained')
    fig.set_facecolor('#faf9f6')
    base=new['flat_LCDM_calibrated']['chi_squared']
    vals=[new['quiet_de_sitter']['chi_squared']-base,
          old['parent_calibrated_epoch']['chi_squared']-base,
          new['parent_fixed_Ao_2']['chi_squared']-base,
          new['parent_calibrated_epoch']['chi_squared']-base,0]
    labels=['Quiet cap','Parent: future only','Parent: fixed $A_o=2$','Parent: earlier continuation','Flat $\\Lambda$CDM: calibrated']
    colors=['#9c4d43','#c2853d','#c8aa7f','#267e82','#375770']
    y=np.arange(len(vals))
    axs[0].barh(y,vals,color=colors,height=.58)
    axs[0].set_yticks(y,labels);axs[0].invert_yaxis()
    axs[0].set_xlim(0,710)
    for yi,v in enumerate(vals):axs[0].text(v+9,yi,f'{v:.2f}',va='center',fontsize=9)
    axs[0].set_xlabel(r'$\chi^2-\chi^2_{\mathrm{flat,cal}}$')
    axs[0].set_title('1322 light curves: same covariance',fontsize=11)
    axs[0].text(.02,-.24,'Parent source and curvature fixed; only observer epoch calibrated.\nAll histories retain a common fitted magnitude intercept.',
                transform=axs[0].transAxes,fontsize=8.5)
    z=np.array(bao['redshifts'])
    for model,offset,color,label in [('parent',-.018,'#267e82','Inherited parent'),('flat_LCDM',.018,'#375770',r'Flat $\Lambda$CDM')]:
        row=bao['models'][model]
        axs[1].plot(z+offset,row['standardized_marginal_residuals'],'o-',color=color,
                    label=label+f" ($\\chi^2={row['chi_squared']:.2f}$)",lw=1.3,ms=5)
    axs[1].axhline(0,color='#666666',lw=.7)
    axs[1].axhspan(-1,1,color='#888888',alpha=.08)
    axs[1].set_xlabel('DESI effective redshift');axs[1].set_ylabel('Linear AP residual / marginal uncertainty')
    axs[1].set_title('SN-calibrated shapes; zero BAO parameters fitted',fontsize=10.5)
    axs[1].set_ylim(-2.6,2.4);axs[1].legend(frameon=False,fontsize=8.5)
    axs[1].text(.01,-.24,'Six contrasts; joint statistic retains correlations.\nConditional on SN point estimates; acoustic origin not derived.',
                transform=axs[1].transAxes,fontsize=8.5)
    for ax in axs:
        ax.set_facecolor('#faf9f6');ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='x',alpha=.12)
    out=ROOT/'docs/figures/nsc-atlas-parent-supernova.png'
    fig.savefig(out,dpi=180,facecolor=fig.get_facecolor(),bbox_inches='tight')
    print(out)


if __name__=='__main__':main()

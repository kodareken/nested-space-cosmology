"""Conditional CMB geometric diagnostic without adding a baryon/radiation model."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_cmb_geometry import diagnostic,PRIOR_MEAN,PRIOR_COVARIANCE

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('new output required')
    source=ROOT/'results/summaries/nsc-atlas-shared-ruler-v1.json'
    previous=json.loads(source.read_text(encoding='utf-8'))
    zstar=1089.92
    rows={name:diagnostic(name,row['joint_best']['theta'],row['joint_best']['rd_h_Mpc'],zstar)
          for name,row in previous['models'].items()}
    result=dict(schema='nsc-atlas-cmb-geometric-diagnostic-v1',
        shared_ruler=dict(path=source.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(source.read_bytes()).hexdigest()),
        models=rows,supplemental_probe_redshift=zstar,
        Planck_compression=dict(source='https://arxiv.org/html/1808.05724v1',table='Table 1, LambdaCDM+Omega_k',
            mean=PRIOR_MEAN.tolist(),marginal_covariance=PRIOR_COVARIANCE.tolist(),
            rounded_table_reconstruction=True,common_compression_for_both_models=True),
        exact_CMB_likelihood_status='undefined for the fixed parent and shape-only control: no predicted recombination surface, baryon/photon state or sound horizon',
        conditional_scope='Geometry-only probe at supplied zstar, dust/scalar continuation, unknown acoustic ratio profiled. Optimistic conditional lower bound, NOT a full CMB likelihood or a fitted acoustic mechanism; not added to the primary SN+BAO score.',
        additional_fitted_diagnostic_parameters=1,
        new_source_couplings=0,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'src/recursive_horizons/nsc_atlas_cmb_geometry.py','scripts/check_nsc_atlas_cmb_geometry.py')})
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))


if __name__=='__main__':main()

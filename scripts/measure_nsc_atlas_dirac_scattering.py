"""One frozen curved Dirac transfer/metric-response comparison; compact output."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from numba import get_num_threads, set_num_threads
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_scattering import (
    batch_null_transfer, exterior_prepared_mode, occupation_response,
    ordered_metric_response, reference_profiles)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = (
    'src/recursive_horizons/nsc_atlas_dirac_scattering.py',
    'src/recursive_horizons/nsc_atlas_dirac_radar.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'src/recursive_horizons/nsc_boundary_state.py',
    'src/recursive_horizons/nsc_influence.py',
    'tests/test_nsc_atlas_dirac_scattering.py',
    'scripts/measure_nsc_atlas_dirac_scattering.py',
)


def digest(array):
    return hashlib.sha256(np.asarray(array,dtype='<c16').tobytes(order='C')).hexdigest()


def level(chart,cells):
    s = chart.duration
    labels,profiles = reference_profiles(s,cells)
    u,v = np.meshgrid(labels,labels,indexing='ij')
    nested = chart.fields(u,v)['potential']
    cap = chart.fields(u,v,region='cap')['potential']
    inputs = np.asarray([k*base for base in (nested,cap) for k in range(1,9)])
    histories,vertices = batch_null_transfer(inputs,s/cells,s/cells,profiles=profiles)
    mode = exterior_prepared_mode(s,cells)
    summaries = []
    for k in range(1,9):
        rows,matrices = {},{}
        for name,index in (('nested',k-1),('smooth_cap',k+7)):
            h,j = histories[index],vertices[index]
            response = ordered_metric_response(j,profiles)
            matrices[name] = response
            identity_error = float(np.max(np.abs(h.conj().T@h-np.eye(2*cells))))
            hermitian_error = float(np.max(np.abs(j-j.conj().transpose(0,2,1))))
            if max(identity_error,hermitian_error) > 2e-11:
                raise RuntimeError('flux or Hermiticity failed on the actual radar chart')
            rows[name] = dict(
                flux_unitarity_error=identity_error,
                metric_vertex_Hermiticity_error=hermitian_error,
                response_operator_norm=float(np.max(np.abs(np.linalg.eigvalsh(response)))),
                response_trace=float(np.trace(response).real),
                prepared_excess_per_unit_occupation=occupation_response(
                    response,mode,angular_k=k),
                prepared_exterior_to_cap_flux=float(np.linalg.norm((h@mode)[cells:])**2),
                transfer_sha256=digest(h),metric_vertices_sha256=digest(j),
                response_sha256=digest(response))
        difference = matrices['nested']-matrices['smooth_cap']
        rows.update(angular_k=k,degeneracy=2*k,
            response_difference_operator_norm=float(np.max(np.abs(np.linalg.eigvalsh(difference)))),
            prepared_excess_nesting_change=(
                rows['nested']['prepared_excess_per_unit_occupation']-
                rows['smooth_cap']['prepared_excess_per_unit_occupation']))
        summaries.append(rows)
    return dict(cells_per_null_direction=cells,characteristic_cell_width=s/cells,
        sectors=summaries,metric_cells_per_sector=cells*cells,
        source_cells=[int(np.count_nonzero(f)) for f in profiles],
        norm='sum abs(sqrt(du)*eta_plus)^2 + sum abs(sqrt(dv)*eta_minus)^2')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('use an unused summary path with an existing parent')
    commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings = {}
    for name in SOURCES:
        raw = (ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT) != raw:
            raise ValueError('freeze the transfer method before its application: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    workers = min(8,os.cpu_count() or 1,get_num_threads())
    set_num_threads(workers)
    started = time.process_time()
    chart = WallRadarChart()
    levels = []
    with threadpool_limits(limits=1):
        for cells in (64,128):
            levels.append(level(chart,cells))
            elapsed = time.process_time()-started
            if elapsed > 180:
                raise RuntimeError('Dirac response CPU allowance exhausted; stop and diagnose')
            print(json.dumps(dict(completed_cells=cells,CPU_seconds=elapsed)),flush=True)
    comparison = []
    for coarse,fine in zip(levels[0]['sectors'],levels[1]['sectors']):
        comparison.append(dict(angular_k=fine['angular_k'],
            prepared_excess_nesting_change=fine['prepared_excess_nesting_change'],
            absolute_64_128_spread=abs(fine['prepared_excess_nesting_change']-
                                     coarse['prepared_excess_nesting_change']),
            response_difference_operator_norm=fine['response_difference_operator_norm'],
            operator_norm_absolute_64_128_spread=abs(
                fine['response_difference_operator_norm']-
                coarse['response_difference_operator_norm'])))
    report = dict(schema='nsc-atlas-dirac-scattering-v1',
        status='finite_curved_transfer_and_ordered_metric_response_evaluated',
        source_commit=commit,source_sha256=bindings,
        geometry=dict(H_ell=1.,k_wall_ell=2.,M_over_ell=.4,
            wall_radius_over_ell=[.45,.55],wall_clock_duration_over_ell=chart.duration,
            null_domain='0<=U,V<=duration',control='smooth continuation of the same cap'),
        metric_insertions=dict(deformation='delta ln N=delta ln a=lambda*f, delta r=0',
            later_cap_center_UV_over_duration=[.80,.70],
            earlier_exterior_center_UV_over_duration=[.35,.55],
            half_width_over_duration=.04,bump='exp(1-1/(1-z^2)), abs(z)<1',
            response='R=-i[J_later,J_earlier], single-copy integrated generators / hbar',
            cross_contacts='zero for these strictly ordered disjoint supports'),
        preparation_difference=dict(
            selection='declared excess-mode modulation, conditional on an admissible inherited parent',
            incoming_leg='eta_minus on U=0, exterior',
            envelope_center_over_duration=.6,envelope_half_width_over_duration=.3,
            phase='exp(-12i*V/duration)',normalization='unit canonical characteristic flux',
            covariance='C_empty=P*C_parent*P, C_nu=C_empty+nu*|f><f|, P=I-|f><f|',
            reported='2*abs(k)*f^dagger*R*f per unit occupation, dimensionless',
            root_vacuum_selected=False,positive_k_multiplets_are_separate_preparations=True,
            angular_vacuum_sum=False),
        levels=levels,comparison=comparison,
        uncertainty='64/128 spread is a numerical diagnostic, not a proven continuum error bound',
        runtime=dict(process_CPU_seconds=time.process_time()-started,Numba_workers=workers,
            BLAS_threads=1,wall_ODE_rhs=chart.wall_rhs_calls,
            P001_new_grids=0,P001_time_steps=0,baseline_vacuum_polarizations=0),
        limits=dict(process_CPU_seconds=180,working_matrix_bytes=268435456,
            estimated_main_matrix_bytes=16*256*256*16*3,summary_bytes=24000),
        nested_Psi_over_Phi=None,
        scope='Actual curved finite Dirac propagation and nonlocal metric-vertex response; '
              'conditional finite-rank state difference. Not the inherited baseline vacuum '
              'polarization, angular UV limit, gravitational/junction solve or observable Psi/Phi.')
    raw = json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode()) > 24000:
        raise RuntimeError('compact Dirac response summary exceeds its admitted size')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=report['status'],comparison=comparison,
        CPU_seconds=report['runtime']['process_CPU_seconds'],workers=workers)),flush=True)


if __name__ == '__main__':
    main()

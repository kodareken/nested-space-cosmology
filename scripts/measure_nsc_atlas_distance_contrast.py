"""One pinned supernova distance-shape contrast; no cosmological parameter fit.

Raw publisher files are cached only under .build. --fetch permits their
download if missing. Source methods must be committed before evaluation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
from urllib.parse import quote
from urllib.request import urlopen

import numpy as np
import scipy
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_cosmological_distance import (
    HIGH_BAND, LOW_BAND, distance_contrast, quiet_cap_domain,
)

ROOT = Path(__file__).resolve().parents[1]
RELEASE = 'c447f0fea703fcd0fff57de5000947b5ca81286b'
PREFIX = 'Pantheon+_Data/4_DISTANCES_AND_COVAR/'
FILES = {
    'distances.dat': dict(path=PREFIX+'Pantheon+SH0ES.dat', bytes=579283,
                          git_blob='cce857db0c15e9ce7a0e0ce77452b6ff62af969a'),
    'stat-plus-systematic.cov': dict(path=PREFIX+'Pantheon+SH0ES_STAT+SYS.cov', bytes=33284960,
                                   git_blob='d1a1498154e7ba826df14bdbef35ebcb7f5efba1'),
}
SOURCES = ('src/recursive_horizons/nsc_atlas_cosmological_distance.py',
           'scripts/measure_nsc_atlas_distance_contrast.py',
           'tests/test_nsc_atlas_cosmological_distance.py')


def cached_file(cache, name, fetch):
    definition = FILES[name]
    path = cache/name
    url = ('https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/'
           +RELEASE+'/'+quote(definition['path'], safe='/'))
    if path.is_symlink():
        raise ValueError('input cache must contain regular files')
    if path.exists():
        raw = path.read_bytes()
    elif fetch:
        with urlopen(url, timeout=60) as stream:
            raw = stream.read(definition['bytes']+1)
    else:
        raise ValueError('missing pinned input; use --fetch: '+str(path))
    blob = hashlib.sha1(b'blob '+str(len(raw)).encode('ascii')+b'\0'+raw).hexdigest()
    if len(raw) != definition['bytes'] or blob != definition['git_blob']:
        raise ValueError('publisher Git-object or byte-count mismatch: '+name)
    if not path.exists():
        with path.open('xb') as stream:
            stream.write(raw)
    return dict(**definition, sha256=hashlib.sha256(raw).hexdigest(), url=url,
                release_commit=RELEASE)


def measure(cache, fetch, threads):
    started = time.process_time()
    def check():
        if time.process_time()-started > 120:
            raise RuntimeError('distance-contrast calculation exceeded 120 CPU seconds')
    with ThreadPoolExecutor(max_workers=2) as executor:
        bindings = dict(zip(FILES, executor.map(lambda name: cached_file(cache, name, fetch), FILES)))
    check()
    table = np.genfromtxt(cache/'distances.dat', names=True, dtype=None, encoding='utf-8')
    required = {'CID', 'zHD', 'zHEL', 'm_b_corr', 'IS_CALIBRATOR'}
    if not required.issubset(table.dtype.names or ()) or table.size != 1701:
        raise ValueError('unexpected publisher table format or row count')
    cov_path = cache/'stat-plus-systematic.cov'
    with cov_path.open('r', encoding='ascii') as stream:
        dimension = int(stream.readline())
    values = np.loadtxt(cov_path, skiprows=1)
    if dimension != table.size or values.size != dimension*dimension:
        raise ValueError('covariance and table dimensions differ')
    C = values.reshape(dimension, dimension)
    check()
    with threadpool_limits(limits=threads):
        indices, _, statistics = distance_contrast(
            table['zHD'], table['zHEL'], table['m_b_corr'], C, table['IS_CALIBRATOR'])
    check()
    return dict(
        schema='nsc-atlas-distance-contrast-v1', status='completed_observational_benchmark',
        sources=bindings,
        selection=dict(non_calibrators_only=True, low_zHD=list(LOW_BAND), low_upper_inclusive=False,
                       high_zHD=list(HIGH_BAND), high_upper_inclusive=True,
                       released_row_order_preserved=True, repeated_light_curves_retained=True,
                       selected_index_sha256=hashlib.sha256(indices.astype('<i8').tobytes()).hexdigest(),
                       distinct_selected_CID_values=int(np.unique(table['CID'][indices]).size)),
        statistics=statistics, quiet_cap_domain=quiet_cap_domain(),
        interpretation=dict(
            observable='Two-band linear contrast of m_b_corr - 5 log10(zHD*(1+zHEL)), using the joint published covariance.',
            exponential_contrast='Exp[ln(10)*Delta/5]; transformed central value and one-standard-deviation interval, not an unbiased mean ratio.',
            candidate='Existing quiet de Sitter cap with comoving standard-candle probes; source amplitude zero.',
            common_magnitude_and_distance_scale='eliminated by zero-sum contrast weights; no Cepheid distance modulus used',
            error_model='Published statistical plus systematic covariance, including publisher standardization and fiducial peculiar-velocity corrections.',
            validation_scope='Retrospective quiet-limit consistency check, not withheld validation of the full NSC cosmology or nonzero-source P001.',
            Lambda_scope='No induced Lambda is computed by fitting these observations; the supplied cap volume term remains an input.'),
        numerical_budget=dict(CPU_seconds_limit=120, BLAS_thread_limit=threads,
                              parameter_grids=0, new_field_evolutions=0),
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
    ), time.process_time()-started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=ROOT/'.build/data/pantheonplus-c447f0fea703')
    parser.add_argument('--fetch', action='store_true')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cache = args.cache if args.cache.is_absolute() else ROOT/args.cache
    cache = cache.resolve()
    output = args.output if args.output.is_absolute() else ROOT/args.output
    output = output.resolve()
    if (not cache.is_relative_to((ROOT/'.build').resolve()) or not cache.is_relative_to(ROOT.resolve())
            or not output.is_relative_to(ROOT.resolve()) or output.suffix != '.json'
            or output.exists() or not output.parent.is_dir() or not 1 <= args.threads <= 8):
        raise ValueError('use an ignored .build cache, a new repository JSON output, and 1..8 BLAS threads')
    dirty = subprocess.check_output(['git', 'status', '--porcelain', '--', *SOURCES], cwd=ROOT, text=True)
    if dirty:
        raise ValueError('commit the method sources before evaluating published observations')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    cache.mkdir(parents=True, exist_ok=True)
    record, elapsed = measure(cache, args.fetch, args.threads)
    record['method_commit'] = commit
    record['method_owner'] = 'docs/schur-memory.md#a-distance-observable-for-the-lambda-target'
    record['source_sha256'] = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}
    record['reproduce'] = ('python scripts/run_python.py scripts/measure_nsc_atlas_distance_contrast.py '
                           '--fetch --output .build/nsc-distance-contrast-reproduced.json')
    with output.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps(dict(summary=output.relative_to(ROOT).as_posix(), CPU_seconds=elapsed,
                          **record['statistics'])))


if __name__ == '__main__':
    main()

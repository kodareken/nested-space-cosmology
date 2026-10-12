"""Bind the imported actual18 rows to the local archive without a new solve.

The nine imported negative covariances differ from the local batch only by a
negative zero in ``C[0, 1].real``. Both raw digests are retained. The approved
per-row dyadics stay the unweighted errors. This module does not call the
vacuum solver.
"""
import hashlib, json, os, subprocess
from pathlib import Path
import numpy as np
from .nsc_evolved_incoming_state import _digest_arrays
from .nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from .nsc_ks_signed_state import s3_conjugate, source_complement_residual

IMPORT_RELATIVE = 'results/development/nsc-threshold18-import-v1'
CHECKPOINT_RELATIVE = IMPORT_RELATIVE + '/results/development/nsc-threshold18-actual18-v1'
IMPORT_TREE_SHA256 = 'bc7c5c2424810650acbeae20c4bedfc18548737f112725c63702eb3f9920a4c7'
BASE_COMMIT = 'eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9'
PROPAGATOR = 'src/recursive_horizons/nsc_ks_energy_propagator.py'
PROPAGATOR_BLOB = '31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227'
PROPAGATOR_SHA256 = '5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94'
PANEL = 'group14/low32_1'
ROWS = tuple(range(32, 41))
ADJUDICATION = 'c82b600b'


def _sha256(data):
    return hashlib.sha256(data).hexdigest()


def _float_view(array):
    array = np.ascontiguousarray(array)
    if np.iscomplexobj(array):
        return array.view(np.float64)
    return array


def _clear_negative_zeros(array):
    copied = np.array(array, copy=True, order='C')
    view = copied.view(np.float64) if np.iscomplexobj(copied) else copied
    view[np.signbit(view) & (view == 0)] = 0.0
    return copied


def require_signed_zero_only(local, imported, label):
    """Accept only value-identical arrays. A nonzero difference is not a signed zero."""
    local = np.ascontiguousarray(local)
    imported = np.ascontiguousarray(imported)
    if local.shape != imported.shape or local.dtype != imported.dtype:
        raise ValueError('physical source parameters differ: ' + label)
    if not np.array_equal(local, imported):
        difference = np.max(np.abs(local - imported))
        raise ValueError('numerical source arrays differ: %s max_abs=%s' % (label, difference))
    local_view, imported_view = _float_view(local), _float_view(imported)
    sign_diff = np.signbit(local_view) != np.signbit(imported_view)
    both_zero = (local_view == 0) & (imported_view == 0)
    nonzero = int((sign_diff & ~both_zero).sum())
    if nonzero:
        raise ValueError('numerical source arrays differ: ' + label)
    signed_zeros = int((sign_diff & both_zero).sum())
    if signed_zeros and not (local.tobytes() != imported.tobytes()):
        raise ValueError('signed-zero difference was not a raw byte difference: ' + label)
    return {'raw_equal': local.tobytes() == imported.tobytes(), 'signed_zero_differences': signed_zeros,
            'max_abs': 0.0, 'dtype': local.dtype.str, 'shape': list(local.shape)}


def _git_bytes(repo, *args):
    env = os.environ.copy()
    env['GIT_ALTERNATE_OBJECT_DIRECTORIES'] = str(repo / 'lab' / '.source-history' / 'objects')
    proc = subprocess.run(['git', '-C', str(repo), *args], capture_output=True, env=env)
    if proc.returncode != 0:
        raise ValueError('pinned git object unavailable: ' + ' '.join(args))
    return proc.stdout


def _file_sha(lab, relative):
    path = lab / relative
    if not path.is_file():
        raise ValueError('closure file missing: ' + relative)
    return _sha256(path.read_bytes())


def verify_closure(lab, record):
    """Check source, archive, baseline and successor hashes. The propagator is not HEAD."""
    repo = lab.parent
    if record.get('base_commit') != BASE_COMMIT:
        raise ValueError('authoritative base changed')
    source, archive_digests = record.get('source_hashes'), record.get('archive_digests')
    baseline, successor = record.get('baseline_family_source_hashes'), record.get('successor_hashes')
    reviewed = record.get('reviewed_snapshot_sha256')
    if not all(isinstance(item, dict) and item for item in (source, archive_digests, baseline, successor, reviewed)):
        raise ValueError('closure hash group missing')
    if len(source) != 67 or len(archive_digests) != 64:
        raise ValueError('closure hash census changed')
    for relative, expected in source.items():
        file_sha = _file_sha(lab, relative)
        git_sha = _sha256(_git_bytes(repo, 'show', 'HEAD:lab/' + relative))
        if file_sha != expected or git_sha != expected:
            raise ValueError('source closure hash mismatch: ' + relative)
    for relative, expected in archive_digests.items():
        if _file_sha(lab, relative) != expected:
            raise ValueError('archive closure hash mismatch: ' + relative)
    manifest = json.loads((lab / '.source-history' / 'manifest.json').read_text())
    pins = {(item.get('path'), item.get('sha256')): item for item in manifest.get('pins') or []}
    historical = []
    for relative, expected in baseline.items():
        file_sha = _file_sha(lab, relative)
        if file_sha == expected:
            continue
        pin = pins.get((relative, expected))
        if not isinstance(pin, dict):
            raise ValueError('baseline hash matches neither HEAD nor a pinned git object: ' + relative)
        blob = pin.get('blob')
        if _sha256(_git_bytes(repo, 'cat-file', 'blob', blob)) != expected:
            raise ValueError('pinned git object hash mismatch: ' + relative)
        historical.append({'path': relative, 'blob': blob, 'sha256': expected, 'head_sha256': file_sha})
    if historical != [{'path': PROPAGATOR, 'blob': PROPAGATOR_BLOB, 'sha256': PROPAGATOR_SHA256,
                       'head_sha256': _file_sha(lab, PROPAGATOR)}]:
        raise ValueError('baseline propagator must resolve as source-history blob 31818fdd, not HEAD')
    if _file_sha(lab, PROPAGATOR) == PROPAGATOR_SHA256:
        raise ValueError('HEAD propagator unexpectedly equals the historical hash')
    root = lab / IMPORT_RELATIVE
    for relative, expected in successor.items():
        if _sha256((root / relative).read_bytes()) != expected:
            raise ValueError('imported successor hash mismatch: ' + relative)
    for relative, expected in reviewed.items():
        if _file_sha(lab, relative) != expected or source.get(relative) != expected:
            raise ValueError('reviewed snapshot hash mismatch: ' + relative)
    published = _file_sha(lab, 'results/development/nsc-ks-source-operator-majorant-v3.json')
    if published != record.get('published_coverage_digest'):
        raise ValueError('published coverage digest mismatch')
    return {'source_hashes_checked': len(source), 'archive_digests_checked': len(archive_digests),
            'baseline_hashes_checked': len(baseline), 'historical_objects': historical,
            'successor_hashes_checked': len(successor), 'reviewed_snapshots_checked': len(reviewed),
            'published_coverage_digest': published, 'base_commit': BASE_COMMIT}


def import_tree_sha256(lab):
    root = lab / IMPORT_RELATIVE
    digest = hashlib.sha256()
    files = sorted(path for path in root.rglob('*') if path.is_file())
    for path in files:
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b'\0')
        digest.update(path.read_bytes())
        digest.update(b'\0')
    if digest.hexdigest() != IMPORT_TREE_SHA256 or len(files) != 25:
        raise ValueError('imported threshold18 tree changed')
    return digest.hexdigest()


def _negative_zero_at_corner(local, imported):
    """The only raw difference is local -0 versus imported +0 at C[0, 1].real."""
    if local.shape != (3, 3) or not np.iscomplexobj(local):
        raise ValueError('negative covariance fiber shape changed')
    local_real, imported_real = local[0, 1].real, imported[0, 1].real
    if not (local_real == 0 and imported_real == 0 and np.signbit(local_real)
            and not np.signbit(np.asarray(imported_real))):
        raise ValueError('negative covariance raw difference is not C[0, 1].real signed zero')
    view_local, view_imported = _float_view(local).ravel(), _float_view(imported).ravel()
    indices = np.flatnonzero(view_local.view(np.uint64) != view_imported.view(np.uint64))
    if indices.tolist() != [2]:
        raise ValueError('negative covariance has a raw difference outside C[0, 1].real')


def _batch(archive):
    positive = negative = None
    for batch, _channel in archive.family_entries((14, 1)):
        if batch.original_panel != PANEL or tuple(batch.rows) != (32, 48):
            continue
        if batch.energy_sign > 0 and batch.angular_sign > 0:
            positive = batch
        elif batch.energy_sign < 0 and batch.angular_sign < 0:
            negative = batch
    if positive is None or negative is None:
        raise ValueError('local group14/low32_1 batch 32:48 missing')
    return positive, negative


def replay_threshold18(lab, archive=None):
    """Compare all nine pairs and bind both negative raw digests. No ODE is run."""
    lab = Path(lab)
    if import_tree_sha256(lab) != IMPORT_TREE_SHA256:
        raise ValueError('imported threshold18 tree changed')
    archive = archive or RetainedUpstreamArchive(lab)
    positive, negative = _batch(archive)
    rows_root = lab / CHECKPOINT_RELATIVE / 'rows'
    loaded_rows = []
    closure = None
    witnessed_signed_zeros = 0
    for index in ROWS:
        directory = rows_root / f'row-{PANEL.replace("/", "__")}-{index:05d}'
        record_bytes = (directory / 'record.json').read_bytes()
        witness_bytes = (directory / 'witness.npz').read_bytes()
        record = json.loads(record_bytes)
        payload = record.get('payload') if isinstance(record, dict) else None
        if (not isinstance(payload, dict) or payload.get('sha256') != _sha256(witness_bytes)
                or payload.get('bytes') != len(witness_bytes)):
            raise ValueError('missing payload')
        if closure is None:
            closure = verify_closure(lab, record)
        elif loaded_rows:
            previous = loaded_rows[0][0]
            for field in ('source_hashes', 'archive_digests', 'baseline_family_source_hashes',
                          'successor_hashes', 'reviewed_snapshot_sha256', 'published_coverage_digest',
                          'base_commit', 'negative_source_digest', 'negative_preparation_digest',
                          'positive_source_digest', 'positive_preparation_digest'):
                if record.get(field) != previous.get(field):
                    raise ValueError('row closure hashes disagree')
        with np.load(directory / 'witness.npz', allow_pickle=False) as data:
            witness = {name: np.array(data[name]) for name in data.files}
        local_index = index - 32
        sl = slice(3 * local_index, 3 * local_index + 3)
        local_columns = positive.initial_columns[:, sl]
        local_negative_columns = negative.initial_columns[:, sl]
        local_covariance = positive.source.covariance[sl, sl]
        local_negative_covariance = negative.source.covariance[sl, sl]
        pairs = (
            ('positive_columns', local_columns, witness['positive_columns']),
            ('negative_columns', local_negative_columns, witness['negative_columns']),
            ('positive_covariance', local_covariance, witness['positive_covariance']),
            ('negative_covariance', local_negative_covariance, witness['negative_covariance']),
        )
        compared = {name: require_signed_zero_only(local, imported, f'{index}:{name}')
                    for name, local, imported in pairs}
        if any(compared[name]['signed_zero_differences'] for name in (
                'positive_columns', 'negative_columns', 'positive_covariance')):
            raise ValueError('signed zero appeared outside the negative covariance')
        if compared['negative_covariance']['signed_zero_differences'] != 1:
            raise ValueError('each negative covariance fiber must have one signed zero')
        _negative_zero_at_corner(local_negative_covariance, witness['negative_covariance'])
        witnessed_signed_zeros += 1
        if local_negative_columns.tobytes() != np.ascontiguousarray(s3_conjugate(local_columns)).tobytes():
            raise ValueError('negative columns are not the S3 conjugate')
        if source_complement_residual(witness['positive_covariance'], witness['negative_covariance']) != 0.0:
            raise ValueError('negative source differs from the fixed horizon source law')
        fiber = np.asarray(positive.source.energies[sl], float)
        negative_fiber = np.asarray(negative.source.energies[sl], float)
        weight = float(archive._arrays[PANEL + '/weights'][index])
        column_weight = float(positive.source.column_weights[sl][0])
        if (float(fiber[0]).hex() != record.get('energy_fiber_hex', [None])[0]
                or float(negative_fiber[0]).hex() != record.get('negative_energy_fiber_hex', [None])[0]
                or weight.hex() != record.get('weight_hex') or column_weight.hex() != record.get('column_weight_hex')
                or record.get('negative_endpoint_adjudication') != ADJUDICATION
                or record.get('negative_density_endpoint') != '(nx,-ny,-nz)'
                or record.get('negative_covariance_law') != 'Cminus=I-conj(Cplus)'):
            raise ValueError('physical source parameters differ')
        if (record.get('positive_source_digest') != positive.source.digest
                or record.get('positive_preparation_digest') != positive.preparation_digest):
            raise ValueError('positive source digest is not the local archive digest')
        if record.get('weight_applied_to_error') is not False:
            raise ValueError('source weight already applied')
        comparisons = record.get('signed_comparisons')
        if (not isinstance(comparisons, list) or len(comparisons) != 2
                or [item.get('energy_sign') for item in comparisons] != [1, -1]):
            raise ValueError('separate signed source bounds required')
        loaded_rows.append((record, index, fiber, negative_fiber))
    if len(loaded_rows) != 9 or witnessed_signed_zeros != 9:
        raise ValueError('actual18 row census changed')
    first = loaded_rows[0][0]
    if any(record.get('source_hashes') != first.get('source_hashes') for record, *_rest in loaded_rows):
        raise ValueError('row closure hashes disagree')
    covariance = negative.source.covariance
    weights = negative.source.column_weights
    energies = negative.source.energies
    if _digest_arrays(covariance, weights, energies) != negative.source.digest:
        raise ValueError('local negative digest is not the archive digest')
    negative_zero_count = int((np.signbit(_float_view(covariance)) & (_float_view(covariance) == 0)).sum())
    if negative_zero_count != 16:
        raise ValueError('negative batch signed-zero census changed')
    if int((np.signbit(weights) & (weights == 0)).sum()) or int((np.signbit(energies) & (energies == 0)).sum()):
        raise ValueError('signed zero appeared outside the negative covariance')
    cleared = _digest_arrays(_clear_negative_zeros(covariance), weights, energies)
    original_source = first['negative_source_digest']
    original_preparation = first['negative_preparation_digest']
    if cleared != original_source or cleared == negative.source.digest:
        raise ValueError('clearing negative zeros did not recover the imported source digest')
    if (first['negative_source_digest'] == negative.source.digest
            or first['negative_preparation_digest'] == negative.preparation_digest):
        raise ValueError('negative raw digests unexpectedly matched')
    rows = []
    for record, index, fiber, negative_fiber in loaded_rows:
        owner = f'{CHECKPOINT_RELATIVE}/rows/row-{PANEL.replace("/", "__")}-{index:05d}/record.json'
        comparisons = []
        for sign, item, energy in ((1, record['signed_comparisons'][0], fiber),
                                   (-1, record['signed_comparisons'][1], negative_fiber)):
            if item.get('signed_total_upper') is None:
                raise ValueError('missing signed total')
            if item.get('weight_applied_to_error') is not False:
                raise ValueError('source weight already applied')
            source_digest = positive.source.digest if sign > 0 else negative.source.digest
            preparation = positive.preparation_digest if sign > 0 else negative.preparation_digest
            original_source_digest = record['positive_source_digest'] if sign > 0 else original_source
            original_preparation_digest = (record['positive_preparation_digest'] if sign > 0
                                           else original_preparation)
            comparisons.append({
                'energy_sign': sign, 'angular_sign': sign, 'energy_hex': float(energy[0]).hex(),
                'signed_total_upper': item['signed_total_upper'], 'weight_applied_to_error': False,
                'local_raw_source_digest': source_digest, 'local_raw_preparation_digest': preparation,
                'original_raw_source_digest': original_source_digest,
                'original_raw_preparation_digest': original_preparation_digest})
        rows.append({'panel': PANEL, 'row': index, 'checkpoint': owner, 'weight_applied_to_error': False,
                     'signed_comparisons': comparisons})
    equivalence = {
        'relation': 'signed_zero_only',
        'location': 'negative covariance C[0, 1].real, one negative zero per fiber',
        'original_raw_source_digest': original_source,
        'local_raw_source_digest': negative.source.digest,
        'original_raw_preparation_digest': original_preparation,
        'local_raw_preparation_digest': negative.preparation_digest,
        'positive_source_digest': positive.source.digest,
        'positive_preparation_digest': positive.preparation_digest,
        'positive_digests_match': True,
        'batch_signed_zeros': negative_zero_count,
        'witnessed_fiber_signed_zeros': witnessed_signed_zeros,
        'nonzero_numeric_differences': 0,
        'weights_and_energies_signed_zeros': 0,
        'clearing_local_negative_zeros_reproduces_original_source_digest': True,
        'propagator_resolution': closure['historical_objects'][0]}
    return {'schema': 'NSC-THRESHOLD18-REPLAY-EQUIVALENCE-v1', 'mode': 'check', 'new_rows_captured': 0,
            'ode_rerun': False, 'uses_imported_dyadic': True, 'weight_applied_to_error': False,
            'coverage_complete': True, 'completed_positive_rows': 9, 'completed_signed_rows': 18,
            'physical_local_gate': 'OPEN', 'physical_upstream_budget_component': None,
            'negative_endpoint_adjudication': ADJUDICATION,
            'negative_density_endpoint': '(nx,-ny,-nz)', 'rows': rows, 'equivalence': equivalence,
            'closure': closure, 'import_tree_sha256': IMPORT_TREE_SHA256}

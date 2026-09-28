"""Actual retained source columns for the source-fixed KS history evaluator.

This reader reuses physical low fields and the baseline's original middle
producer. It does not turn moments, contours or a frozen covariance into modes.
"""
from contextlib import ExitStack
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.integrate import DOP853

from .nsc_common_ks_trace import restrict_resolved_modes
from .nsc_evolved_incoming_state import FixedSourcePreparation, _finite_array
from .nsc_incoming_source_assembly import read_panels, selected_pieces
from .nsc_incoming_spectral_panels import RECORDS, record, sha
from .nsc_paired_horizon_preparation import source_covariance
from .nsc_pg_archived_high_energy_modes import archived_middle_boundary_modes
from .nsc_retarded_radial_response import ks_generator
from .nsc_transmitting_dirac_domain import S3


@dataclass(frozen=True)
class ReferenceSourcePanel:
    name: str
    group: int
    angular_sign: int
    mass: float
    angular_magnitude: float
    energies: object
    weights: object
    covariance: object
    amplitudes_at_one: object
    provenance: object

    def __post_init__(self):
        for name, dtype in (('energies', float), ('weights', float),
                            ('covariance', complex), ('amplitudes_at_one', complex)):
            raw = getattr(self, name)
            if dtype is float and np.iscomplexobj(raw):
                raise ValueError('physical source energies and weights must be real')
            object.__setattr__(self, name, _finite_array(raw, dtype, name))
        n = len(self.energies)
        if (self.energies.shape != (n,) or self.weights.shape != (n,) or n == 0
                or np.any(self.energies == 0) or np.any(self.weights <= 0)
                or self.covariance.shape != (n, 3, 3) or self.amplitudes_at_one.shape != (n, 2, 3)):
            raise ValueError('complete nonzero real three-column source fibers required')
        if self.angular_sign not in (-1, 1) or self.mass < 0 or self.angular_magnitude < 0:
            raise ValueError('actual signed channel required')
        if np.max(abs(self.covariance - self.covariance.swapaxes(-1, -2).conj())) > 3e-13:
            raise ValueError('Hermitian source covariance required')
        eig = np.linalg.eigvalsh(self.covariance)
        if eig.min() < -3e-13 or eig.max() > 1 + 3e-13:
            raise ValueError('source CAR interval violated')

    @property
    def angular(self):
        return self.angular_sign * self.angular_magnitude

    def batches(self, energy_count):
        """Keep each three-column coherent fiber intact when limiting memory."""
        if isinstance(energy_count, bool) or int(energy_count) != energy_count or energy_count < 1:
            raise ValueError('positive integer energy batch size required')
        for start in range(0, len(self.energies), int(energy_count)):
            end = min(start + int(energy_count), len(self.energies))
            yield ReferenceSourcePanel(self.name + f'/rows{start}:{end}', self.group,
                self.angular_sign, self.mass, self.angular_magnitude, self.energies[start:end],
                self.weights[start:end], self.covariance[start:end], self.amplitudes_at_one[start:end],
                {'parent_panel': self.name, 'rows': [start, end], 'source': self.provenance})

    def negative_partner(self, config):
        """Owned opposite-angular antiunitary map; no signed-energy folding."""
        if np.any(self.energies <= 0):
            raise ValueError('positive input panel required for the negative partner')
        blocks = np.array([source_covariance(-e, config['surface_gravity'], config['omega'], self.mass)
                           for e in self.energies])
        if np.max(abs(blocks - (np.eye(3) - self.covariance.conj()))) > 3e-13:
            raise ValueError('panel differs from the fixed horizon source law')
        return ReferenceSourcePanel(self.name + '/negative-partner', self.group,
            -self.angular_sign if self.angular_magnitude else 1, self.mass, self.angular_magnitude,
            -self.energies, self.weights, blocks, S3 @ self.amplitudes_at_one.conj(),
            {'positive_panel': self.name, 'map': 'S3 conjugation with opposite angular sign',
             'source': self.provenance, 'energy_folding_factor': 1})

    def fixed_upstream(self, *, rho_up, rtol, atol, max_step):
        """Reference continuation once, then immutable Cup for all trial histories.

        This is a numerical mode continuation, not a certificate of the
        source approximation. Prefer small complete-energy batches.
        """
        if not np.isfinite([rho_up, rtol, atol, max_step]).all() or rho_up < 1.03 or min(rtol, atol, max_step) <= 0:
            raise ValueError('explicit rho_up>=1.03 and positive continuation controls required')
        initial = self.amplitudes_at_one

        def rhs(rho, state):
            return (ks_generator(rho, self.energies, self.mass, self.angular)
                    @ state.reshape(initial.shape)).ravel()

        solver = DOP853(rhs, 1., initial.ravel(), float(rho_up), rtol=rtol, atol=atol, max_step=max_step)
        while solver.status == 'running':
            message = solver.step()
            if solver.status == 'failed':
                raise ArithmeticError(message)
        amplitudes = solver.y.reshape(initial.shape).transpose(1, 0, 2).reshape(2, -1)
        source = FixedSourcePreparation.from_signed_blocks(self.energies, self.weights, self.covariance)
        return source, amplitudes, {
            'panel': self.name, 'rho_up': float(rho_up), 'nfev': int(solver.nfev),
            'source_digest': source.digest, 'source_provenance': self.provenance,
            'continuation_error_bound': None, 'physical_source_error_bound': None,
            'physical_gate': 'OPEN', 'incoming_state_frozen': False,
        }


class RetainedSourceInventory:
    """Context-managed authenticated reader; no field or history solve on load."""
    def __init__(self, root):
        self.root = Path(root)
        self._stack = ExitStack()
        try:
            self.panels, self.metadata = read_panels(self.root)
            receipts = {p: record(self.root, p) for p in RECORDS}
            self.provenance = []

            def opened(spec):
                if sha(self.root, spec['path']) != spec['sha256']:
                    raise ValueError('source artifact changed: ' + spec['path'])
                self.provenance.append(spec)
                return self._stack.enter_context(np.load(self.root / spec['path'], allow_pickle=False))

            self.recovery, self.retained, self.g13, self.mode = [opened(receipts[p]['payload']) for p in RECORDS[:4]]
            self.g14 = opened(json.loads(self.retained['metadata_json'].tobytes())['group14_input'])
            self.channels = receipts[RECORDS[4]]['channels']
            self.config = receipts[RECORDS[5]]['scattering_provenance']['config']
            self._g13meta = json.loads(self.g13['metadata_json'].tobytes())
        except BaseException:
            self._stack.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self._stack.close()

    def _middle(self, name):
        if name.startswith(('mid/', 'mid_ref/')):
            return json.loads(self.retained[name + '/metadata_json'].tobytes())
        if name.startswith('group13/mid'):
            m = self._g13meta
            return {'sources': m['source_hashes'], 'channel': m['channel'], 'angular_sign': 1,
                    'left': 16., 'right': 160., 'config': m['config']}
        if name.startswith('group14/mid'):
            return json.loads(self.g14[name.split('/')[1] + '/metadata_json'].tobytes())
        return None

    def positive_panels(self, groups=None):
        """Disjoint original quadratures; exact ell=0 response is handled separately."""
        groups = range(1, 33) if groups is None else tuple(groups)
        if any(g not in range(1, 33) for g in groups):
            raise ValueError('retained real-field groups are 1 through 32; LLL has its analytic owner')
        for group in groups:
            channel = self.channels[group]
            for sign in ((1,) if channel['angular_eigenvalue'] == 0 else (-1, 1)):
                pieces = selected_pieces(self.panels, self.metadata['panels'], group, sign)
                for name, mask in pieces:
                    E, weights = (self.panels[name + '/' + k][mask] for k in ('energies', 'weights'))
                    middle = self._middle(name)
                    if middle is not None:
                        phi = archived_middle_boundary_modes(E, middle, repo_root=self.root).mode_at_one
                        C = np.array([source_covariance(e, self.config['surface_gravity'], self.config['omega'],
                                                       channel['compact_mass']) for e in E])
                        method = 'same original middle producer as the incoming baseline'
                    elif name.startswith('subgap8/'):
                        label = self.mode['label']
                        ix = np.flatnonzero((label[:, 0] == group) & (label[:, 1] == sign)
                            & (label[:, 5] > 1) & (label[:, 5] < label[:, 3]))
                        at_one = int(np.flatnonzero(self.mode['rho_interior'] == 1.)[0])
                        if len(ix) != 8 or not np.array_equal(label[ix, 5], E):
                            raise ValueError('coarse subgap energy inventory changed')
                        phi, C = self.mode['interior'][ix, at_one], self.mode['source_covariance'][ix]
                        method = 'eight inherited real subgap rows; accuracy OPEN'
                    else:
                        if not np.array_equal(self.recovery[name + '/energies'][mask], E):
                            raise ValueError('recovered source energies differ from selected panel')
                        if not self.recovery[name + '/accepted'][mask].all():
                            raise ValueError('unresolved physical low source fields')
                        phi, C = (self.recovery[name + '/' + k][mask] for k in ('mode_at_one', 'source'))
                        method = 'accepted physical low-mode reconstruction'
                    canonical, _ = restrict_resolved_modes(1., E, phi)
                    yield ReferenceSourcePanel(name, group, sign, channel['compact_mass'],
                        channel['angular_eigenvalue'], E, weights, C, canonical,
                        {'input_payloads': self.provenance, 'selection': name,
                         'selection_mask_sha256': sha256(mask.tobytes()).hexdigest(),
                         'method': method, 'source_error_bound': None, 'infinite_tail_included': False})

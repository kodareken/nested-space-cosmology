"""Restore the original retained upstream batches without repeating preparation."""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .nsc_evolved_incoming_state import FixedSourcePreparation
from .nsc_ks_batched_constraints import KSUpstreamBatch, map_negative_batch
from .nsc_ks_source_inventory import ReferenceSourcePanel


class RetainedUpstreamArchive:
    """Authenticated original A_up/C_src/weights, independent of trial history.

    The archived errors remain the archived errors. Restoring these arrays
    supplies no new preparation or low/subgap accuracy certificate.
    """
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.input_hashes = {}
        inventory = self._record('results/development/nsc-ks-source-inventory.json')
        aggregate = self._record('results/development/nsc-ks-retained-response-sum.json')
        if not aggregate['all_retained_nonzero_angular_families_included']:
            raise ValueError('complete retained family archive required')
        with np.load(self._payload(inventory['payload']), allow_pickle=False) as f:
            self._arrays = {k: np.array(f[k], copy=True) for k in f.files}
        self.meta = json.loads(self._arrays['metadata_json'].tobytes())
        self._families = {}
        for item in aggregate['family_records']:
            record = self._record(item['path'], item['sha256'])
            key = tuple(record['positive_family'])
            if key in self._families:
                raise ValueError('duplicate original family archive')
            self._families[key] = record
        if len(self._families) != aggregate['required_positive_families']:
            raise ValueError('missing original family archive')
        self._panels = {}

    def _file(self, name, expected=None):
        path = (self.root/name).resolve()
        if not path.is_relative_to(self.root):
            raise ValueError('source archive input must be inside its declared repository')
        value = sha256(path.read_bytes()).hexdigest()
        if expected is not None and value != expected:
            raise ValueError('source archive hash mismatch: '+str(name))
        self.input_hashes[str(path.relative_to(self.root))] = value
        return path

    def _record(self, name, expected=None):
        return json.loads(self._file(name, expected).read_text())

    def _payload(self, description):
        return self._file(description['path'], description['sha256'])

    @property
    def family_keys(self):
        return tuple(sorted(self._families))

    def _panel_pair(self, name):
        if name not in self._panels:
            p = self.meta['panels'][name]
            panel = ReferenceSourcePanel(name, p['group'], p['angular_sign'], p['mass'],
                p['angular_magnitude'], *(self._arrays[name+'/'+key] for key in
                ('energies', 'weights', 'covariance', 'amplitudes_at_one')), p['provenance'])
            self._panels[name] = panel, panel.negative_partner(self.meta['config'])
        return self._panels[name]

    def family_entries(self, key):
        """Original positive family plus explicit opposite-angular negative partner."""
        key = tuple(key)
        if key not in self._families:
            raise ValueError('retained positive angular family required')
        group, sign = key
        record = self._families[key]
        channel = self.meta['channels'][group]
        if channel['index'] != group:
            raise ValueError('original channel ledger index changed')
        with np.load(self._payload(record['payload']), allow_pickle=False) as f:
            receipts = json.loads(f['metadata_json'].tobytes())['upstream_batches']
            entries = []
            for receipt in receipts:
                if (receipt['group'], receipt['angular_sign'], receipt['energy_sign']) != (group, sign, 1):
                    raise ValueError('original upstream family binding mismatch')
                positive_panel, negative_panel = self._panel_pair(receipt['original_panel'])
                lo, hi = receipt['rows']
                def source(panel):
                    return FixedSourcePreparation.from_signed_blocks(panel.energies[lo:hi],
                        panel.weights[lo:hi], panel.covariance[lo:hi])
                positive = KSUpstreamBatch(**receipt, source=source(positive_panel),
                    initial_columns=f['upstream/'+receipt['panel_name']])
                negative = map_negative_batch(positive, source(negative_panel))
                entries.extend(((positive, {**channel, 'angular_sign': sign}),
                                (negative, {**channel, 'angular_sign': -sign})))
        if sum(len(batch.source.energies)//3 for batch, _ in entries if batch.energy_sign > 0) != record['source_rows']:
            raise ValueError('restored rows differ from the original family coverage')
        return tuple(entries)

    def all_entries(self):
        return tuple(entry for key in self.family_keys for entry in self.family_entries(key))

    def analytic_zero_channels(self):
        return tuple({**self.meta['channels'][group], 'angular_sign': 1}
                     for group in (0, 13, 23))

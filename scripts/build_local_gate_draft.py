#!/usr/bin/env python3
"""Build two identical PDFs and a source archive for the OPEN research draft."""
import argparse
import json
from local_gate_draft import (ROOT, EVIDENCE, SOURCE, PDF, ARCHIVE, MANIFEST,
                              authenticate_snapshot, source_files, compile_twice,
                              packed_source, digest, verify)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='compare with tracked artifacts without writing')
    args = parser.parse_args()
    snapshot = authenticate_snapshot(ROOT)
    files = source_files(ROOT)
    pdf, bbl, compiler = compile_twice(files)
    archive = packed_source({**files, 'main.bbl': bbl})
    record = {
        'schema': 'NSC-OPEN-DRAFT-BUILD-v1', 'status': 'OPEN', 'submission_ready': False,
        'lab_commit': snapshot['lab_commit'], 'complete_scientific_dependency_closure': False,
        'compiler': compiler, 'source_date_epoch': 0, 'clean_builds': 2,
        'byte_identical_clean_builds': True, 'pdf_sha256': digest(pdf),
        'bbl_sha256': digest(bbl), 'source_archive_sha256': digest(archive),
        'texlive_container_digest': None, 'arxiv_runtime_verified': False,
        'license_selected_for_arxiv': None,
        'inputs': {**{(SOURCE / n).as_posix(): digest(b) for n, b in files.items()},
                   (EVIDENCE / 'snapshot.json').as_posix(): digest((ROOT / EVIDENCE / 'snapshot.json').read_bytes()),
                   'scripts/local_gate_draft.py': digest((ROOT / 'scripts/local_gate_draft.py').read_bytes())},
    }
    if args.check:
        verify()
        if pdf != (ROOT / PDF).read_bytes() or archive != (ROOT / ARCHIVE).read_bytes():
            raise ValueError('rebuilt artifact differs; use the recorded compiler environment')
    else:
        (ROOT / PDF).write_bytes(pdf)
        (ROOT / ARCHIVE).write_bytes(archive)
        (ROOT / MANIFEST).write_text(json.dumps(record, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'ok': True, 'status': 'OPEN', 'pdf_sha256': digest(pdf),
                      'identical_clean_builds': True, 'submission_ready': False}, indent=2))


if __name__ == '__main__':
    main()

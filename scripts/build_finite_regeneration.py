#!/usr/bin/env python3
"""Build the finite companion from an explicit immutable science snapshot."""
import argparse
import json

from finite_regeneration import build, compile_twice, ROOT, source_files, verify

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--science-commit")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--proof-only", action="store_true")
    args = parser.parse_args()
    if args.check:
        print(json.dumps(verify(), indent=2))
    elif args.proof_only:
        pdf, _, _, compiler = compile_twice(source_files())
        path = ROOT / ".build/finite-regeneration/proof.pdf"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(pdf)
        print(path, compiler)
    else:
        if not args.science_commit:
            parser.error("--science-commit requires a full immutable commit")
        print(json.dumps(build(args.science_commit), indent=2))

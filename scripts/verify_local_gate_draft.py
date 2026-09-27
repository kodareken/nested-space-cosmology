#!/usr/bin/env python3
"""Check draft hashes and OPEN scope without rerunning scientific campaigns."""
import json
from local_gate_draft import verify

if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))

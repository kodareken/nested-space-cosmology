#!/usr/bin/env python3
"""Check draft hashes and OPEN scope without rerunning scientific campaigns."""
import json
from replay_history import run

if __name__ == '__main__':
    print(json.dumps(run("companion"), indent=2, sort_keys=True))

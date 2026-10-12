#!/usr/bin/env python3
"""Read-only authentication of the focused finite-regeneration publication."""
import json
from replay_history import run

if __name__ == "__main__":
    print(json.dumps(run("finite"), indent=2))

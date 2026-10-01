#!/usr/bin/env python3
"""Read-only authentication of the focused finite-regeneration publication."""
import json
from finite_regeneration import verify

if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))

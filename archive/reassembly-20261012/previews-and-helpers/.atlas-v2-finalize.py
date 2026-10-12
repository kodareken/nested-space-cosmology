from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent / "BlackHoles-Infinity"

def names():
    changed = subprocess.check_output(["git", "diff", "--name-only", "--", "."], cwd=ROOT, text=True).splitlines()
    added = subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard"], cwd=ROOT, text=True).splitlines()
    return sorted(n for n in set(changed + added) if not n.startswith(".pi/"))

if sys.argv[1] == "normalize":
    for name in names():
        p = ROOT / name
        if not name.startswith("lab/") and p.suffix in (".md", ".json", ".py"):
            p.write_bytes(p.read_bytes().replace(b"\r\n", b"\n"))
    print("Normalized task-owned publication text to LF")
elif sys.argv[1] == "snapshot":
    lab_paths = [n for n in names() if n.startswith("lab/")]
    subprocess.run([sys.executable, "scripts/update_lab_snapshot.py", *lab_paths], cwd=ROOT, check=True)
    p = ROOT / "docs/lab-snapshot.json"
    p.write_bytes(p.read_bytes().replace(b"\r\n", b"\n"))
    print(json.dumps({"registered_lab_paths": lab_paths}, indent=2))
elif sys.argv[1] == "links":
    sys.path.insert(0, str(ROOT / "scripts"))
    from check_publication import check_markdown_links
    errors = []
    for name in names():
        if name.endswith(".md"):
            p = ROOT / name
            check_markdown_links(errors, p, p.read_text(encoding="utf-8"))
    print(json.dumps({"link_errors": errors}, indent=2))
    raise SystemExit(bool(errors))
elif sys.argv[1] == "scope":
    paths = names()
    print(json.dumps({"task_paths": paths}, indent=2))

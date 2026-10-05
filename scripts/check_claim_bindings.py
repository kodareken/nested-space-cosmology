"""Check the paper's claims against their records, and the findings register.

Two questions, both made mechanical instead of remembered.

First, structure. Every <!-- nsc-claim:ID --> marker in the foundation paper
must have an entry in paper/claim-bindings.json, and every entry must name a
marker that exists. Every record and note the paper links from a claim must
resolve to a file. Every status must be one the binding file declares.

Second, coverage. Every finding the laboratory registers must name a status,
a note, and either a record or an explicit null with the reason carried by
post_preprint. Findings whose paper_section is null are reported as
post-preprint, which is the answer to "is this in the paper".

Exit status is non-zero only for a structural failure. An unassigned status
or a post-preprint finding is reported as a count, not as an error, because
neither is wrong on its own.

The record links are read from the paper, not from the binding file, so the
two cannot silently disagree.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

HEADING = re.compile(r"^#{2,4}\s+(.*)$")
MARKER = re.compile(r"<!--\s*nsc-claim:([a-z0-9-]+)\s*-->")
LINK = re.compile(r"\]\(([^)]+)\)")


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def read_paper(path: pathlib.Path):
    """Return {claim_id: {...}} with the links that follow each marker."""
    lines = path.read_text(encoding="utf-8").splitlines()
    marks: list[tuple[int, str, str]] = []
    heading = "(preamble)"
    for number, line in enumerate(lines):
        head = HEADING.match(line)
        if head:
            heading = head.group(1).strip()
        for found in MARKER.finditer(line):
            marks.append((number, found.group(1), heading))

    claims: dict[str, dict] = {}
    for index, (start, claim_id, heading) in enumerate(marks):
        if claim_id in claims:
            raise ValueError(f"duplicate claim marker: {claim_id}")
        stop = marks[index + 1][0] if index + 1 < len(marks) else len(lines)
        body = "\n".join(lines[start:stop])
        records, notes = [], []
        for target in LINK.findall(body):
            if target.startswith(("http://", "https://")):
                continue
            target = target.split("#")[0]
            if target.endswith((".json", ".npz")):
                records.append(target)
            elif target.endswith(".md"):
                notes.append(target)
        claims[claim_id] = {
            "line": start + 1,
            "heading": heading,
            "records": sorted(set(records)),
            "notes": sorted(set(notes)),
        }
    return claims


def resolve(root: pathlib.Path, link: str) -> pathlib.Path:
    """Paper links are written relative to the paper directory."""
    return (root / "paper" / link).resolve()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=None)
    args = parser.parse_args(argv)

    root = pathlib.Path(args.root).resolve() if args.root else repo_root()
    paper_path = root / "paper/nested-space-cosmology.md"
    binding_path = root / "paper/claim-bindings.json"
    if not paper_path.is_file():
        raise SystemExit(f"paper not found: {paper_path}")
    if not binding_path.is_file():
        raise SystemExit(f"binding file not found: {binding_path}")

    binding = json.loads(binding_path.read_text(encoding="utf-8"))
    allowed = set(binding["status_values"])
    claims = read_paper(paper_path)
    declared = binding["claims"]

    errors: list[str] = []

    missing = sorted(set(claims) - set(declared))
    extra = sorted(set(declared) - set(claims))
    for claim_id in missing:
        errors.append(f"paper marker has no binding entry: {claim_id}")
    for claim_id in extra:
        errors.append(f"binding entry has no paper marker: {claim_id}")

    by_status: dict[str, int] = {}
    no_record: list[str] = []
    record_count = 0
    for claim_id, row in sorted(claims.items()):
        status = declared.get(claim_id, {}).get("status", "<absent>")
        by_status[status] = by_status.get(status, 0) + 1
        if status not in allowed:
            errors.append(f"invalid status for {claim_id}: {status!r}")
        if not row["records"]:
            no_record.append(claim_id)
        for link in row["records"]:
            record_count += 1
            if not resolve(root, link).is_file():
                errors.append(f"{claim_id}: record does not exist: {link}")
        for link in row["notes"]:
            if not resolve(root, link).is_file():
                errors.append(f"{claim_id}: note does not exist: {link}")

    findings = binding.get("findings", [])
    post_preprint: list[str] = []
    in_paper: list[str] = []
    pending: list[str] = []
    for item in findings:
        fid = item.get("id", "<unnamed>")
        if item.get("status") not in allowed:
            errors.append(f"finding {fid}: invalid status {item.get('status')!r}")
        note = item.get("note")
        if note and not (root / note).is_file():
            errors.append(f"finding {fid}: note does not exist: {note}")
        record = item.get("record")
        if record and not (root / record).is_file():
            errors.append(f"finding {fid}: record does not exist: {record}")
        draft = item.get("pending_section")
        if draft:
            if not (root / draft).is_file():
                errors.append(f"finding {fid}: pending section does not exist: {draft}")
            pending.append(fid)
        if item.get("paper_section"):
            in_paper.append(fid)
        else:
            post_preprint.append(fid)

    print(f"paper   {paper_path.relative_to(root)}")
    print(f"binding {binding_path.relative_to(root)}")
    print()
    print(f"claims                {len(claims)}")
    print(f"record links resolved {record_count}")
    for status, count in sorted(by_status.items(), key=lambda kv: -kv[1]):
        print(f"    {count:3d}  {status}")
    print()
    print(f"claims with no record link  {len(no_record)}")
    for claim_id in no_record:
        print(f"    {claim_id}")
    print()
    print(f"findings registered   {len(findings)}")
    print(f"    in a paper section  {len(in_paper)}")
    print(f"    post-preprint       {len(post_preprint)}")
    for fid in post_preprint:
        mark = " (draft written)" if fid in pending else ""
        print(f"        {fid}{mark}")

    if errors:
        print()
        print(f"FAILURES ({len(errors)})")
        for message in errors:
            print(f"    {message}")
        return 1

    print()
    print("claim bindings check passed")
    if by_status.get("unspecified"):
        print(f"note: {by_status['unspecified']} claims still have an unassigned status")
    return 0


if __name__ == "__main__":
    sys.exit(main())
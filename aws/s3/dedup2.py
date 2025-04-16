#!/usr/bin/env python3

import os
from pathlib import Path
from collections import defaultdict

base_dir = Path.home() / "data"
dupes_dir = base_dir / "duplicates"
rnaseq_dir = base_dir / "rnaseq"
tapestri_dir = base_dir / "tapestri"

def get_file_set(directory):
    return {f.name: f.stat().st_size for f in directory.rglob("*") if f.is_file()}

report = []

for dup_dir in dupes_dir.iterdir():
    if not dup_dir.is_dir():
        continue

    target_dir = tapestri_dir / dup_dir.name if dup_dir.name.endswith("_Mission_Bio") else rnaseq_dir / dup_dir.name
    dup_files = get_file_set(dup_dir)

    if target_dir.exists():
        target_files = get_file_set(target_dir)

        match_files = sorted(set(dup_files) & set(target_files))
        extra_in_dup = sorted(set(dup_files) - set(target_files))
        extra_in_target = sorted(set(target_files) - set(dup_files))

        size_mismatch = [
            f for f in match_files if dup_files[f] != target_files[f]
        ]
        match_files = [f for f in match_files if f not in size_mismatch]

        report.append(f"✓ {dup_dir.name} matches {target_dir.name}" if not extra_in_dup and not extra_in_target and not size_mismatch else f"⚠ {dup_dir.name} differs from {target_dir.name}")
        report.append(f"  ✔ Matched files: {len(match_files)}")
        if size_mismatch:
            report.append(f"  ✗ Size mismatches: {', '.join(size_mismatch)}")
        if extra_in_dup:
            report.append(f"  ✗ Extra in duplicates: {', '.join(extra_in_dup)}")
        if extra_in_target:
            report.append(f"  ✗ Extra in target: {', '.join(extra_in_target)}")
    else:
        report.append(f"✗ {dup_dir.name} not found in target directory")

print("\n".join(report))

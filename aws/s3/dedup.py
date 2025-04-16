#!/usr/bin/env python3

import os
import shutil
from pathlib import Path

base_dir = Path.home() / "data"
rnaseq_dir = base_dir / "rnaseq"
tapestri_dir = base_dir / "tapestri"
dupes_dir = base_dir / "duplicates"

dupes_dir.mkdir(exist_ok=True)


def get_sample_dirs(root):
    return {p.name: p for p in root.iterdir() if p.is_dir()}


def move_and_cleanup(src, dest_root):
    dest = dest_root / src.name
    print(f"Moving {src} to {dest}")
    shutil.move(str(src), dest)
    parent = src.parent
    if not any(parent.iterdir()):
        print(f"Removing empty directory {parent}")
        parent.rmdir()


rnaseq_samples = get_sample_dirs(rnaseq_dir)
tapestri_samples = get_sample_dirs(tapestri_dir)

# From rnaseq
for name, path in rnaseq_samples.items():
    if name.endswith("_Mission_Bio") and name in tapestri_samples:
        move_and_cleanup(path, dupes_dir)

# From tapestrI
for name, path in tapestri_samples.items():
    if not name.endswith("_Mission_Bio") and name in rnaseq_samples:
        move_and_cleanup(path, dupes_dir)


#!/usr/bin/env python3
"""Regression test: LIIC failure should trigger DM fallback (when enabled).

We intentionally force LIIC to fail deterministically by setting a huge
--ev-thresh so no DLC eigenvectors are selected.

We then assert that the program prints the DM-fallback message and still
produces an output trajectory.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: test_dm_fallback.py <neb_interpolator> <init.xyz> <final.xyz>")
        return 2

    exe = Path(sys.argv[1])
    init_xyz = Path(sys.argv[2])
    final_xyz = Path(sys.argv[3])

    for p in (exe, init_xyz, final_xyz):
        if not p.exists():
            print(f"FAIL: missing path: {p}")
            return 2

    prefix = "dmfb_"
    out_xyz = Path(prefix + "trajectory.xyz")
    if out_xyz.exists():
        out_xyz.unlink()

    cmd = [
        str(exe),
        "--no-align",
        "-o",
        "multiframe",
        "-p",
        prefix,
        "-n",
        "3",
        "-m",
        "liic",
        "--ev-thresh",
        "1e9",
        str(init_xyz),
        str(final_xyz),
    ]

    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = p.stdout

    if p.returncode != 0:
        print(out)
        print(f"FAIL: neb_interpolator returned {p.returncode}")
        return 1

    if "Falling back to distance-matrix interpolation (DM)" not in out:
        print(out)
        print("FAIL: did not observe DM fallback message")
        return 1

    if not out_xyz.exists():
        print(out)
        print(f"FAIL: expected output trajectory not found: {out_xyz}")
        return 1

    print("PASS: DM fallback triggered and output written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

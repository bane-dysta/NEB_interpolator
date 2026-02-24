#!/usr/bin/env python3
"""Regression test: external engine interface (I/O + parsing + assembly).

We run NEB with --engine-cmd pointing to mock_engine_constant.py.
The mock produces zero gradients so the NEB should converge immediately.

We assert that the run succeeds and writes an output trajectory.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "Usage: test_external_engine.py <neb_interpolator> <mock_engine.py> <init.xyz> <final.xyz>"
        )
        return 2

    exe = Path(sys.argv[1])
    engine = Path(sys.argv[2])
    init_xyz = Path(sys.argv[3])
    final_xyz = Path(sys.argv[4])

    for p in (exe, engine, init_xyz, final_xyz):
        if not p.exists():
            print(f"FAIL: missing path: {p}")
            return 2

    prefix = "engine_"
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
        "2",
        "-m",
        "neb",
        "-i",
        "50",
        "-c",
        "1e-6",
        "--engine-cmd",
        f"python3 {engine}",
        str(init_xyz),
        str(final_xyz),
    ]

    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = p.stdout

    if p.returncode != 0:
        print(out)
        print(f"FAIL: neb_interpolator returned {p.returncode}")
        return 1

    if "[engine]" not in out:
        print(out)
        print("FAIL: did not observe engine invocation logs")
        return 1

    if not out_xyz.exists():
        print(out)
        print(f"FAIL: expected output trajectory not found: {out_xyz}")
        return 1

    print("PASS: external engine mode ran and output written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

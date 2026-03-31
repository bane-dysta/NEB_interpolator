#!/usr/bin/env python3
"""Regression test: LIIC failure should trigger DM fallback and be reported explicitly.

We intentionally force LIIC to fail deterministically by setting a huge
--ev-thresh so no DLC eigenvectors are selected.

The CLI should then:
1. generate a path via DM,
2. print a run summary exposing the fallback, and
3. optionally return exit code 2 when --strict-path-method is requested.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run_case(exe: Path, init_xyz: Path, final_xyz: Path, prefix: str, extra_args: list[str]) -> subprocess.CompletedProcess[str]:
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
        *extra_args,
        str(init_xyz),
        str(final_xyz),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if not out_xyz.exists():
        print(proc.stdout)
        print(f"FAIL: expected output trajectory not found: {out_xyz}")
        raise SystemExit(1)
    return proc


def require_contains(output: str, needle: str, label: str) -> None:
    if needle not in output:
        print(output)
        print(f"FAIL: missing {label}: {needle}")
        raise SystemExit(1)


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

    normal = run_case(exe, init_xyz, final_xyz, "dmfb_", [])
    if normal.returncode != 0:
        print(normal.stdout)
        print(f"FAIL: neb_interpolator returned {normal.returncode} for fallback run")
        return 1

    require_contains(normal.stdout, "Status: GENERATED_WITH_FALLBACK", "fallback status")
    require_contains(normal.stdout, "Requested path method: LIIC", "requested path method")
    require_contains(normal.stdout, "Actual path method: DM", "actual path method")
    require_contains(normal.stdout, "Attempted path methods: LIIC -> DM", "attempt chain")

    strict = run_case(exe, init_xyz, final_xyz, "dmfb_strict_", ["--strict-path-method"])
    if strict.returncode != 2:
        print(strict.stdout)
        print(f"FAIL: neb_interpolator returned {strict.returncode} for strict fallback run (expected 2)")
        return 1

    require_contains(strict.stdout, "Status: GENERATED_WITH_FALLBACK", "strict fallback status")
    require_contains(strict.stdout, "Actual path method: DM", "strict actual path method")
    require_contains(strict.stdout, "Requested path method was not realized", "strict exit explanation")

    print("PASS: fallback reporting and strict-path-method semantics are correct")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

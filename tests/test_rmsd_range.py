#!/usr/bin/env python3
"""Regression test for calc_rmsd_xyz range arguments.

This catches a historical bug where the code reused the same variable for:
- argument count (iargc)
- index('.') when detecting file extension

As a result, 4-argument invocations could silently ignore the ranges.

The test constructs two 3-atom structures that differ only at atom 3.
With range=1-2, RMSD should be ~0.
Without range, RMSD must be noticeably > 0.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


def _run(cmd: list[str]) -> tuple[int, str]:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.returncode, p.stdout


def _extract_rmsd(output: str) -> float:
    m = re.search(r"RMSD\s*=\s*([0-9eE+\-.]+)", output)
    if not m:
        raise ValueError("RMSD line not found")
    return float(m.group(1))


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: test_rmsd_range.py <calc_rmsd_xyz> <ref.xyz> <mobile.xyz>")
        return 2

    exe = Path(sys.argv[1])
    ref = Path(sys.argv[2])
    mob = Path(sys.argv[3])

    if not exe.exists():
        print(f"FAIL: executable not found: {exe}")
        return 2
    if not ref.exists() or not mob.exists():
        print(f"FAIL: input file(s) missing: {ref} {mob}")
        return 2

    # 1) With range: should be ~0
    rc, out = _run([str(exe), str(ref), str(mob), "1-2", "1-2"])
    if rc != 0:
        print(out)
        print(f"FAIL: calc_rmsd_xyz returned {rc} (range mode)")
        return 1
    try:
        rmsd_range = _extract_rmsd(out)
    except Exception as e:
        print(out)
        print(f"FAIL: {e}")
        return 1

    # 2) Without range: should be clearly non-zero for this pair
    rc2, out2 = _run([str(exe), str(ref), str(mob)])
    if rc2 != 0:
        print(out2)
        print(f"FAIL: calc_rmsd_xyz returned {rc2} (no-range mode)")
        return 1
    try:
        rmsd_all = _extract_rmsd(out2)
    except Exception as e:
        print(out2)
        print(f"FAIL: {e}")
        return 1

    # Thresholds: keep loose to tolerate BLAS/LAPACK differences
    if rmsd_range > 1e-6:
        print(out)
        print(f"FAIL: RMSD with range too large: {rmsd_range}")
        return 1
    if rmsd_all < 1e-2:
        print(out2)
        print(f"FAIL: RMSD without range unexpectedly small: {rmsd_all}")
        return 1

    print(f"PASS: range RMSD={rmsd_range:.6e}, no-range RMSD={rmsd_all:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

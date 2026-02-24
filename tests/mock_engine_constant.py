#!/usr/bin/env python3
"""A tiny mock NEB external engine.

It reads the engine input file (XYZ blocks) and writes an output file with
matching blocks but constant vectors (here: all zeros, type=gradient).

Usage:
  python3 mock_engine_constant.py <infile> <outfile>

The neb_interpolator will typically call it like:
  --engine-cmd "python3 mock_engine_constant.py"
(and will append infile/outfile automatically).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def parse_blocks(lines: list[str]):
    blocks = []
    i = 0
    seq = 1
    while i < len(lines):
        # skip empty lines
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break
        nat = int(lines[i].strip())
        i += 1
        if i >= len(lines):
            raise ValueError("Unexpected EOF reading header")
        header = lines[i].strip()
        i += 1
        coords = lines[i : i + nat]
        if len(coords) != nat:
            raise ValueError("Unexpected EOF reading coordinates")
        i += nat

        m = re.search(r"\bimage\s*=\s*(\d+)", header)
        img = int(m.group(1)) if m else seq
        m = re.search(r"\bunits\s*=\s*([^\s]+)", header)
        units = m.group(1) if m else "Angstrom"
        blocks.append((nat, img, units, coords))
        seq += 1
    return blocks


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: mock_engine_constant.py <infile> <outfile>")
        return 2

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    lines = in_path.read_text().splitlines()
    blocks = parse_blocks(lines)

    out_lines: list[str] = []
    for nat, img, units, coords in blocks:
        out_lines.append(str(nat))
        out_lines.append(f"image={img} units={units} type=gradient")
        for line in coords:
            sym = line.split()[0]
            out_lines.append(f"{sym}  0.0  0.0  0.0")

    out_path.write_text("\n".join(out_lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

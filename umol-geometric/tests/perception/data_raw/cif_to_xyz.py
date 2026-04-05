#!/usr/bin/env python3
"""Convert small-molecule CIF files to XYZ by expanding the asymmetric unit
and extracting individual molecules via connectivity analysis.

Usage:
    python cif_to_xyz.py input.cif [output_dir]
    python cif_to_xyz.py *.cif -o output_dir

Handles:
- Non-CIF preamble lines (e.g., journal headers)
- Latin-1 encoded files (copyright symbols, etc.)
- Malformed CIF values (e.g., "0.5 mm" instead of "0.5")
- Molecules spanning periodic boundaries (minimum image unwrapping)
- Multiple molecules per unit cell (outputs the largest by default, or all with --all)
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import gemmi

COVALENT_RADII = {
    "H": 0.32, "He": 0.46, "Li": 1.33, "Be": 1.02, "B": 0.85, "C": 0.77,
    "N": 0.71, "O": 0.73, "F": 0.64, "Ne": 0.67, "Na": 1.55, "Mg": 1.39,
    "Al": 1.26, "Si": 1.16, "P": 1.07, "S": 1.05, "Cl": 0.99, "Ar": 0.96,
    "K": 1.96, "Ca": 1.71, "Sc": 1.48, "Ti": 1.36, "V": 1.34, "Cr": 1.22,
    "Mn": 1.19, "Fe": 1.16, "Co": 1.11, "Ni": 1.10, "Cu": 1.12, "Zn": 1.18,
    "Ga": 1.24, "Ge": 1.21, "As": 1.21, "Se": 1.16, "Br": 1.14, "Kr": 1.17,
    "Rb": 2.10, "Sr": 1.85, "Y": 1.63, "Zr": 1.54, "Nb": 1.47, "Mo": 1.38,
    "Tc": 1.28, "Ru": 1.25, "Rh": 1.25, "Pd": 1.20, "Ag": 1.28, "Cd": 1.36,
    "In": 1.42, "Sn": 1.40, "Sb": 1.40, "Te": 1.36, "I": 1.33, "Xe": 1.31,
}
BOND_TOLERANCE = 0.4


def read_cif(path: Path) -> gemmi.cif.Document:
    raw = path.read_bytes()
    text = raw.decode("latin-1")
    # Skip non-CIF preamble
    idx = text.find("data_")
    if idx < 0:
        raise ValueError(f"No data_ block found in {path}")
    text = text[idx:]
    # Fix common CIF formatting issues
    text = re.sub(r"(_exptl_crystal_size_\w+)\s+([\d.]+)\s+mm", r"\1  \2", text)
    # Quote unquoted multi-word values on data lines (not loop_ lines)
    def quote_multiword(line):
        m = re.match(r"^(_\S+)\s{2,}([^'\";#\n][^\n]*\S)\s*$", line)
        if m and " " in m.group(2).strip():
            return f"{m.group(1)}  '{m.group(2).strip()}'\n"
        return line
    text = "".join(quote_multiword(line) for line in text.splitlines(keepends=True))
    return gemmi.cif.read_string(text)


def expand_unit_cell(st):
    """Expand asymmetric unit by space group operations. Returns [(element, fx, fy, fz)]."""
    ops = st.spacegroup.operations()
    atoms = []
    for site in st.sites:
        seen = set()
        for op in ops:
            fract = op.apply_to_xyz([site.fract.x, site.fract.y, site.fract.z])
            fract = [x % 1.0 for x in fract]
            key = tuple(round(x, 4) for x in fract)
            if key not in seen:
                seen.add(key)
                atoms.append((site.element.name, fract[0], fract[1], fract[2]))
    return atoms


def min_image_dist(f1, f2, cell):
    df = [f1[i] - f2[i] for i in range(3)]
    df = [d - round(d) for d in df]
    cart = cell.orthogonalize(gemmi.Fractional(*df))
    return (cart.x ** 2 + cart.y ** 2 + cart.z ** 2) ** 0.5


def max_bond_dist(e1, e2):
    return COVALENT_RADII.get(e1, 1.5) + COVALENT_RADII.get(e2, 1.5) + BOND_TOLERANCE


def extract_molecules(atoms, cell):
    """Find connected molecules and unwrap coordinates across periodic boundaries."""
    n = len(atoms)

    # Build adjacency
    adj = defaultdict(set)
    for i in range(n):
        for j in range(i + 1, n):
            d = min_image_dist(atoms[i][1:4], atoms[j][1:4], cell)
            if 0.1 < d < max_bond_dist(atoms[i][0], atoms[j][0]):
                adj[i].add(j)
                adj[j].add(i)

    # BFS with coordinate unwrapping
    all_visited = set()
    molecules = []
    for start in range(n):
        if start in all_visited:
            continue
        visited = {start}
        unwrapped = {start: list(atoms[start][1:4])}
        queue = [start]
        while queue:
            cur = queue.pop(0)
            fc = unwrapped[cur]
            for nb in adj[cur]:
                if nb in visited:
                    continue
                visited.add(nb)
                fn = list(atoms[nb][1:4])
                for k in range(3):
                    fn[k] -= round(fn[k] - fc[k])
                unwrapped[nb] = fn
                queue.append(nb)

        mol = []
        for idx in visited:
            f = unwrapped[idx]
            pos = cell.orthogonalize(gemmi.Fractional(*f))
            mol.append((atoms[idx][0], pos.x, pos.y, pos.z))
        all_visited |= visited
        molecules.append(mol)

    # Sort by size descending
    molecules.sort(key=len, reverse=True)
    return molecules


def formula(mol):
    ec = defaultdict(int)
    for e, *_ in mol:
        ec[e] += 1
    return "".join(f"{e}{c}" for e, c in sorted(ec.items()))


def write_xyz(mol, path, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(mol)}\n")
        f.write(f"{comment}\n")
        for e, x, y, z in mol:
            f.write(f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def process_cif(cif_path: Path, output_dir: Path, all_molecules=False):
    doc = read_cif(cif_path)
    results = []

    for block in doc:
        try:
            st = gemmi.make_small_structure_from_block(block)
        except Exception as e:
            print(f"  Warning: skipping block '{block.name}': {e}", file=sys.stderr)
            continue

        atoms = expand_unit_cell(st)
        molecules = extract_molecules(atoms, st.cell)

        if not molecules:
            print(f"  Warning: no molecules found in block '{block.name}'", file=sys.stderr)
            continue

        stem = cif_path.stem
        if all_molecules:
            for i, mol in enumerate(molecules):
                name = f"{stem}_mol{i}" if len(molecules) > 1 else stem
                out = output_dir / f"{name}.xyz"
                comment = f"{formula(mol)} from {cif_path.name} block={block.name}"
                write_xyz(mol, out, comment)
                results.append((out, len(mol), formula(mol)))
        else:
            mol = molecules[0]
            out = output_dir / f"{stem}.xyz"
            comment = f"{formula(mol)} from {cif_path.name} block={block.name}"
            write_xyz(mol, out, comment)
            results.append((out, len(mol), formula(mol)))

    return results


def main():
    parser = argparse.ArgumentParser(description="Convert CIF to XYZ")
    parser.add_argument("cif_files", nargs="+", type=Path)
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("."))
    parser.add_argument("--all", action="store_true", help="Output all molecules, not just the largest")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for cif in args.cif_files:
        print(f"Processing {cif}...")
        results = process_cif(cif, args.output_dir, args.all)
        for path, natoms, form in results:
            print(f"  {path}: {natoms} atoms, {form}")


if __name__ == "__main__":
    main()

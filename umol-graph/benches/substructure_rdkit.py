"""RDKit VF2 substructure matching baselines."""
import time, pathlib, glob
from rdkit import Chem, RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)

smiles_files = glob.glob(str(pathlib.Path(__file__).resolve().parent.parent / "tests/smiles_parsing/data/basic_opensmiles/**/*.smiles"), recursive=True)
smiles_list = []
for f in sorted(smiles_files):
    lines = open(f).readlines()
    if len(lines) >= 2 and lines[1].strip():
        smiles_list.append(lines[1].strip())

mols = [Chem.MolFromSmiles(s) for s in smiles_list]
mols = [m for m in mols if m is not None]
print(f"{len(mols)} molecules")

patterns = {
    "branched":  Chem.MolFromSmarts("[#6](~[#6])~[#6](~[#6])~[#7]"),
    "phenol":    Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8]"),
    "bicyclic":  Chem.MolFromSmarts("[#6]1~[#6](~[#6]~[#6]~[#6]2)~[#6]2~[#6]~[#6]~[#6]1"),
}

N = 10

for name, pat in patterns.items():
    match_count = sum(1 for m in mols if m.HasSubstructMatch(pat))

    # Warmup
    for m in mols:
        m.GetSubstructMatches(pat)

    t0 = time.perf_counter()
    for _ in range(N):
        for m in mols:
            m.GetSubstructMatches(pat)
    elapsed = (time.perf_counter() - t0) / N
    print(f"{name:10s}: {elapsed*1000:7.1f} ms / {len(mols)} mols ({match_count} hits) = {elapsed/len(mols)*1e6:.1f} us/mol")

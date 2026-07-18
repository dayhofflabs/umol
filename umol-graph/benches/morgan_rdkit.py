"""RDKit SMILES parsing and Morgan fingerprint baselines."""
import time, pathlib, glob
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.logger().setLevel(RDLogger.ERROR)

smiles_files = glob.glob(str(pathlib.Path(__file__).resolve().parents[2] / "umol-io/tests/smiles_parsing/data/opensmiles/**/*.smiles"), recursive=True)
smiles_list = []
for f in sorted(smiles_files):
    lines = open(f).readlines()
    if len(lines) >= 2 and lines[1].strip():
        smiles_list.append(lines[1].strip())

print(f"{len(smiles_list)} SMILES strings")

# --- SMILES parsing (no sanitization) ---
for s in smiles_list:
    Chem.MolFromSmiles(s, sanitize=False)

N = 10
t0 = time.perf_counter()
for _ in range(N):
    for s in smiles_list:
        Chem.MolFromSmiles(s, sanitize=False)
elapsed = (time.perf_counter() - t0) / N
parsed_nosanit = sum(1 for s in smiles_list if Chem.MolFromSmiles(s, sanitize=False) is not None)
print(f"SMILES parse (no sanitize): {elapsed*1000:.1f} ms / {parsed_nosanit} mols = {elapsed/parsed_nosanit*1e6:.1f} us/mol")

# --- SMILES parsing (with sanitization) ---
t0 = time.perf_counter()
for _ in range(N):
    for s in smiles_list:
        Chem.MolFromSmiles(s)
elapsed = (time.perf_counter() - t0) / N
mols = [Chem.MolFromSmiles(s) for s in smiles_list]
mols = [m for m in mols if m is not None]
print(f"SMILES parse (sanitized):   {elapsed*1000:.1f} ms / {len(mols)} mols = {elapsed/len(mols)*1e6:.1f} us/mol")

# --- Morgan fingerprint ---
for m in mols:
    AllChem.GetMorganFingerprint(m, radius=2)

t0 = time.perf_counter()
for _ in range(N):
    for m in mols:
        AllChem.GetMorganFingerprint(m, radius=2)
elapsed = (time.perf_counter() - t0) / N
print(f"ECFP4 fingerprint only:     {elapsed*1000:.1f} ms / {len(mols)} mols = {elapsed/len(mols)*1e6:.1f} us/mol")

# --- ECFP6 (radius 3) ---
for m in mols:
    AllChem.GetMorganFingerprint(m, radius=3)

t0 = time.perf_counter()
for _ in range(N):
    for m in mols:
        AllChem.GetMorganFingerprint(m, radius=3)
elapsed = (time.perf_counter() - t0) / N
print(f"ECFP6 fingerprint only:     {elapsed*1000:.1f} ms / {len(mols)} mols = {elapsed/len(mols)*1e6:.1f} us/mol")

# MOL Compliance Suite

## Data sources for Mol V2000 files for compliance suite

- [RDKit](https://github.com/rdkit/rdkit)
- [CDK](https://github.com/cdk/cdk)
- [Indigo](https://github.com/epam/Indigo/tree/master)
- [Ketcher](https://github.com/epam/ketcher)
- [OpenBabel](https://github.com/openbabel/openbabel)
- [ChEBI](https://www.ebi.ac.uk/chebi/init.do)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov)
- [KEGG](https://www.genome.jp/kegg/kegg2.html)
- [SciFinder](https://scifinder-n.cas.org/)
- [ChEMBL](https://www.ebi.ac.uk/chembl/)
- [ChemSpider](https://www.chemspider.com/)
- [Reaxys](https://www.reaxys.com/)
- [Jmol/JSmol](https://sourceforge.net/projects/jmol/)
- [ChemAxon/MarvinSketch](https://freetrial.marvin.cxn.io)
- [NIST WebBook](https://webbook.nist.gov/cgi/cbook.cgi)

Missing:

- ACD/ChemSketch
- ChemBioDraw
- ChemDoodle

Additionally retrieved MOL3k, SDF, RXN files, but not touching them for now.

## MOL File Classification Results

(Using umol-models-graphs/bin/mol_classifier)

Classification based on basic parser success:

- **Basic**: Files that parse successfully with the molecule parser
- **Extended**: Files that require the full (moleculelike) parser
  (queries, S-groups, etc.)

  | Source     | Total    | Molecule | MoleculeLike | Invalid | Valid %   |
  | ---------- | -------- | -------- | ------------ | ------- | --------- |
  | cdk        | 219      | 164      | 28           | 27      | 87.7%     |
  | chebi      | 16       | 7        | 9            | 0       | 100.0%    |
  | chembl     | 2        | 0        | 0            | 2       | 0.0%      |
  | chemspider | 10       | 9        | 0            | 1       | 90.0%     |
  | indigo     | 448      | 312      | 49           | 87      | 80.6%     |
  | jmol       | 112      | 112      | 0            | 0       | 100.0%    |
  | ketcher    | 269      | 147      | 52           | 70      | 74.0%     |
  | marvin     | 24       | 14       | 5            | 5       | 79.2%     |
  | nist       | 9        | 9        | 0            | 0       | 100.0%    |
  | openbabel  | 219      | 216      | 1            | 2       | 99.1%     |
  | rdkit      | 387      | 297      | 67           | 23      | 94.1%     |
  | reaxys     | 1        | 1        | 0            | 0       | 100.0%    |
  | rhea       | 11       | 3        | 7            | 1       | 90.9%     |
  | scifinder  | 14       | 12       | 2            | 0       | 100.0%    |
  | **Total**  | **1741** | **1303** | **220**      | **218** | **87.5%** |

**Summary:**

- Processed: 1741/1741 files
- 74.8% of files are basic (can use optimized parser)

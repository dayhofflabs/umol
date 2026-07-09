from umol import MoleculeAst


def test_molecule_ast_empty():
    assert len(MoleculeAst().atoms) == 0


def test_molecule_ast_eq():
    assert MoleculeAst() == MoleculeAst()


def test_molecule_ast_repr():
    assert repr(MoleculeAst()) == "MoleculeAst(atoms=0)"

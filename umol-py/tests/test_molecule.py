import pytest

from umol import AtomAst, BondAst, Element, MoleculeAst


def test_molecule_ast_empty():
    assert len(MoleculeAst().atoms) == 0
    assert len(MoleculeAst().bonds) == 0


def test_molecule_ast_eq():
    assert MoleculeAst() == MoleculeAst()


def test_molecule_ast_repr():
    assert repr(MoleculeAst()) == "MoleculeAst(atoms=0, bonds=0)"


def test_molecule_ast_from_atoms_and_bonds():
    mol = MoleculeAst.from_atoms_and_bonds(
        [AtomAst(Element("C")), AtomAst(Element("C"))],
        [(0, 1, BondAst(2))],
    )
    assert len(mol.atoms) == 2
    assert len(mol.bonds) == 1
    assert repr(mol) == "MoleculeAst(atoms=2, bonds=1)"


def test_molecule_ast_from_atoms_and_bonds_default_bonds():
    mol = MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))])
    assert len(mol.atoms) == 1
    assert len(mol.bonds) == 0


def test_molecule_ast_bonds_out_of_range():
    with pytest.raises(IndexError):
        MoleculeAst.from_atoms_and_bonds([AtomAst(Element("C"))]).bonds[0]

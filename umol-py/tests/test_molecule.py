import pytest

from umol import (
    AromaticSystemAst,
    AtomAst,
    BondAst,
    Element,
    MoleculeAst,
    NoncovalentBondAst,
    NoncovalentBondKind,
)


def test_molecule_ast_empty():
    assert len(MoleculeAst().atoms) == 0
    assert len(MoleculeAst().bonds) == 0


def test_molecule_ast_eq():
    assert MoleculeAst() == MoleculeAst()


def test_molecule_ast_repr():
    assert repr(MoleculeAst()) == "MoleculeAst(atoms=0, bonds=0)"


def test_molecule_ast_repr_includes_all_entities():
    # atoms + bonds always; the other families surface only when present
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("C")), AtomAst(Element("C")), AtomAst(Element("O"))],
        bonds=[(0, 1, BondAst(1))],
        aromatic=[([0, 1, 2], AromaticSystemAst([1, 1, 1]))],
        noncovalent=[([0, 2], NoncovalentBondAst(NoncovalentBondKind.HydrogenBond))],
    )
    assert repr(mol) == "MoleculeAst(atoms=3, bonds=1, aromatic=1, noncovalent=1)"


def test_molecule_ast_from_parts():
    mol = MoleculeAst.from_parts(
        [AtomAst(Element("C")), AtomAst(Element("C"))],
        bonds=[(0, 1, BondAst(2))],
    )
    assert len(mol.atoms) == 2
    assert len(mol.bonds) == 1
    assert repr(mol) == "MoleculeAst(atoms=2, bonds=1)"


def test_molecule_ast_from_parts_default_bonds():
    mol = MoleculeAst.from_parts([AtomAst(Element("C"))])
    assert len(mol.atoms) == 1
    assert len(mol.bonds) == 0


def test_molecule_ast_bonds_out_of_range():
    with pytest.raises(IndexError):
        MoleculeAst.from_parts([AtomAst(Element("C"))]).bonds[0]

import pytest

from umol import (
    AtomAst,
    CisTransStereo,
    CisTransStereoAst,
    Element,
    MoleculeAst,
    Permutation,
    StereoAtomAst,
    StereoAtomConstraintAst,
    StereoAtomConstraintsAst,
    StereoBondAst,
    StereoConfigurationAst,
    StereoCosetAst,
    StereoKind,
    StereoLigand,
    StereoLigandKind,
    StereoTerm,
    Stereogenicity,
    StereogenicityAst,
    TetrahedralStereo,
    TetrahedralStereoAst,
)


def test_permutation_image_degree():
    p = Permutation([1, 0, 2, 3])
    assert p.image() == [1, 0, 2, 3]
    assert p.degree == 4


@pytest.mark.parametrize(
    ("image", "message"),
    [
        ([0, 1, 2, 3, 4, 5, 6], "permutation image length 7 exceeds maximum 6"),
        ([0, 1, 3], "permutation image value 3 at position 2 is outside 0..3"),
        ([0, 1, 1], "permutation image value 1 occurs more than once"),
    ],
)
def test_permutation_error(image, message):
    with pytest.raises(ValueError) as error:
        Permutation(image)
    assert str(error.value) == message


def test_permutation_identity():
    assert Permutation.identity(4).image() == [0, 1, 2, 3]


def test_permutation_eq_hash():
    assert Permutation([1, 0, 2]) == Permutation([1, 0, 2])
    assert Permutation([1, 0, 2]) != Permutation([0, 1, 2])
    assert len({Permutation([1, 0, 2]), Permutation([1, 0, 2])}) == 1


def test_stereoterm_lit():
    assert StereoTerm.Lit(1) == StereoTerm.Lit(1)


def test_stereoterm_apply():
    term = StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))
    assert term == StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))


def test_stereocosetast_term():
    assert StereoCosetAst.Term(StereoTerm.Lit(2)) == StereoCosetAst.Term(StereoTerm.Lit(2))


def test_tetrahedralstereoast_stereo():
    assert TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1)) == TetrahedralStereoAst.Stereo(
        StereoCosetAst.Lit(1)
    )


def test_tetrahedralstereoast_not_stereo():
    assert TetrahedralStereoAst.NotStereo() == TetrahedralStereoAst.NotStereo()


def test_stereocosetast_eq_hash_repr():
    assert StereoCosetAst.Lit(1) == StereoCosetAst.Lit(1)
    assert StereoCosetAst.Lit(1) != StereoCosetAst.Lit(0)
    assert len({StereoCosetAst.Lit(1), StereoCosetAst.Lit(1)}) == 1
    assert repr(StereoCosetAst.Lit(1)) == "StereoCosetAst.Lit(1)"


def test_stereoterm_eq_repr():
    assert StereoTerm.Lit(0) == StereoTerm.Lit(0)
    assert StereoTerm.Lit(0) != StereoTerm.Swap(StereoTerm.Lit(0))
    assert repr(StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))) == (
        "StereoTerm.Apply(StereoTerm.Lit(0), Permutation([1, 0, 2, 3]))"
    )


def test_tetrahedralstereoast_stereo_repr():
    assert repr(TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1))) == (
        "TetrahedralStereoAst.Stereo(StereoCosetAst.Lit(1))"
    )


def test_cistransstereoast_stereo():
    assert CisTransStereoAst.Stereo(StereoCosetAst.Lit(1)) == CisTransStereoAst.Stereo(
        StereoCosetAst.Lit(1)
    )


def test_cistransstereoast_not_stereo():
    assert CisTransStereoAst.NotStereo() == CisTransStereoAst.NotStereo()


def test_cistransstereoast_stereo_repr():
    assert repr(CisTransStereoAst.Stereo(StereoCosetAst.Lit(1))) == (
        "CisTransStereoAst.Stereo(StereoCosetAst.Lit(1))"
    )


def test_cistransstereo_enum():
    assert CisTransStereo.Z == CisTransStereo.Z
    assert CisTransStereo.Z != CisTransStereo.E
    assert len({CisTransStereo.E, CisTransStereo.E, CisTransStereo.Z}) == 2


def stereo_atom_molecule():
    # a tetrahedral stereocenter on atom 0 with four atom ligands (atoms 1-4)
    return MoleculeAst.from_parts(
        [AtomAst(Element("C")) for _ in range(5)],
        stereo_atoms=[
            (
                0,
                [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)],
                StereoAtomAst(TetrahedralStereo.Ccw),
            )
        ],
    )


def test_stereoatomast_new():
    atom = StereoAtomAst(TetrahedralStereo.Ccw)
    assert str(atom) == "Th0"
    assert atom.configuration == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
    )
    assert len(atom.constraints) == 0


@pytest.mark.parametrize("dsl", ["Th0", "Th*", "Sp2"])
def test_stereoatomast_parse_roundtrip(dsl):
    atom = StereoAtomAst.parse(dsl)
    assert str(atom) == dsl
    assert repr(atom) == f"StereoAtomAst.parse('{dsl}')"


def test_stereoatomast_configuration_setter():
    atom = StereoAtomAst(TetrahedralStereo.Ccw)
    atom.configuration = TetrahedralStereo.Cw
    assert str(atom) == "Th1"


def test_stereoatomast_constraints_kwarg():
    atom = StereoAtomAst(
        TetrahedralStereo.Ccw,
        constraints=StereoAtomConstraintsAst(
            [
                StereoAtomConstraintAst.Stereogenicity(
                    StereogenicityAst.Lit(Stereogenicity.Stereogenic)
                )
            ]
        ),
    )
    assert atom.constraints.stereogenicity() == StereogenicityAst.Lit(
        Stereogenicity.Stereogenic
    )


def test_stereoatomast_asdict():
    atom = StereoAtomAst(TetrahedralStereo.Ccw)
    d = atom.asdict()
    assert set(d.keys()) == {"configuration", "constraints"}
    assert d["configuration"] == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
    )
    assert d["constraints"] == []


def test_stereobondast_parse_roundtrip():
    bond = StereoBondAst.parse("Ct0")
    assert str(bond) == "Ct0"
    assert bond.configuration == StereoConfigurationAst.Kinded(
        StereoKind.CisTrans, StereoCosetAst.Lit(0)
    )


def test_molecule_stereo_atoms_from_parts():
    views = stereo_atom_molecule().stereo_atoms
    assert len(views) == 1
    view = views[0]
    assert view.id == 0
    assert view.site_id == 0
    assert view.ligands == [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)]
    assert view.kind == StereoKind.Tetrahedral
    assert view.coset == StereoCosetAst.Lit(0)
    assert view.configuration == StereoConfigurationAst.Kinded(
        StereoKind.Tetrahedral, StereoCosetAst.Lit(0)
    )


def test_stereoatomview_configuration_write_through():
    mol = stereo_atom_molecule()
    mol.stereo_atoms[0].configuration = TetrahedralStereo.Cw
    # a fresh view re-reads the molecule, proving the write landed on it
    assert mol.stereo_atoms[0].coset == StereoCosetAst.Lit(1)


def test_stereoatomview_constraints_write_through():
    mol = stereo_atom_molecule()
    mol.stereo_atoms[0].constraints.set(
        StereoAtomConstraintAst.Stereogenicity(
            StereogenicityAst.Lit(Stereogenicity.Stereogenic)
        )
    )
    assert mol.stereo_atoms[0].constraints.stereogenicity() == StereogenicityAst.Lit(
        Stereogenicity.Stereogenic
    )


def test_stereoatomviews_at():
    views = stereo_atom_molecule().stereo_atoms
    assert views.at(0).id == 0
    assert views.at(1) is None


def test_stereoatomviews_of():
    views = stereo_atom_molecule().stereo_atoms
    ligands = [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)]
    # order-independent full-ligand-set match
    assert views.of(0, list(reversed(ligands))).id == 0
    # a partial ligand set does not match
    assert views.of(0, ligands[:2]) is None


def test_stereoatomviews_setitem():
    mol = stereo_atom_molecule()
    mol.stereo_atoms[0] = StereoAtomAst.parse("Th1")
    view = mol.stereo_atoms[0]
    assert view.coset == StereoCosetAst.Lit(1)
    # site topology preserved
    assert view.site_id == 0


def test_molecule_repr_stereo():
    assert repr(stereo_atom_molecule()) == (
        "MoleculeAst(atoms=5, bonds=0, stereo_atoms=1)"
    )

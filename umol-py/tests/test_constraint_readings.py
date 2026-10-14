import sys

import pytest

from umol import (
    AromaticSystemConstraintForm,
    AromaticSystemConstraintsForm,
    AromaticSystemForm,
    AromaticValenceForm,
    AtomConstraintForm,
    AtomConstraintKey,
    AtomConstraintsForm,
    AtomForm,
    BondConstraintForm,
    BondConstraintKey,
    BondConstraintsForm,
    BondForm,
    BooleanForm,
    CisTransConfiguration,
    CisTransStereoForm,
    DativeBondConstraintForm,
    DativeBondConstraintKey,
    DativeBondConstraintsForm,
    DativeBondForm,
    Edits,
    Element,
    ElementForm,
    FluxionalityForm,
    InvalidatedViewError,
    LigandPermutation,
    LigandSymmetryForm,
    Molecule,
    MulticenterBondConstraintForm,
    MulticenterBondConstraintsForm,
    MulticenterBondForm,
    NoncovalentBondConstraintForm,
    NoncovalentBondConstraintsForm,
    NoncovalentBondForm,
    NoncovalentBondKind,
    NumForm,
    Orientation,
    OrientedLigandPermutation,
    Permutation,
    RingConfig,
    RingMembershipForm,
    RingModel,
    RingScope,
    StereoAtomConstraintForm,
    StereoAtomConstraintsForm,
    StereoAtomForm,
    StereoBondConstraintForm,
    StereoBondConstraintsForm,
    StereoBondForm,
    StereoCoset,
    StereoLigand,
    StereoLigandKind,
    StereoLigandPair,
    Stereogenicity,
    StereogenicityForm,
    TetrahedralConfiguration,
    TetrahedralStereoForm,
    Topicity,
    TopicityForm,
    TopicityRelationForm,
)


READINGS = (
    "asserted", "asserted_complete", "derived", "derived_complete",
    "evidenced", "evidenced_complete",
)

ENTITY_CASES = [
    pytest.param(
        "atoms", AtomConstraintsForm,
        AtomConstraintForm.Valence(NumForm.LitSet({5})),
        AtomConstraintForm.Valence(NumForm.Lit(5)),
        AtomConstraintForm.Valence(NumForm.Undetermined()),
        AtomConstraintForm.Valence(NumForm.Lit(5)),
        AtomConstraintForm.Valence(NumForm.Undetermined()),
        id="atom",
    ),
    pytest.param(
        "bonds", BondConstraintsForm,
        BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        BondConstraintForm.Aromatic(BooleanForm.Undetermined()),
        BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        BondConstraintForm.Aromatic(BooleanForm.Lit(False)),
        id="bond",
    ),
    pytest.param(
        "dative_bonds", DativeBondConstraintsForm,
        DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        DativeBondConstraintForm.Aromatic(BooleanForm.Undetermined()),
        DativeBondConstraintForm.Aromatic(BooleanForm.Lit(True)),
        DativeBondConstraintForm.Aromatic(BooleanForm.Lit(False)),
        id="dative",
    ),
    pytest.param(
        "aromatic_systems", AromaticSystemConstraintsForm,
        AromaticSystemConstraintForm.ElectronCount(NumForm.LitSet({3})),
        AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(3)),
        AromaticSystemConstraintForm.ElectronCount(NumForm.Undetermined()),
        AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(3)),
        AromaticSystemConstraintForm.ElectronCount(NumForm.Undetermined()),
        id="aromatic",
    ),
    pytest.param(
        "multicenter_bonds", MulticenterBondConstraintsForm,
        MulticenterBondConstraintForm.ElectronCount(NumForm.LitSet({3})),
        MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(3)),
        MulticenterBondConstraintForm.ElectronCount(NumForm.Undetermined()),
        MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(3)),
        MulticenterBondConstraintForm.ElectronCount(NumForm.Undetermined()),
        id="multicenter",
    ),
    pytest.param(
        "noncovalent_bonds", NoncovalentBondConstraintsForm,
        NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True)),
        NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True)),
        NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Undetermined()),
        NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(True)),
        NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Undetermined()),
        id="noncovalent",
    ),
    pytest.param(
        "stereo_atoms", StereoAtomConstraintsForm,
        StereoAtomConstraintForm.Stereogenicity(
            StereogenicityForm.LitSet({Stereogenicity.Stereogenic})
        ),
        StereoAtomConstraintForm.Stereogenicity(
            StereogenicityForm.Lit(Stereogenicity.Stereogenic)
        ),
        StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        id="stereo_atom",
    ),
    pytest.param(
        "stereo_bonds", StereoBondConstraintsForm,
        StereoBondConstraintForm.Stereogenicity(
            StereogenicityForm.LitSet({Stereogenicity.Symmetric})
        ),
        StereoBondConstraintForm.Stereogenicity(
            StereogenicityForm.Lit(Stereogenicity.Symmetric)
        ),
        StereoBondConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        StereoBondConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        StereoBondConstraintForm.Stereogenicity(StereogenicityForm.Undetermined()),
        id="stereo_bond",
    ),
]


@pytest.fixture
def molecule():
    return Molecule.from_entries(
        [AtomForm(Element("C"), implicit_hydrogens=0) for _ in range(5)],
        bonds=[
            (0, 1, BondForm(2)), (0, 2, BondForm(1)),
            (0, 3, BondForm(1)), (0, 4, BondForm(1)), (1, 3, BondForm(1)),
        ],
        dative_bonds=[([2], 1, DativeBondForm(1))],
        aromatic_systems=[([0, 1, 2], AromaticSystemForm([1, 1, 1]))],
        multicenter_bonds=[([0, 1, 2], MulticenterBondForm([1, 1, 1]))],
        noncovalent_bonds=[
            ([0, 2], NoncovalentBondForm(NoncovalentBondKind.HydrogenBond))
        ],
        stereo_atoms=[
            (0, [StereoLigand(i, StereoLigandKind.Atom) for i in range(1, 5)],
             StereoAtomForm(TetrahedralConfiguration.Ccw))
        ],
        stereo_bonds=[
            (0, [StereoLigand(2, StereoLigandKind.Atom),
                 StereoLigand(0, StereoLigandKind.ImplicitHydrogen),
                 StereoLigand(3, StereoLigandKind.Atom),
                 StereoLigand(1, StereoLigandKind.ImplicitHydrogen)],
             StereoBondForm(CisTransConfiguration.Z))
        ],
    )


@pytest.fixture
def forms():
    return {
        "atoms": AtomForm(Element("C")),
        "bonds": BondForm(1),
        "dative_bonds": DativeBondForm(1),
        "aromatic_systems": AromaticSystemForm([1, 1, 1]),
        "multicenter_bonds": MulticenterBondForm([1, 1, 1]),
        "noncovalent_bonds": NoncovalentBondForm(NoncovalentBondKind.HydrogenBond),
        "stereo_atoms": StereoAtomForm(TetrahedralConfiguration.Ccw),
        "stereo_bonds": StereoBondForm(CisTransConfiguration.Z),
    }


@pytest.mark.parametrize(
    "namespace,container_type,raw,expected,top,derived,closed", ENTITY_CASES
)
def test_entity_constraints_form_asserted(
    namespace, container_type, raw, expected, top, derived, closed,
):
    constraints = container_type([])
    key = raw.key
    assert constraints.asserted(key) == top
    assert constraints.asserted_complete(key) == closed
    assert constraints.get(key) is None
    constraints.set(raw)
    assert constraints.asserted(key) == expected
    assert constraints.asserted_complete(key) == expected
    assert constraints.get(key) == raw
    saved = constraints.asserted(key)
    constraints.pop(key)
    assert saved == expected
    assert saved.normalize() == expected


@pytest.mark.parametrize(
    "namespace,container_type,raw,expected,top,derived,closed", ENTITY_CASES
)
@pytest.mark.parametrize("backing", ["standalone", "molecule"])
def test_entity_constraints_view_asserted(
    molecule, forms, namespace, container_type, raw, expected, top, derived, closed, backing,
):
    entity = forms[namespace] if backing == "standalone" else getattr(molecule, namespace)[0]
    view = entity.constraints
    key = raw.key
    assert view.asserted(key) == top
    assert view.asserted_complete(key) == closed
    entity.constraints = container_type([raw])
    assert view.asserted(key) == expected
    assert view.asserted_complete(key) == expected
    assert view.get(key) == raw
    saved = view.asserted(key)
    view.pop(key)
    assert view.asserted(key) == top
    assert saved == expected


@pytest.mark.parametrize(
    "namespace,container_type,raw,expected,top,derived,closed", ENTITY_CASES
)
def test_entity_constraints_view_evidenced(
    molecule, namespace, container_type, raw, expected, top, derived, closed,
):
    entity = getattr(molecule, namespace)[0]
    view = entity.constraints
    key = raw.key
    for method in READINGS[2:]:
        assert getattr(view, method)(key) == derived
    entity.constraints = container_type([raw])
    assert view.derived(key) == derived
    assert view.derived_complete(key) == derived
    assert view.evidenced(key) == expected
    assert view.evidenced_complete(key) == expected
    assert view.get(key) == raw
    assert view.satisfies(container_type([expected]))
    assert view.is_compatible(container_type([expected]))
    assert view.satisfies(container_type([]))
    assert view.is_compatible(container_type([]))
    saved = view.evidenced(key)
    view.pop(key)
    assert view.evidenced(key) == derived
    assert saved == expected


@pytest.mark.parametrize(
    "namespace,container_type,raw,expected,top,derived,closed", ENTITY_CASES
)
def test_entity_constraints_view_derived_error(
    forms, namespace, container_type, raw, expected, top, derived, closed,
):
    view = forms[namespace].constraints
    for method in READINGS[2:]:
        with pytest.raises(ValueError, match="requires molecule context"):
            getattr(view, method)(raw.key)
    for method in ("satisfies", "is_compatible"):
        for entries in ([], [expected]):
            with pytest.raises(ValueError, match="requires molecule context"):
                getattr(view, method)(container_type(entries))


@pytest.mark.parametrize(
    "namespace,container_type,raw,expected,top,derived,closed", ENTITY_CASES
)
@pytest.mark.parametrize("operation", ["transact", "tracked_transact", "combine"])
def test_entity_constraints_view_evidenced_invalidation(
    molecule, namespace, container_type, raw, expected, top, derived, closed, operation,
):
    view = getattr(molecule, namespace)[0].constraints
    view.set(raw)
    saved = view.evidenced(raw.key)
    if operation == "combine":
        molecule.combine_from(Molecule())
    else:
        getattr(molecule, operation)([Edits()])
    for method in READINGS:
        with pytest.raises(InvalidatedViewError):
            getattr(view, method)(raw.key)
    for method in ("satisfies", "is_compatible"):
        with pytest.raises(InvalidatedViewError):
            getattr(view, method)(container_type([]))
    assert saved == expected


@pytest.mark.parametrize(
    "namespace,container_type,raw,expected,top,derived,closed", ENTITY_CASES
)
@pytest.mark.parametrize("method", ["edit", "apply", "tracked_apply"])
def test_entity_constraints_view_evidenced_edit(
    molecule, namespace, container_type, raw, expected, top, derived, closed, method,
):
    view = getattr(molecule, namespace)[0].constraints
    view.set(raw)
    before = [getattr(view, reading)(raw.key) for reading in READINGS]
    if method == "edit":
        assert molecule.edit().finish() == molecule
    else:
        getattr(molecule, method)(Edits())
    assert [getattr(view, reading)(raw.key) for reading in READINGS] == before
    for reading in ("satisfies", "is_compatible"):
        assert getattr(view, reading)(container_type([]))
    view.set(raw)
    assert getattr(molecule, namespace)[0].constraints.get(raw.key) == raw


@pytest.mark.parametrize(
    "namespace,container_type,bottom",
    [
        ("atoms", AtomConstraintsForm, AtomConstraintForm.Valence(NumForm.LitSet(set()))),
        ("bonds", BondConstraintsForm, BondConstraintForm.RingMembership(
            RingMembershipForm(RingScope.All(), NumForm.LitSet(set())))),
        ("dative_bonds", DativeBondConstraintsForm, DativeBondConstraintForm.RingMembership(
            RingMembershipForm(RingScope.All(), NumForm.LitSet(set())))),
        ("aromatic_systems", AromaticSystemConstraintsForm,
         AromaticSystemConstraintForm.ElectronCount(NumForm.LitSet(set()))),
        ("multicenter_bonds", MulticenterBondConstraintsForm,
         MulticenterBondConstraintForm.ElectronCount(NumForm.LitSet(set()))),
        ("stereo_atoms", StereoAtomConstraintsForm,
         StereoAtomConstraintForm.Stereogenicity(StereogenicityForm.LitSet(set()))),
        ("stereo_bonds", StereoBondConstraintsForm,
         StereoBondConstraintForm.Stereogenicity(StereogenicityForm.LitSet(set()))),
    ],
)
def test_entity_constraints_view_asserted_bottom(molecule, forms, namespace, container_type, bottom):
    constraints = container_type([bottom])
    forms[namespace].constraints = constraints
    entity = getattr(molecule, namespace)[0]
    entity.constraints = constraints
    view = entity.constraints
    if namespace in ("atoms", "bonds"):
        rings = molecule.rings(RingModel(), RingConfig()).into_ring_set()
        view = view.with_rings(rings)
    for source in (constraints, forms[namespace].constraints, view):
        assert source.asserted(bottom.key) is None
        assert source.asserted_complete(bottom.key) is None
        assert source.get(bottom.key) == bottom
    assert view.evidenced(bottom.key) is None
    assert view.evidenced_complete(bottom.key) is None
    assert not view.satisfies(constraints)
    assert not view.is_compatible(constraints)


@pytest.mark.parametrize(
    "namespace,container_type,assertion",
    [
        ("atoms", AtomConstraintsForm, AtomConstraintForm.Valence(NumForm.Lit(4))),
        ("bonds", BondConstraintsForm, BondConstraintForm.Aromatic(BooleanForm.Lit(False))),
        ("dative_bonds", DativeBondConstraintsForm,
         DativeBondConstraintForm.Aromatic(BooleanForm.Lit(False))),
        ("aromatic_systems", AromaticSystemConstraintsForm,
         AromaticSystemConstraintForm.ElectronCount(NumForm.Lit(4))),
        ("multicenter_bonds", MulticenterBondConstraintsForm,
         MulticenterBondConstraintForm.ElectronCount(NumForm.Lit(4))),
        ("noncovalent_bonds", NoncovalentBondConstraintsForm,
         NoncovalentBondConstraintForm.Intramolecular(BooleanForm.Lit(False))),
    ],
)
def test_entity_constraints_view_evidenced_conflict(molecule, namespace, container_type, assertion):
    view = getattr(molecule, namespace)[0].constraints
    view.set(assertion)
    assert view.asserted(assertion.key) == assertion
    assert view.evidenced(assertion.key) is None
    assert view.evidenced_complete(assertion.key) is None
    assert not view.satisfies(container_type([assertion]))
    assert not view.is_compatible(container_type([assertion]))


@pytest.mark.parametrize(
    "source,namespace,expected",
    [
        pytest.param(
            '{:atoms ["C#h#T1" "F" "Cl" "Br"] '
            ':bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]] '
            ':stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs "Th0"}]}',
            "atoms",
            AtomConstraintForm.TetrahedralStereo(TetrahedralStereoForm.Stereo(StereoCoset.Lit(1))),
            id="tetrahedral",
        ),
        pytest.param(
            '{:atoms ["C#h" "C#h" "F" "Cl"] '
            ':bonds [[0 1 "2#C1"] [0 2 "1"] [1 3 "1"]] '
            ':stereo-bonds [{:site 0 :ligands [[:h 1] 3 2 [:h 0]] :attrs "Ct0"}]}',
            "bonds",
            BondConstraintForm.CisTransStereo(CisTransStereoForm.Stereo(StereoCoset.Lit(1))),
            id="cis_trans",
        ),
    ],
)
@pytest.mark.parametrize("method", ["derived", "derived_complete", "evidenced", "evidenced_complete"])
def test_entity_constraints_view_derived_frame(source, namespace, expected, method):
    molecule = Molecule.parse(source)
    original = molecule.copy()
    view = getattr(molecule, namespace)[0].constraints
    assert getattr(view, method)(expected.key) == expected
    assert molecule == original


@pytest.mark.parametrize(
    "namespace,container_type,constraint_type",
    [
        ("stereo_atoms", StereoAtomConstraintsForm, StereoAtomConstraintForm),
        ("stereo_bonds", StereoBondConstraintsForm, StereoBondConstraintForm),
    ],
)
@pytest.mark.parametrize(
    "variant,raw,expected,top",
    [
        ("LigandSymmetry",
         LigandSymmetryForm(OrientedLigandPermutation(
             LigandPermutation(Permutation([1, 0, 2, 3])), Orientation.Proper), True),
         LigandSymmetryForm(OrientedLigandPermutation(
             LigandPermutation(Permutation([1, 0, 2, 3])), Orientation.Proper), True),
         LigandSymmetryForm(OrientedLigandPermutation(
             LigandPermutation(Permutation([1, 0, 2, 3])), Orientation.Proper), BooleanForm.Undetermined())),
        ("Fluxionality",
         FluxionalityForm(LigandPermutation(Permutation([1, 0, 2, 3])), False),
         FluxionalityForm(LigandPermutation(Permutation([1, 0, 2, 3])), False),
         FluxionalityForm(LigandPermutation(Permutation([1, 0, 2, 3])), BooleanForm.Undetermined())),
        ("Topicity",
         TopicityForm(StereoLigandPair(0, 1), TopicityRelationForm.LitSet({Topicity.Homotopic})),
         TopicityForm(StereoLigandPair(0, 1), TopicityRelationForm.Lit(Topicity.Homotopic)),
         TopicityForm(StereoLigandPair(0, 1), TopicityRelationForm.Undetermined())),
    ],
)
def test_stereo_constraints_view_asserted(
    molecule, namespace, container_type, constraint_type, variant, raw, expected, top,
):
    raw = getattr(constraint_type, variant)(raw)
    expected = getattr(constraint_type, variant)(expected)
    top = getattr(constraint_type, variant)(top)
    constraints = container_type([raw])
    entity = getattr(molecule, namespace)[0]
    entity.constraints = constraints
    view = entity.constraints
    assert constraints.asserted(raw.key) == expected
    assert constraints.asserted_complete(raw.key) == expected
    for method in ("asserted", "asserted_complete", "evidenced", "evidenced_complete"):
        assert getattr(view, method)(raw.key) == expected
    assert view.derived(raw.key) == top
    assert view.derived_complete(raw.key) == top
    assert view.get(raw.key) == raw


@pytest.mark.parametrize(
    "assertion,open_value,complete_value",
    [
        (None, AromaticValenceForm.Undetermined(), AromaticValenceForm.NotAromatic()),
        (BooleanForm.Lit(True), AromaticValenceForm.Aromatic(NumForm.Undetermined()),
         AromaticValenceForm.Aromatic(NumForm.Undetermined())),
        (BooleanForm.Lit(False), AromaticValenceForm.Undetermined(), AromaticValenceForm.NotAromatic()),
    ],
)
def test_atom_constraints_view_evidenced_closure(assertion, open_value, complete_value):
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")), AtomForm(Element("C"))], bonds=[(0, 1, BondForm(1))]
    )
    if assertion is not None:
        molecule.bonds[0].constraints.aromatic = assertion
    view = molecule.atoms[0].constraints
    key = AtomConstraintKey.AromaticValence()
    assert view.asserted(key) == AtomConstraintForm.AromaticValence(AromaticValenceForm.Undetermined())
    assert view.asserted_complete(key) == AtomConstraintForm.AromaticValence(AromaticValenceForm.NotAromatic())
    assert view.derived(key) == AtomConstraintForm.AromaticValence(AromaticValenceForm.Undetermined())
    assert view.derived_complete(key) == AtomConstraintForm.AromaticValence(AromaticValenceForm.NotAromatic())
    assert view.evidenced(key) == AtomConstraintForm.AromaticValence(open_value)
    assert view.evidenced_complete(key) == AtomConstraintForm.AromaticValence(complete_value)


@pytest.mark.parametrize("method", ["satisfies", "is_compatible"])
def test_atom_constraints_view_satisfies_closure(method):
    molecule = Molecule.from_entries([AtomForm(Element("C"))])
    pattern = AtomConstraintsForm([
        AtomConstraintForm.AromaticValence(AromaticValenceForm.NotAromatic())
    ])
    assert getattr(molecule.atoms[0].constraints, method)(pattern)
    positive = AtomConstraintsForm([
        AtomConstraintForm.AromaticValence(AromaticValenceForm.Aromatic(NumForm.Lit(1)))
    ])
    assert getattr(molecule.atoms[0].constraints, method)(positive) == (method == "is_compatible")


@pytest.fixture
def triangle():
    return Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(3)],
        bonds=[(0, 1, BondForm(1)), (1, 2, BondForm(1)), (2, 0, BondForm(1))],
    )


@pytest.mark.parametrize(
    "namespace,container_type,expected",
    [
        ("atoms", AtomConstraintsForm, AtomConstraintForm.RingDegree(NumForm.Lit(2))),
        ("atoms", AtomConstraintsForm, AtomConstraintForm.RingValence(NumForm.Lit(2))),
        ("atoms", AtomConstraintsForm, AtomConstraintForm.RingMembership(
            RingMembershipForm(RingScope.All(), NumForm.Lit(1)))),
        ("atoms", AtomConstraintsForm, AtomConstraintForm.RingMembership(
            RingMembershipForm(RingScope.Size(3), NumForm.Lit(1)))),
        ("bonds", BondConstraintsForm, BondConstraintForm.RingMembership(
            RingMembershipForm(RingScope.All(), NumForm.Lit(1)))),
        ("bonds", BondConstraintsForm, BondConstraintForm.RingMembership(
            RingMembershipForm(RingScope.Size(3), NumForm.Lit(1)))),
    ],
)
def test_entity_constraints_view_with_rings(triangle, namespace, container_type, expected):
    rings = triangle.rings(RingModel(), RingConfig()).into_ring_set()
    source = getattr(triangle, namespace)[0].constraints
    key = expected.key
    references = sys.getrefcount(rings)
    view = source.with_rings(rings)
    assert sys.getrefcount(rings) == references + 1
    pattern = container_type([expected])
    for method in READINGS[2:]:
        assert getattr(view, method)(key) == expected
        with pytest.raises(ValueError, match="requires ring context"):
            getattr(source, method)(key)
    assert source.asserted(key) == container_type([]).asserted(key)
    assert source.asserted_complete(key) == container_type([]).asserted_complete(key)
    assert view.satisfies(pattern)
    assert view.is_compatible(pattern)
    for method in ("satisfies", "is_compatible"):
        with pytest.raises(ValueError, match="requires ring context"):
            getattr(source, method)(pattern)
        assert getattr(source, method)(container_type([]))
    view.set(expected)
    assert source.get(key) == expected
    assert view.asserted(key) == expected
    del rings
    assert view.derived(key) == expected


@pytest.mark.parametrize("namespace", ["atoms", "bonds"])
def test_entity_constraints_view_with_rings_error(triangle, forms, namespace):
    rings = triangle.rings(RingModel(), RingConfig()).into_ring_set()
    with pytest.raises(ValueError, match="requires molecule context"):
        forms[namespace].constraints.with_rings(rings)
    view = getattr(triangle, namespace)[0].constraints
    triangle.transact([])
    with pytest.raises(InvalidatedViewError):
        view.with_rings(rings)


@pytest.mark.parametrize(
    "namespace,container_type,conflict,ring",
    [
        ("atoms", AtomConstraintsForm, AtomConstraintForm.Valence(NumForm.Lit(1)),
         AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.LitSet(set())))),
        ("bonds", BondConstraintsForm, BondConstraintForm.Aromatic(BooleanForm.Lit(True)),
         BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.LitSet(set())))),
    ],
)
def test_entity_constraints_view_derived_ring_error(triangle, namespace, container_type, conflict, ring):
    view = getattr(triangle, namespace)[0].constraints
    view.set(ring)
    assert view.asserted(ring.key) is None
    for method in READINGS[2:]:
        with pytest.raises(ValueError, match="requires ring context"):
            getattr(view, method)(ring.key)
    pattern = container_type([conflict, ring])
    for method in ("satisfies", "is_compatible"):
        with pytest.raises(ValueError, match="requires ring context"):
            getattr(view, method)(pattern)


@pytest.mark.parametrize("namespace", ["atoms", "bonds"])
@pytest.mark.parametrize("operation", ["transact", "tracked_transact"])
def test_entity_constraints_view_with_rings_invalidation(triangle, namespace, operation):
    rings = triangle.rings(RingModel(), RingConfig()).into_ring_set()
    view = getattr(triangle, namespace)[0].constraints.with_rings(rings)
    key = (AtomConstraintKey.RingMembership(RingScope.All()) if namespace == "atoms"
           else BondConstraintKey.RingMembership(RingScope.All()))
    saved = view.derived(key)
    getattr(triangle, operation)([Edits()])
    for method in READINGS:
        with pytest.raises(InvalidatedViewError):
            getattr(view, method)(key)
    assert saved == (AtomConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(1)))
                     if namespace == "atoms" else
                     BondConstraintForm.RingMembership(RingMembershipForm(RingScope.All(), NumForm.Lit(1))))
    assert list(rings.ids()) == [0]


@pytest.mark.parametrize("namespace", ["atoms", "bonds"])
@pytest.mark.parametrize("method", ["edit", "apply", "tracked_apply"])
def test_entity_constraints_view_with_rings_edit(triangle, namespace, method):
    rings = triangle.rings(RingModel(), RingConfig()).into_ring_set()
    view = getattr(triangle, namespace)[0].constraints.with_rings(rings)
    key = (AtomConstraintKey.RingMembership(RingScope.All()) if namespace == "atoms"
           else BondConstraintKey.RingMembership(RingScope.All()))
    before = [getattr(view, reading)(key) for reading in READINGS]
    if method == "edit":
        assert triangle.edit().finish() == triangle
    else:
        getattr(triangle, method)(Edits())
    assert [getattr(view, reading)(key) for reading in READINGS] == before
    assert list(rings.ids()) == [0]


def test_atom_constraints_view_with_rings_mutation(triangle):
    rings = triangle.rings(RingModel(), RingConfig()).into_ring_set()
    view = triangle.atoms[0].constraints.with_rings(rings)
    key = AtomConstraintKey.RingValence()
    saved = view.derived(key)
    triangle.bonds[0].order = NumForm.LitSet({2})
    assert view.derived(key) == AtomConstraintForm.RingValence(NumForm.Lit(3))
    assert triangle.bonds[0].order == NumForm.LitSet({2})
    triangle.bonds[0].order = NumForm.LitSet(set())
    assert view.derived(key) is None
    assert saved == AtomConstraintForm.RingValence(NumForm.Lit(2))
    empty = Molecule().rings(RingModel(), RingConfig()).into_ring_set()
    replacement = view.with_rings(empty)
    assert replacement.derived(key) == AtomConstraintForm.RingValence(NumForm.Lit(0))
    assert view.derived(key) is None


def test_atom_constraints_view_with_rings_lifetime():
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(3)],
        bonds=[(0, 1, BondForm(1)), (1, 2, BondForm(1)), (2, 0, BondForm(1))],
    )
    rings = molecule.rings(RingModel(), RingConfig()).into_ring_set()
    atom = molecule.atoms[0].constraints.with_rings(rings)
    bond = molecule.bonds[0].constraints.with_rings(rings)
    del molecule, rings
    assert atom.derived(AtomConstraintKey.RingDegree()) == AtomConstraintForm.RingDegree(NumForm.Lit(2))
    assert bond.derived(BondConstraintKey.RingMembership(RingScope.All())) == BondConstraintForm.RingMembership(
        RingMembershipForm(RingScope.All(), NumForm.Lit(1))
    )


@pytest.mark.parametrize("method", ["derived", "derived_complete", "evidenced", "evidenced_complete"])
def test_dative_bond_constraints_view_derived_ring(molecule, method):
    key = DativeBondConstraintKey.RingMembership(RingScope.All())
    expected = DativeBondConstraintForm.RingMembership(
        RingMembershipForm(RingScope.All(), NumForm.Undetermined())
    )
    assert getattr(molecule.dative_bonds[0].constraints, method)(key) == expected


@pytest.mark.parametrize(
    "name,expected",
    [
        ("valence", 2), ("donated_pairs", 0), ("accepted_pairs", 0),
        ("total_degree", 2), ("total_hydrogens", 1), ("total_valence", 3),
        ("covalence", 3), ("heavy_atom_degree", 1), ("heavy_atom_valence", 2),
    ],
)
def test_atom_view_projection(name, expected):
    molecule = Molecule.from_entries(
        [AtomForm(Element("C"), implicit_hydrogens=1), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(NumForm.LitSet({2})))],
    )
    atom = molecule.atoms[0]
    value = getattr(atom, name)
    assert value == NumForm.Lit(expected)
    molecule.atoms[0].charge = NumForm.LitSet(set())
    assert getattr(atom, name) == NumForm.Lit(expected)
    assert molecule.bonds[0].order == NumForm.LitSet({2})
    assert molecule.edit().finish() == molecule
    assert getattr(atom, name) == value
    molecule.transact([])
    with pytest.raises(InvalidatedViewError):
        getattr(atom, name)
    assert value == NumForm.Lit(expected)


@pytest.mark.parametrize("name", ["valence", "total_valence", "covalence", "heavy_atom_valence"])
@pytest.mark.parametrize("order,expected", [(NumForm.Undetermined(), NumForm.Undetermined()), (NumForm.LitSet(set()), None)])
def test_atom_view_valence(name, order, expected):
    molecule = Molecule.from_entries(
        [AtomForm(Element("C"), implicit_hydrogens=0), AtomForm(Element("N"))],
        bonds=[(0, 1, BondForm(order))],
    )
    assert getattr(molecule.atoms[0], name) == expected
    assert molecule.bonds[0].order == order


@pytest.mark.parametrize("name", ["total_degree", "total_hydrogens", "total_valence", "covalence"])
@pytest.mark.parametrize("hydrogens,expected", [(NumForm.Undetermined(), NumForm.Undetermined()), (NumForm.LitSet(set()), None)])
def test_atom_view_total_hydrogens(name, hydrogens, expected):
    molecule = Molecule.from_entries([AtomForm(Element("C"), implicit_hydrogens=hydrogens)])
    assert getattr(molecule.atoms[0], name) == expected
    assert molecule.atoms[0].implicit_hydrogens == hydrogens


@pytest.mark.parametrize(
    "element,degree,valence,hydrogens",
    [
        (ElementForm.LitSet({Element("H")}), NumForm.Lit(0), NumForm.Lit(0), NumForm.Lit(1)),
        (ElementForm.Undetermined(), NumForm.Lit(1), None, NumForm.Undetermined()),
        (ElementForm.LitSet(set()), None, None, None),
    ],
)
def test_atom_view_heavy_atom_valence(element, degree, valence, hydrogens):
    molecule = Molecule.from_entries(
        [AtomForm(Element("C"), implicit_hydrogens=0), AtomForm(element)],
        bonds=[(0, 1, BondForm(NumForm.LitSet(set())))],
    )
    atom = molecule.atoms[0]
    assert atom.heavy_atom_degree == degree
    assert atom.heavy_atom_valence == valence
    assert atom.total_hydrogens == hydrogens


@pytest.mark.parametrize("order,expected", [
    (NumForm.LitSet({2}), NumForm.Lit(2)),
    (NumForm.Undetermined(), NumForm.Undetermined()),
    (NumForm.LitSet(set()), None),
])
def test_atom_view_donated_pairs(order, expected):
    molecule = Molecule.from_entries(
        [AtomForm(Element("N")), AtomForm(Element("B"))],
        dative_bonds=[([0], 1, DativeBondForm(order))],
    )
    assert molecule.atoms[0].donated_pairs == expected
    assert molecule.atoms[1].accepted_pairs == expected
    assert molecule.atoms[0].accepted_pairs == NumForm.Lit(0)
    assert molecule.atoms[1].donated_pairs == NumForm.Lit(0)


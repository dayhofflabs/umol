import sys

import pytest

from umol import (
    AtomForm,
    BondForm,
    ConsumedError,
    Edits,
    Element,
    InvalidatedViewError,
    Molecule,
    NumForm,
    RelevantCycleEnumerationAlgorithm,
    RingAtomView,
    RingBondView,
    RingConfig,
    RingConnection,
    RingGraph,
    RingModel,
    RingRelation,
    RingScope,
    RingSet,
    RingSetKind,
    RingView,
    RingViews,
    SimpleCycleEnumerationAlgorithm,
)


@pytest.fixture
def fused_rings():
    return Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(5)],
        bonds=[
            (0, 1, BondForm(1)),
            (1, 2, BondForm(1)),
            (2, 0, BondForm(1)),
            (1, 3, BondForm(1)),
            (3, 2, BondForm(1)),
            (0, 4, BondForm(1)),
        ],
    )


@pytest.mark.parametrize("kind", [RingSetKind.Simple, RingSetKind.Relevant])
def test_ring_set_kind_repr(kind):
    assert eval(repr(kind), {"RingSetKind": RingSetKind}) == kind
    assert {kind, kind} == {kind}


def test_ring_model_new():
    model = RingModel()
    assert model.kind == RingSetKind.Relevant
    assert model.max_ring_size == 22
    assert model == RingModel(kind=RingSetKind.Relevant, max_ring_size=22)
    assert eval(repr(model), {"RingModel": RingModel, "RingSetKind": RingSetKind}) == model


def test_ring_model_copy():
    model = RingModel()
    copied = model.copy()
    copied.kind = RingSetKind.Simple
    copied.max_ring_size = 4
    assert copied == RingModel(kind=RingSetKind.Simple, max_ring_size=4)
    assert model == RingModel()
    with pytest.raises(TypeError):
        hash(model)


def test_ring_config_default():
    config = RingConfig()

    assert (
        config.simple_cycle_algorithm
        == SimpleCycleEnumerationAlgorithm.ReadTarjan()
    )
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )
    assert config == RingConfig()


@pytest.mark.parametrize(
    "kwargs",
    [
        {
            "simple_cycle_algorithm": (
                SimpleCycleEnumerationAlgorithm.ReadTarjan()
            ),
        },
        {
            "relevant_cycle_algorithm": (
                RelevantCycleEnumerationAlgorithm.Vismara()
            ),
        },
        {
            "simple_cycle_algorithm": (
                SimpleCycleEnumerationAlgorithm.ReadTarjan()
            ),
            "relevant_cycle_algorithm": (
                RelevantCycleEnumerationAlgorithm.Vismara()
            ),
        },
    ],
)
def test_ring_config_new(kwargs):
    config = RingConfig(**kwargs)

    assert (
        config.simple_cycle_algorithm
        == SimpleCycleEnumerationAlgorithm.ReadTarjan()
    )
    assert (
        config.relevant_cycle_algorithm
        == RelevantCycleEnumerationAlgorithm.Vismara()
    )


def test_ring_config_new_error():
    with pytest.raises(TypeError):
        RingConfig(
            SimpleCycleEnumerationAlgorithm.ReadTarjan(),
            RelevantCycleEnumerationAlgorithm.Vismara(),
        )


def test_ring_config_repr():
    assert repr(RingConfig()) == (
        "RingConfig("
        "simple_cycle_algorithm=SimpleCycleEnumerationAlgorithm.ReadTarjan(), "
        "relevant_cycle_algorithm="
        "RelevantCycleEnumerationAlgorithm.Vismara())"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        (
            "simple_cycle_algorithm",
            SimpleCycleEnumerationAlgorithm.ReadTarjan(),
        ),
        (
            "relevant_cycle_algorithm",
            RelevantCycleEnumerationAlgorithm.Vismara(),
        ),
    ],
)
def test_ring_config_mutation(field, value):
    config = RingConfig()

    with pytest.raises(AttributeError):
        setattr(config, field, value)


@pytest.mark.parametrize(
    "relation",
    [
        RingRelation.Identical, RingRelation.Disjoint, RingRelation.Spiro,
        RingRelation.Fused, RingRelation.Bridged, RingRelation.MultiSpiro,
        RingRelation.Noncontiguous,
    ],
)
def test_ring_connection_new(relation):
    connection = RingConnection(0, 1, relation)
    assert (connection.first, connection.second, connection.relation) == (0, 1, relation)
    assert not connection.readonly
    assert connection == RingConnection(0, 1, relation)
    assert eval(repr(connection), {"RingConnection": RingConnection, "RingRelation": RingRelation}) == connection
    assert {relation, relation} == {relation}


def test_ring_connection_copy():
    connection = RingConnection(0, 1, RingRelation.Fused)
    copied = connection.copy()
    copied.first = 2
    copied.second = 3
    copied.relation = RingRelation.Spiro
    assert copied == RingConnection(2, 3, RingRelation.Spiro)
    assert connection == RingConnection(0, 1, RingRelation.Fused)
    with pytest.raises(TypeError):
        hash(connection)


def test_ring_graph_edges(fused_rings):
    rings = fused_rings.rings(RingModel(), RingConfig()).into_ring_set()
    graph = rings.graph()
    assert graph.neighbors(0) == [(1, RingRelation.Fused)]
    assert graph.neighbors(1) == [(0, RingRelation.Fused)]
    assert graph.neighbors(2) == []
    assert graph.relation(0, 1) == RingRelation.Fused
    assert graph.relation(0, 0) == RingRelation.Identical
    assert graph.relation(0, 2) == RingRelation.Disjoint
    edges = graph.edges
    edge = edges[0]
    assert list(edges) == [RingConnection(0, 1, RingRelation.Fused)]
    assert edges[-1] == edge
    assert edge.readonly
    for field, value in [("first", 1), ("second", 0), ("relation", RingRelation.Spiro)]:
        with pytest.raises(TypeError, match="read-only"):
            setattr(edge, field, value)
    with pytest.raises(TypeError):
        edges[0] = RingConnection(0, 1, RingRelation.Spiro)
    with pytest.raises(IndexError, match="ring connection index out of range"):
        _ = edges[1]
    with pytest.raises(IndexError):
        _ = edges[-2]
    copied_edge = edge.copy()
    copied_edge.relation = RingRelation.Spiro
    assert edge.relation == RingRelation.Fused
    copied_graph = graph.copy()
    iterator = iter(edges)
    del graph, edges, rings, fused_rings
    assert edge == RingConnection(0, 1, RingRelation.Fused)
    assert list(iterator) == [RingConnection(0, 1, RingRelation.Fused)]
    assert list(copied_graph.edges) == [RingConnection(0, 1, RingRelation.Fused)]


def test_ring_graph_edges_iteration(fused_rings):
    graph = fused_rings.rings(RingModel(), RingConfig()).into_ring_set().graph()
    edges = graph.edges
    references = sys.getrefcount(graph)
    iterator = iter(edges)
    assert sys.getrefcount(graph) == references
    edge = next(iterator)
    assert sys.getrefcount(graph) == references + 1
    assert edge == RingConnection(0, 1, RingRelation.Fused)
    with pytest.raises(StopIteration):
        next(iterator)


def test_ring_set_get(fused_rings):
    rings = fused_rings.rings(RingModel(), RingConfig()).into_ring_set()
    copied = rings.copy()
    del fused_rings
    assert rings.kind == RingSetKind.Relevant
    assert rings.max_ring_size == 22
    assert list(rings.ids()) == [0, 1]
    assert [(ring.id, sorted(ring.atoms), sorted(ring.bonds)) for ring in rings] == [
        (0, [0, 1, 2], [0, 1, 2]),
        (1, [1, 2, 3], [1, 3, 4]),
    ]
    assert rings.get(0).id == 0
    assert rings.get(2) is None
    assert 0 in rings
    assert 2 not in rings
    assert rings[-1].id == 1
    with pytest.raises(IndexError):
        _ = rings[2]
    with pytest.raises(IndexError):
        _ = rings[-3]
    with pytest.raises(TypeError):
        rings[0] = rings[1]
    assert rings.shared_atoms(0, 1) == [1, 2]
    assert rings.shared_bonds(0, 1) == [1]
    assert rings.shared_atoms(0, 2) == []
    assert rings.shared_bonds(2, 0) == []
    assert rings.relation(0, 1) == RingRelation.Fused
    assert rings.relation(0, 0) == RingRelation.Identical
    assert rings.relation(0, 2) == RingRelation.Disjoint
    assert rings.are_fused(0, 1)
    assert not rings.are_spiro(0, 1)
    assert not rings.are_bridged(0, 1)
    assert rings.spiro_neighbors(0) == []
    assert rings.fused_neighbors(0) == [1]
    assert rings.bridged_neighbors(0) == []
    assert rings.fused_component(0) == [0, 1]
    assert rings.fused_components() == [[0, 1]]
    assert rings.contains_atom(0)
    assert not rings.contains_atom(4)
    assert rings.atom_smallest_ring_size(0) == 3
    assert rings.atom_smallest_ring_size(4) is None
    assert rings.contains_bond(0)
    assert not rings.contains_bond(5)
    assert rings.bond_smallest_ring_size(0) == 3
    assert rings.bond_smallest_ring_size(5) is None
    del rings
    assert list(copied.ids()) == [0, 1]
    assert copied.shared_atoms(0, 1) == [1, 2]
    assert copied.shared_bonds(0, 1) == [1]


@pytest.mark.parametrize(
    ("atom_count", "bonds", "relation", "neighbors", "components"),
    [
        (6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)],
         RingRelation.Disjoint, [], [[0], [1]]),
        (5, [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)],
         RingRelation.Spiro, [(1, RingRelation.Spiro)], [[0], [1]]),
        (5, [(0, 1), (1, 3), (0, 2), (2, 3), (0, 4), (4, 3)],
         RingRelation.Bridged, [(1, RingRelation.Bridged), (2, RingRelation.Bridged)],
         [[0], [1], [2]]),
    ],
)
def test_ring_set_relation(atom_count, bonds, relation, neighbors, components):
    molecule = Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(atom_count)],
        bonds=[(a, b, BondForm(1)) for a, b in bonds],
    )
    rings = molecule.rings(RingModel(), RingConfig()).into_ring_set()
    assert rings.relation(0, 1) == relation
    assert rings.graph().neighbors(0) == neighbors
    assert rings.are_spiro(0, 1) == (relation == RingRelation.Spiro)
    assert rings.are_fused(0, 1) == (relation == RingRelation.Fused)
    assert rings.are_bridged(0, 1) == (relation == RingRelation.Bridged)
    assert rings.spiro_neighbors(0) == ([1] if relation == RingRelation.Spiro else [])
    assert rings.bridged_neighbors(0) == ([1, 2] if relation == RingRelation.Bridged else [])
    assert rings.fused_components() == components


@pytest.mark.parametrize(
    ("kind", "limit", "expected"),
    [
        (RingSetKind.Relevant, 22, [[0, 1, 2], [1, 2, 3]]),
        (RingSetKind.Simple, 3, [[0, 1, 2], [1, 2, 3]]),
        (RingSetKind.Simple, 4, [[0, 1, 2], [0, 1, 2, 3], [1, 2, 3]]),
        (RingSetKind.Relevant, 2, []),
    ],
)
def test_molecule_rings(fused_rings, kind, limit, expected):
    model = RingModel(kind=kind, max_ring_size=limit)
    rings = fused_rings.rings(model, RingConfig())
    assert rings.kind == kind
    assert rings.max_ring_size == limit
    assert sorted(sorted(ring.atoms) for ring in rings) == expected
    model.kind = RingSetKind.Simple
    model.max_ring_size = 0
    assert rings.kind == kind
    assert rings.max_ring_size == limit


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
def test_ring_views_from_ring_set_error(fused_rings, borrowed):
    source = fused_rings.rings(RingModel(), RingConfig()).into_ring_set()
    fused_rings.resolve()
    with pytest.raises(ConsumedError, match="Molecule has been consumed"):
        if borrowed:
            RingViews.from_ring_set(source, fused_rings)
        else:
            fused_rings.rings(RingModel(), RingConfig())


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
def test_ring_views_from_ring_set(fused_rings, borrowed):
    rings = fused_rings.rings(RingModel(), RingConfig())
    if borrowed:
        source = rings.into_ring_set()
        rings = RingViews.from_ring_set(source, fused_rings)
    assert rings.kind == RingSetKind.Relevant
    assert rings.max_ring_size == 22
    assert len(rings) == 2
    assert list(rings.ids()) == [0, 1]
    assert [(ring.id, sorted(ring.atoms), sorted(ring.bonds)) for ring in rings] == [
        (0, [0, 1, 2], [0, 1, 2]),
        (1, [1, 2, 3], [1, 3, 4]),
    ]
    assert rings.get(0).id == 0
    assert rings.get(2) is None
    assert 0 in rings
    assert 2 not in rings
    assert rings[-1].id == 1
    with pytest.raises(IndexError):
        _ = rings[2]
    with pytest.raises(IndexError):
        _ = rings[-3]
    with pytest.raises(TypeError):
        rings[0] = rings[1]
    atom = rings.atom(1)
    bond = rings.bond(1)
    assert atom.is_in_ring()
    assert bond.is_in_ring()
    assert [ring.id for ring in atom.rings()] == [0, 1]
    assert [ring.id for ring in bond.rings()] == [0, 1]
    assert atom.ring_membership(RingScope.All()) == NumForm.Lit(2)
    assert bond.ring_membership(RingScope.Size(3)) == NumForm.Lit(2)
    assert atom.ring_count == bond.ring_count == NumForm.Lit(2)
    assert atom.ring_size_count(3) == bond.ring_size_count(3) == NumForm.Lit(2)
    assert atom.ring_size_count(4) == bond.ring_size_count(4) == NumForm.Lit(0)
    assert atom.smallest_ring_size == bond.smallest_ring_size == 3
    assert atom.ring_degree == NumForm.Lit(3)
    assert atom.ring_valence == NumForm.Lit(3)
    assert not rings.atom(4).is_in_ring()
    assert not rings.bond(5).is_in_ring()
    assert list(rings.atom(4).rings()) == []
    assert list(rings.bond(5).rings()) == []
    assert rings.atom(4).smallest_ring_size is None
    assert rings.bond(5).smallest_ring_size is None
    assert rings.atom(4).ring_degree == NumForm.Lit(0)
    assert rings.atom(4).ring_valence == NumForm.Lit(0)


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
@pytest.mark.parametrize(
    ("order", "expected"),
    [
        (NumForm.Lit(2), NumForm.Lit(3)),
        (NumForm.Undetermined(), NumForm.Undetermined()),
        (NumForm.LitSet({2}), NumForm.Lit(3)),
        (NumForm.LitSet(set()), None),
    ],
    ids=["literal", "top", "singleton", "bottom"],
)
def test_ring_atom_view_ring_valence(fused_rings, borrowed, order, expected):
    rings = fused_rings.rings(RingModel(), RingConfig())
    if borrowed:
        source = rings.into_ring_set()
        rings = RingViews.from_ring_set(source, fused_rings)
    atom = rings.atom(0)
    ring = rings[0]
    assert atom.ring_valence == NumForm.Lit(2)
    fused_rings.bonds[0].order = order
    assert atom.ring_valence == expected
    assert atom.ring_degree == NumForm.Lit(2)
    assert atom.ring_count == NumForm.Lit(1)
    assert sorted(ring.atoms) == [0, 1, 2]
    assert sorted(ring.bonds) == [0, 1, 2]
    assert fused_rings.bonds[0].order == order


@pytest.mark.parametrize(
    ("atom_count", "bonds"),
    [(0, []), (2, [(0, 1)]), (5, [(0, 4)]), (6, [(0, 5)])],
    ids=["empty", "smaller", "different_topology", "larger"],
)
def test_ring_views_from_ring_set_context(fused_rings, atom_count, bonds):
    source = fused_rings.rings(RingModel(), RingConfig()).into_ring_set()
    context = Molecule.from_entries(
        [AtomForm(Element("C")) for _ in range(atom_count)],
        bonds=[(a, b, BondForm(1)) for a, b in bonds],
    )
    rings = RingViews.from_ring_set(source, context)
    for atom_id in range(max(5, atom_count) + 1):
        atom = rings.atom(atom_id)
        _ = atom.ring_degree, atom.ring_valence
    for bond_id in range(7):
        bond = rings.bond(bond_id)
        _ = bond.ring_count, bond.smallest_ring_size, list(bond.rings())
    assert [(ring.id, sorted(ring.atoms), sorted(ring.bonds)) for ring in rings] == [
        (0, [0, 1, 2], [0, 1, 2]),
        (1, [1, 2, 3], [1, 3, 4]),
    ]


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
def test_ring_views_into_ring_set(fused_rings, borrowed):
    rings = fused_rings.rings(RingModel(), RingConfig())
    if borrowed:
        source = rings.into_ring_set()
        rings = RingViews.from_ring_set(source, fused_rings)
    expected = [(ring.id, list(ring.atoms), list(ring.bonds)) for ring in rings]
    alias = rings
    extracted = rings.into_ring_set()
    if borrowed:
        assert source.shared_atoms(0, 1) == [1, 2]
        del source
    del fused_rings
    assert [(ring.id, list(ring.atoms), list(ring.bonds)) for ring in extracted] == expected
    assert list(extracted.graph().edges) == [RingConnection(0, 1, RingRelation.Fused)]
    with pytest.raises(ConsumedError, match="RingViews has been consumed"):
        len(alias)


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
def test_ring_views_from_ring_set_lifetime(fused_rings, borrowed):
    rings = fused_rings.rings(RingModel(), RingConfig())
    if borrowed:
        source = rings.into_ring_set()
        rings = RingViews.from_ring_set(source, fused_rings)
        del source
    ring = rings[0]
    atom = rings.atom(1)
    bond = rings.bond(1)
    iterator = iter(rings)
    del rings, fused_rings
    assert sorted(ring.atoms) == [0, 1, 2]
    assert atom.ring_valence == NumForm.Lit(3)
    assert bond.ring_count == NumForm.Lit(2)
    assert [ring.id for ring in iterator] == [0, 1]


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
@pytest.mark.parametrize("method", ["edit", "apply", "tracked_apply"])
def test_ring_views_from_ring_set_edit(fused_rings, borrowed, method):
    original = fused_rings.copy()
    source = fused_rings.rings(RingModel(), RingConfig()).into_ring_set()
    rings = (RingViews.from_ring_set(source, fused_rings) if borrowed
             else fused_rings.rings(RingModel(), RingConfig()))
    ring = rings[0]
    members = ring.atoms
    expected_members = list(members)
    atom = rings.atom(1)
    bond = rings.bond(1)
    iterator = iter(rings)
    edits = Edits()
    edits.remove_topology([4], [])
    if method == "edit":
        editor = fused_rings.edit()
        editor.apply(edits)
        product = editor.finish()
    else:
        result = getattr(fused_rings, method)(edits)
        product = result[0] if method == "tracked_apply" else result

    assert fused_rings == original
    assert product == Molecule.parse(
        '{:atoms ["C" "C" "C" "C"] '
        ':bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"] [1 3 "1"] [3 2 "1"]]}'
    )
    assert list(members) == expected_members
    assert ring.id == 0
    assert sorted(ring.copy().atoms) == [0, 1, 2]
    assert sorted(ring.copy().bonds) == [0, 1, 2]
    assert atom.ring_count == NumForm.Lit(2)
    assert bond.ring_count == NumForm.Lit(2)
    assert [item.id for item in iterator] == [0, 1]
    assert list(rings.ids()) == [0, 1]


@pytest.mark.parametrize("borrowed", [False, True], ids=["owned", "borrowed"])
@pytest.mark.parametrize("operation", ["extract", "resolve", "transact", "tracked_transact", "combine"])
def test_ring_views_from_ring_set_invalidation(fused_rings, borrowed, operation):
    source = fused_rings.rings(RingModel(), RingConfig()).into_ring_set()
    rings = (
        RingViews.from_ring_set(source, fused_rings) if borrowed
        else fused_rings.rings(RingModel(), RingConfig())
    )
    ring = rings[0]
    members = ring.atoms
    atom = rings.atom(1)
    bond = rings.bond(1)
    iterators = [
        iter(rings), rings.ids(), iter(members), atom.rings(), bond.rings(),
    ]
    iterated = next(iterators[0])
    exhausted = [
        iter(rings), rings.ids(), iter(members), atom.rings(), bond.rings(),
    ]
    for iterator in exhausted:
        list(iterator)
    iterators.extend(exhausted)
    if operation == "extract":
        extracted = rings.into_ring_set()
        assert extracted.shared_bonds(0, 1) == [1]
        root_error = ConsumedError
    else:
        if operation == "resolve":
            fused_rings.resolve()
        elif operation == "combine":
            fused_rings.combine_from(Molecule())
        else:
            getattr(fused_rings, operation)([Edits()])
        root_error = InvalidatedViewError
    for access in [
        lambda: len(rings), lambda: repr(rings), lambda: rings.kind,
        lambda: rings.max_ring_size, lambda: rings.get(0), lambda: rings[0],
        lambda: 0 in rings, lambda: rings.ids(), lambda: iter(rings),
        lambda: rings.atom(0), lambda: rings.bond(0), lambda: rings.into_ring_set(),
    ]:
        with pytest.raises(root_error):
            access()
    for access in [
        lambda: ring.id, lambda: iterated.id, lambda: setattr(ring, "id", 1), lambda: repr(ring),
        lambda: ring.atoms, lambda: ring.bonds, lambda: ring.is_empty(),
        lambda: ring.copy(), lambda: len(ring), lambda: len(members),
        lambda: members[0], lambda: iter(members),
        lambda: atom.is_in_ring(), lambda: atom.rings(), lambda: atom.ring_count,
        lambda: atom.ring_membership(RingScope.All()), lambda: atom.ring_size_count(3),
        lambda: atom.smallest_ring_size, lambda: atom.ring_degree, lambda: atom.ring_valence,
        lambda: repr(atom), lambda: bond.is_in_ring(), lambda: bond.rings(),
        lambda: bond.ring_count, lambda: bond.ring_membership(RingScope.All()),
        lambda: bond.ring_size_count(3), lambda: bond.smallest_ring_size, lambda: repr(bond),
    ]:
        with pytest.raises(InvalidatedViewError):
            access()
    for iterator in iterators:
        with pytest.raises(InvalidatedViewError):
            iter(iterator)
        with pytest.raises(InvalidatedViewError):
            next(iterator)
    assert source.shared_atoms(0, 1) == [1, 2]
    assert list(source.ids()) == [0, 1]


def test_ring_view_atoms(fused_rings):
    ring = fused_rings.rings(RingModel(), RingConfig())[0]
    atoms, bonds = ring.atoms, ring.bonds
    assert len(ring) == 3
    assert not ring.is_empty()
    assert sorted(atoms) == [0, 1, 2]
    assert sorted(bonds) == [0, 1, 2]
    assert atoms[-1] == atoms[2]
    assert bonds[-1] == bonds[2]
    with pytest.raises(IndexError):
        _ = atoms[3]
    with pytest.raises(IndexError):
        _ = bonds[-4]
    with pytest.raises(TypeError):
        atoms[0] = 3
    with pytest.raises(TypeError):
        bonds[0] = 3
    with pytest.raises(AttributeError):
        ring.atoms = [0, 1, 2]
    copied = ring.copy()
    ring.id = 7
    assert ring.id == 7
    assert copied.id == 0
    assert sorted(ring.atoms) == [0, 1, 2]
    iterator = iter(atoms)
    del atoms, ring, fused_rings
    assert sorted(iterator) == [0, 1, 2]
    assert sorted(bonds) == [0, 1, 2]


@pytest.mark.parametrize("standalone", [False, True], ids=["views", "set"])
def test_ring_views_iter(fused_rings, standalone):
    rings = fused_rings.rings(RingModel(), RingConfig())
    if standalone:
        rings = rings.into_ring_set()
    references = sys.getrefcount(rings)
    iterator = iter(rings)
    assert sys.getrefcount(rings) == references + 1
    first = next(iterator)
    assert sys.getrefcount(rings) == references + 2
    assert first.id == 0
    del first
    assert sys.getrefcount(rings) == references + 1
    assert [ring.id for ring in iterator] == [1]
    ids = rings.ids()
    references = sys.getrefcount(rings)
    assert next(ids) == 0
    assert sys.getrefcount(rings) == references
    assert list(ids) == [1]


@pytest.mark.parametrize("entity", ["atom", "bond"])
def test_ring_atom_view_rings(fused_rings, entity):
    rings = fused_rings.rings(RingModel(), RingConfig())
    view = getattr(rings, entity)(1)
    references = sys.getrefcount(rings)
    iterator = view.rings()
    assert sys.getrefcount(rings) == references + 1
    first = next(iterator)
    assert sys.getrefcount(rings) == references + 2
    assert first.id == 0
    assert [ring.id for ring in iterator] == [1]


@pytest.mark.parametrize("ring_type", [RingGraph, RingSet, RingViews, RingView, RingAtomView, RingBondView])
def test_ring_views_new_error(ring_type):
    with pytest.raises(TypeError):
        ring_type()

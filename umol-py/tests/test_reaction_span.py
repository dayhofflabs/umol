from umol import MoleculeRemapping
import pytest
from umol import (
    AtomForm,
    BondForm,
    BondDelta,
    BondFieldChange,
    ContradictionError,
    Delta,
    Deltas,
    Element,
    Entity,
    MetadataError,
    Molecule,
    MoleculeCorrespondence,
    MoleculeDefaults,
    MoleculeMetadata,
    ParseError,
    Reaction,
    ReactionSpan,
    NumForm,
)


def test_reaction_span_parse():
    span = ReactionSpan.parse(r'{:atoms ["C" {:add "O"}]}')

    assert span == ReactionSpan.from_entries(
        [
            (AtomForm(Element("C")), AtomForm(Element("C"))),
            (None, AtomForm(Element("O"))),
        ]
    )


def test_reaction_span_parse_error():
    with pytest.raises(
        ParseError,
        match="^EDN parse: unexpected token 'n' at byte 0$",
    ):
        ReactionSpan.parse("not edn")


def test_reaction_span_parse_defaults():
    span = ReactionSpan.parse(
        r'{:atoms ["C#h4#v0#d0#t0#a!#m!"]}',
        defaults=MoleculeDefaults.concrete(),
    )

    atom = AtomForm.parse("C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!")
    assert span == ReactionSpan.from_entries([(atom, atom)])


def test_reaction_span_parse_with_metadata_roundtrip():
    span, metadata = ReactionSpan.parse_with_metadata(
        r'{:atoms [[:carbon :x] {:add "O"}] '
        r':bonds [{:add [0 1 :single]}] :atom-aliases [:x "C"]}'
    )

    rendered = span.render_with_metadata(metadata)
    reparsed, reparsed_metadata = ReactionSpan.parse_with_metadata(rendered)

    assert rendered == (
        r'{:atom-aliases [:x "C"] :atoms [[:carbon :x] {:add "O"}] '
        r':bonds [{:add [:carbon 1 :single]}]}'
    )
    assert (reparsed, reparsed_metadata) == (span, metadata)


def test_reaction_span_render():
    span = ReactionSpan.from_entries(
        [
            (AtomForm(Element("C")), AtomForm(Element("C"))),
            (None, AtomForm(Element("O"))),
        ]
    )

    assert span.render() == r'{:atoms ["C" {:add "O"}]}'
    assert str(span) == span.render()


def test_reaction_span_render_with_metadata_error():
    span = ReactionSpan.from_entries(
        [(AtomForm(Element("C")), AtomForm(Element("C")))]
    )
    metadata = MoleculeMetadata()
    metadata.set_keyword(Entity.Atom(1), "outside")

    with pytest.raises(
        MetadataError,
        match="^metadata entity is out of range: atom 1$",
    ):
        span.render_with_metadata(metadata)


def test_reaction_span_from_entries():
    lhs = AtomForm(Element("C"), charge=NumForm.Lit(1))
    rhs = AtomForm(Element("C"), charge=NumForm.LitSet({1}))

    span = ReactionSpan.from_entries([(lhs, rhs)])

    assert span == ReactionSpan.from_entries([(lhs, lhs)])


def test_reaction_span_from_entries_error():
    with pytest.raises(
        ValueError,
        match="^reaction span entry is absent from both sides$",
    ):
        ReactionSpan.from_entries([(None, None)])


def test_reaction_span_from_entries_reference_error():
    with pytest.raises(
        ValueError,
        match="^reaction span entries reference unavailable atom 1$",
    ):
        ReactionSpan.from_entries(
            [(AtomForm(Element("C")), AtomForm(Element("C")))],
            bonds=[(0, 1, (BondForm(1), BondForm(1)))],
        )


def test_reaction_span_lhs():
    span = ReactionSpan.parse(
        r'{:atoms ["C" {:remove "O"} {:add "N"}] '
        r':bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}'
    )

    assert span.lhs() == Molecule.parse(
        r'{:atoms ["C" "O"] :bonds [[0 1 :single]]}'
    )


def test_reaction_span_rhs():
    span = ReactionSpan.parse(
        r'{:atoms ["C" {:remove "O"} {:add "N"}] '
        r':bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}'
    )

    assert span.rhs() == Molecule.parse(
        r'{:atoms ["C" "N"] :bonds [[0 1 :double]]}'
    )


def test_reaction_span_correspondence():
    span = ReactionSpan.parse(
        r'{:atoms ["C" {:remove "O"} {:add "N"}] '
        r':bonds [{:remove [0 1 :single]} {:add [0 2 :double]}]}'
    )

    correspondence = span.correspondence()

    assert correspondence.atoms.matched_pairs == [(0, 0)]
    assert correspondence.atoms.left_count == 2
    assert correspondence.atoms.right_count == 2
    assert correspondence.atoms.left_unmatched == [1]
    assert correspondence.atoms.right_unmatched == [1]
    assert correspondence.bonds.matched_pairs == []
    assert correspondence.bonds.left_count == 1
    assert correspondence.bonds.right_count == 1
    assert correspondence.bonds.left_unmatched == [0]
    assert correspondence.bonds.right_unmatched == [0]


def test_reaction_span_to_reaction_roundtrip():
    span = ReactionSpan.parse(
        r'{:atoms ["C" {:add "O"}] :bonds [{:add [0 1 :single]}]}'
    )

    reaction = span.to_reaction()

    assert reaction.lhs == span.lhs()
    assert reaction.to_reaction_span() == span


def test_reaction_span_canonicalize():
    source = ReactionSpan.parse(
        '{:atoms [{:add "O"} {:modify ["C" "N"]} {:remove "F"} "Cl"] '
        ':bonds [{:remove [2 3 :single]} {:add [0 1 :double]} '
        '{:modify [1 3 [:single :double]]}]}'
    )
    expected = ReactionSpan.parse(
        '{:atoms ["Cl" {:remove "F"} {:modify ["C" "N"]} {:add "O"}] '
        ':bonds [{:remove [0 1 :single]} {:modify [0 2 [:single :double]]} '
        '{:add [2 3 :double]}]}'
    )

    canonical = source.canonicalize()

    assert canonical is not source
    assert canonical == expected
    assert source != expected
    assert source.canonical_eq(expected)


def test_reaction_span_tracked_canonicalize():
    source = ReactionSpan.parse(
        '{:atoms [{:add "O"} {:modify ["C" "N"]} {:remove "F"} "Cl"] '
        ':bonds [{:remove [2 3 :single]} {:add [0 1 :double]} '
        '{:modify [1 3 [:single :double]]}]}'
    )

    canonical, remapping = source.tracked_canonicalize()

    assert canonical == source.canonicalize()
    assert isinstance(remapping, MoleculeRemapping)
    assert remapping.to_correspondence().is_total()
    assert remapping.atoms.images == [3, 2, 1, 0]
    assert remapping.bonds.images == [0, 2, 1]


def test_reaction_span_canonicalize_error():
    span = ReactionSpan.from_entries(
        [
            (
                AtomForm(Element("C"), charge=NumForm.LitSet(set())),
                AtomForm(Element("C"), charge=NumForm.LitSet(set())),
            )
        ]
    )

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        span.canonicalize()
    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        span.tracked_canonicalize()


def test_reaction_to_reaction_span_error():
    reaction = Reaction(
        Molecule.from_entries(
            [AtomForm(Element("C")), AtomForm(Element("C"))],
            bonds=[(0, 1, BondForm(1))],
        ),
        Deltas(
            [
                Delta.Bond(
                    BondDelta.ModifyField(
                        id=0,
                        change=BondFieldChange.Order(
                            old=NumForm.Lit(2),
                            new=NumForm.Lit(3),
                        ),
                    )
                )
            ]
        ),
    )

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        reaction.to_reaction_span()

import pytest
from umol import AtomAst, BondAst, Element, ReactionSpanAst, ValueAst


def test_reaction_span_ast_from_entries():
    lhs = AtomAst(Element("C"), charge=ValueAst.Lit(1))
    rhs = AtomAst(Element("C"), charge=ValueAst.LitSet({1}))

    span = ReactionSpanAst.from_entries([(lhs, rhs)])

    assert span == ReactionSpanAst.from_entries([(lhs, lhs)])


def test_reaction_span_ast_from_entries_error():
    with pytest.raises(
        ValueError,
        match="^reaction span entry is absent from both sides$",
    ):
        ReactionSpanAst.from_entries([(None, None)])


def test_reaction_span_ast_from_entries_reference_error():
    with pytest.raises(
        ValueError,
        match="^reaction span entries reference unavailable atom 1$",
    ):
        ReactionSpanAst.from_entries(
            [(AtomAst(Element("C")), AtomAst(Element("C")))],
            bonds=[(0, 1, (BondAst(1), BondAst(1)))],
        )

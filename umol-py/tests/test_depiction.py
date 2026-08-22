import pytest

from umol import _native as _native_module

if not hasattr(_native_module, "Svg"):
    pytest.skip("umol-py was built without the depiction feature", allow_module_level=True)

from umol import (  # noqa: E402
    BondDelta,
    BondFieldChange,
    ContradictionError,
    Delta,
    Deltas,
    Molecule,
    MoleculeLayoutAlgorithm,
    NumForm,
    Reaction,
    Svg,
)


def test_molecule_layout_algorithm():
    algorithm = MoleculeLayoutAlgorithm.CoordGen()

    assert algorithm == MoleculeLayoutAlgorithm.CoordGen()
    assert repr(algorithm) == "MoleculeLayoutAlgorithm.CoordGen()"


def test_svg_constructor_error():
    with pytest.raises(TypeError):
        Svg()


def test_molecule_depict_with():
    svg = Molecule().depict_with(MoleculeLayoutAlgorithm.CoordGen())

    assert isinstance(svg, Svg)
    assert svg._repr_svg_() == (
        '<svg xmlns="http://www.w3.org/2000/svg" class="umol-depiction" '
        'viewBox="-0.5 -0.5 1 1">\n</svg>'
    )


def test_molecule_depict_with_algorithm_error():
    with pytest.raises(TypeError):
        Molecule().depict_with()


def test_reaction_depict_with():
    svg = Reaction.parse('{:lhs {:atoms ["C"]} :deltas []}').depict_with(
        MoleculeLayoutAlgorithm.CoordGen()
    )

    assert isinstance(svg, Svg)
    assert 'data-umol-item="arrow"' in svg._repr_svg_()
    assert "correspondence-pair/0" in svg._repr_svg_()


def test_reaction_depict_with_error():
    lhs = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')
    deltas = Deltas(
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
    )

    with pytest.raises(ContradictionError, match="^reached a contradiction$"):
        Reaction(lhs, deltas).depict_with(MoleculeLayoutAlgorithm.CoordGen())

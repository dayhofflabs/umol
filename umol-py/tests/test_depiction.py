import inspect
import sys
import xml.etree.ElementTree as ET

import pytest

from umol import _native as _native_module

if not hasattr(_native_module, "Depiction"):
    pytest.skip("umol-py was built without the depiction feature", allow_module_level=True)

from umol import (  # noqa: E402
    BondDelta,
    BondFieldChange,
    ContradictionError,
    Delta,
    DepictConfig,
    Depiction,
    Deltas,
    Molecule,
    MoleculeLayout,
    MoleculeLayoutAlgorithm,
    NumForm,
    Reaction,
    ReactionLayout,
)

TRANS_BUTENE = """{:atoms ["C" "C" "C" "C"]
    :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
    :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"""

REACTION = """{:lhs {:atoms ["C" "O"]
           :bonds [{:id :co :atoms [0 1] :attrs "1"}]}
     :deltas [{:bond {:modify [:co "2"]}}
              {:atom {:add [:n "N"]}}
              {:bond {:add [0 :n "1"]}}]}"""


def reflect(point, axis_start, axis_end):
    (px, py), (ax, ay), (bx, by) = point, axis_start, axis_end
    dx, dy = bx - ax, by - ay
    t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
    foot = (ax + t * dx, ay + t * dy)
    return (2 * foot[0] - px, 2 * foot[1] - py)


def test_depict_config_new():
    config = DepictConfig()

    assert config == DepictConfig.default()
    assert config.layout_algorithm == MoleculeLayoutAlgorithm.CoordGen()
    assert repr(config) == "DepictConfig.default()"


def test_depiction_constructor_error():
    with pytest.raises(TypeError):
        Depiction()


def test_molecule_layout_algorithm():
    algorithm = MoleculeLayoutAlgorithm.CoordGen()

    assert algorithm == MoleculeLayoutAlgorithm.CoordGen()
    assert repr(algorithm) == "MoleculeLayoutAlgorithm.CoordGen()"


@pytest.mark.parametrize(
    ("callable_", "expected"),
    [
        pytest.param(
            MoleculeLayout,
            "(positions)",
            marks=pytest.mark.skipif(
                sys.version_info < (3, 10),
                reason="class __text_signature__ requires Python 3.10",
            ),
        ),
        pytest.param(
            ReactionLayout,
            "(lhs, rhs, arrow_start, arrow_end)",
            marks=pytest.mark.skipif(
                sys.version_info < (3, 10),
                reason="class __text_signature__ requires Python 3.10",
            ),
        ),
        (MoleculeLayout.set_position, "(self, /, atom_id, position)"),
        (ReactionLayout.arrange, "(lhs, rhs)"),
        (ReactionLayout.set_lhs_position, "(self, /, atom_id, position)"),
        (ReactionLayout.set_rhs_position, "(self, /, atom_id, position)"),
        (ReactionLayout.set_arrow, "(self, /, start, end)"),
        (Molecule.layout, "(self, /, *, algorithm=Ellipsis)"),
        (Molecule.verify_layout, "(self, /, layout)"),
        (Molecule.depict_layout, "(self, /, layout)"),
        (Reaction.layout, "(self, /, *, algorithm=Ellipsis)"),
        (Reaction.verify_layout, "(self, /, layout)"),
        (Reaction.depict_layout, "(self, /, layout)"),
    ],
)
def test_layout_signature(callable_, expected):
    assert str(inspect.signature(callable_)) == expected


def test_molecule_layout_new():
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    assert layout.positions == ((0.0, 1.0), (2.0, -1.0))
    assert isinstance(layout.positions, tuple)
    assert len(layout) == 2
    assert layout == MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])
    assert layout != MoleculeLayout([(0.0, 1.0), (2.0, 1.0)])


def test_molecule_layout_unhashable():
    with pytest.raises(TypeError):
        hash(MoleculeLayout([(0.0, 0.0)]))


@pytest.mark.parametrize(
    ("positions", "message"),
    [
        pytest.param(
            [(0.0, 0.0), (float("nan"), 1.0)],
            "atom 1 has non-finite position",
            id="nan",
        ),
        pytest.param(
            [(float("inf"), 0.0)],
            "atom 0 has non-finite position",
            id="infinity",
        ),
    ],
)
def test_molecule_layout_new_error(positions, message):
    with pytest.raises(ValueError, match=message):
        MoleculeLayout(positions)


def test_molecule_layout_position():
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    assert layout.position(1) == (2.0, -1.0)
    with pytest.raises(IndexError, match="atom 2 is outside layout frame of size 2"):
        layout.position(2)


def test_molecule_layout_set_position():
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    assert layout.set_position(1, (3.5, 4.0)) is None

    assert layout == MoleculeLayout([(0.0, 1.0), (3.5, 4.0)])
    assert layout.position(1) == (3.5, 4.0)


@pytest.mark.parametrize(
    ("atom_id", "position", "message"),
    [
        pytest.param(2, (0.0, 0.0), "atom 2 is outside layout frame of size 2", id="frame"),
        pytest.param(0, (float("nan"), 0.0), "atom 0 has non-finite position", id="nan"),
        pytest.param(1, (0.0, float("-inf")), "atom 1 has non-finite position", id="infinity"),
    ],
)
def test_molecule_layout_set_position_error(atom_id, position, message):
    layout = MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])

    with pytest.raises(ValueError, match=message):
        layout.set_position(atom_id, position)

    assert layout == MoleculeLayout([(0.0, 1.0), (2.0, -1.0)])


@pytest.mark.parametrize(
    ("source", "item_kind", "reference", "svg_class"),
    [
        pytest.param(
            '{:atoms ["C#i13#c+#h2#u2" "O#c-"] :bonds [[0 1 "2"]]}',
            "atom",
            "molecule/atom/0",
            "umol-atom-right-subscript",
            id="atom-labels",
        ),
        pytest.param(
            """{:atoms ["C" "F" "Cl" "Br" "I"]
                 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                 :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}""",
            "wedge",
            "molecule/stereo-atom/0",
            "umol-wedge-solid",
            id="tetrahedral-stereo",
        ),
        pytest.param(
            """{:atoms ["C" "C" "C" "C" "C" "C"]
                 :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                         [3 4 "1"] [4 5 "1"] [5 0 "1"]]
                 :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*"}]}""",
            "dashed-contour",
            "molecule/aromatic-system/0",
            "umol-dashed-contour",
            id="aromatic-system",
        ),
    ],
)
def test_molecule_depict(source, item_kind, reference, svg_class):
    depiction = Molecule.parse(source).depict()
    text = depiction.render_svg()
    root = ET.fromstring(text)
    groups = list(root.findall("{http://www.w3.org/2000/svg}g"))

    assert isinstance(depiction, Depiction)
    assert depiction._repr_svg_() == text
    assert root.attrib["class"] == "umol-depiction"
    assert any(
        group.attrib.get("data-umol-item") == item_kind
        and reference in group.attrib.get("data-umol-references", "").split()
        for group in groups
    )
    assert f'class="{svg_class}"' in text
    assert 'data-umol-item="marker"' not in text


def test_molecule_depict_with():
    molecule = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "2"]]}')

    assert (
        molecule.depict_with(DepictConfig()).render_svg()
        == molecule.depict().render_svg()
    )


def test_molecule_depict_write_svg(tmp_path):
    depiction = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "2"]]}').depict()
    output = tmp_path / "molecule.svg"

    output.write_text(depiction.render_svg())

    assert output.read_text() == depiction.render_svg()


def test_molecule_depict_with_config_error():
    with pytest.raises(TypeError):
        Molecule().depict_with()


def test_molecule_layout():
    molecule = Molecule.parse('{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}')

    layout = molecule.layout()

    assert isinstance(layout, MoleculeLayout)
    assert len(layout) == 3
    assert layout == molecule.layout(algorithm=MoleculeLayoutAlgorithm.CoordGen())
    assert all(isinstance(coordinate, float) for x, y in layout.positions for coordinate in (x, y))


def test_molecule_layout_keyword_only_error():
    with pytest.raises(TypeError):
        Molecule().layout(MoleculeLayoutAlgorithm.CoordGen())


def test_molecule_depict_layout():
    molecule = Molecule.parse('{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}')

    assert (
        molecule.depict_layout(molecule.layout()).render_svg()
        == molecule.depict().render_svg()
    )


def bond_lines(text):
    root = ET.fromstring(text)
    lines = {}
    for group in root.findall("{http://www.w3.org/2000/svg}g"):
        if group.attrib.get("data-umol-item") != "bond":
            continue
        (line,) = list(group)
        lines[group.attrib["data-umol-references"]] = {
            name: float(line.attrib[name]) for name in ("x1", "y1", "x2", "y2")
        }
    return lines


def test_molecule_depict_layout_roundtrip():
    molecule = Molecule.parse(
        '{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]]}'
    )
    layout = molecule.layout()
    before = bond_lines(molecule.depict_layout(layout).render_svg())

    layout.set_position(1, (2.5, -1.25))
    after = bond_lines(molecule.depict_layout(layout).render_svg())

    assert (after["molecule/bond/0"]["x2"], after["molecule/bond/0"]["y2"]) == (2.5, 1.25)
    assert (after["molecule/bond/1"]["x1"], after["molecule/bond/1"]["y1"]) == (2.5, 1.25)
    assert (after["molecule/bond/0"]["x1"], after["molecule/bond/0"]["y1"]) == (
        before["molecule/bond/0"]["x1"],
        before["molecule/bond/0"]["y1"],
    )
    assert (after["molecule/bond/1"]["x2"], after["molecule/bond/1"]["y2"]) == (
        before["molecule/bond/1"]["x2"],
        before["molecule/bond/1"]["y2"],
    )
    assert after["molecule/bond/2"] == before["molecule/bond/2"]


def test_molecule_depict_layout_frame_error():
    molecule = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')

    with pytest.raises(
        ValueError,
        match="^layout frame: molecule atom count 2 does not match layout atom count 1$",
    ):
        molecule.depict_layout(MoleculeLayout([(0.0, 0.0)]))


def test_molecule_verify_layout():
    molecule = Molecule.from_smiles("C/C=C/C")

    assert molecule.verify_layout(molecule.layout()) is None


def test_molecule_verify_layout_cis_trans_error():
    molecule = Molecule.from_smiles("C/C=C/C")
    layout = molecule.layout()
    flipped = reflect(layout.position(3), layout.position(1), layout.position(2))
    layout.set_position(3, flipped)

    with pytest.raises(
        ValueError, match="^bond 1 is stored as E but its supplied coordinates draw Z$"
    ):
        molecule.verify_layout(layout)
    with pytest.raises(ValueError, match="^bond 1 is stored as E"):
        molecule.depict_layout(layout)


def test_molecule_verify_layout_cis_trans_degenerate_error():
    molecule = Molecule.parse(TRANS_BUTENE)

    with pytest.raises(ValueError, match="^bond 1: ligand atom 0 lies on the cis/trans site axis$"):
        molecule.verify_layout(MoleculeLayout([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 1.0)]))


@pytest.mark.parametrize(
    ("positions", "message"),
    [
        pytest.param(
            [(0.0, 0.0), (1.0, 0.0), (1.0, 0.0)],
            "^bond 1 has non-finite or degenerate derived geometry$",
            id="coincident-bonded-atoms",
        ),
        pytest.param(
            [(1e308, 0.0), (-1e308, 0.0), (0.0, 1.0)],
            "^bond 0 has non-finite or degenerate derived geometry$",
            id="overflowing-bond",
        ),
    ],
)
def test_molecule_verify_layout_geometry_error(positions, message):
    molecule = Molecule.parse('{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}')

    with pytest.raises(ValueError, match=message):
        molecule.verify_layout(MoleculeLayout(positions))


def test_molecule_verify_layout_frame_error():
    molecule = Molecule.parse('{:atoms ["C" "O"] :bonds [[0 1 "1"]]}')

    with pytest.raises(
        ValueError,
        match="^layout frame: molecule atom count 2 does not match layout atom count 3$",
    ):
        molecule.verify_layout(MoleculeLayout([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]))


def test_molecule_depict_layout_tetrahedral_error():
    molecule = Molecule.parse(
        """{:atoms ["C" "F" "Cl" "Br" "I"]
             :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
             :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"""
    )
    collinear = MoleculeLayout([(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)])

    with pytest.raises(
        ValueError,
        match="^tetrahedral geometry cannot establish a display wedge for stereo atom 0$",
    ):
        molecule.depict_layout(collinear)


def test_reaction_layout_new():
    lhs = MoleculeLayout([(0.0, 0.0), (1.0, 0.0)])
    rhs = MoleculeLayout([(4.0, 0.0)])

    layout = ReactionLayout(lhs, rhs, (-1, 0), (1, 0))

    assert layout.lhs == lhs
    assert layout.rhs == rhs
    assert layout.arrow_start == (-1.0, 0.0)
    assert layout.arrow_end == (1.0, 0.0)
    assert layout == ReactionLayout(lhs, rhs, (-1.0, 0.0), (1.0, 0.0))
    assert layout != ReactionLayout(lhs, rhs, (-1.0, 0.0), (2.0, 0.0))


@pytest.mark.parametrize(
    ("start", "end", "message"),
    [
        pytest.param((1.0, 1.0), (1.0, 1.0), "reaction arrow starts and ends at", id="degenerate"),
        pytest.param((0.0, 0.0), (float("nan"), 0.0), "has a non-finite endpoint", id="nan"),
    ],
)
def test_reaction_layout_new_error(start, end, message):
    side = MoleculeLayout([(0.0, 0.0)])

    with pytest.raises(ValueError, match=message):
        ReactionLayout(side, side, start, end)


def test_reaction_layout_unhashable():
    side = MoleculeLayout([(0.0, 0.0)])

    with pytest.raises(TypeError):
        hash(ReactionLayout(side, side, (-1.0, 0.0), (1.0, 0.0)))


def test_reaction_layout_arrange():
    lhs = MoleculeLayout([(0.0, 0.0), (1.0, 0.0)])
    rhs = MoleculeLayout([(4.0, 0.0)])

    layout = ReactionLayout.arrange(lhs, rhs)

    assert isinstance(layout, ReactionLayout)
    assert layout.lhs == MoleculeLayout([(-2.75, 0.0), (-1.75, 0.0)])
    assert layout.rhs == MoleculeLayout([(1.75, 0.0)])
    assert layout.arrow_start == (-0.75, 0.0)
    assert layout.arrow_end == (0.75, 0.0)


def test_reaction_layout_arrange_error():
    with pytest.raises(ValueError, match="^lhs translation: atom 0 has non-finite position"):
        ReactionLayout.arrange(
            MoleculeLayout([(-1e308, 0.0), (1e308, 0.0)]), MoleculeLayout([(0.0, 0.0)])
        )


def test_reaction_layout_side_copies():
    side = MoleculeLayout([(0.0, 0.0), (1.0, 0.0)])
    layout = ReactionLayout(side, side, (-1.0, 0.0), (1.0, 0.0))

    layout.lhs.set_position(0, (5.0, 5.0))

    assert layout.lhs == side

    assert layout.set_lhs_position(0, (5.0, 5.0)) is None
    assert layout.set_rhs_position(1, (6.0, 6.0)) is None

    assert layout.lhs == MoleculeLayout([(5.0, 5.0), (1.0, 0.0)])
    assert layout.rhs == MoleculeLayout([(0.0, 0.0), (6.0, 6.0)])


@pytest.mark.parametrize(
    ("atom_id", "position", "message"),
    [
        pytest.param(2, (0.0, 0.0), "atom 2 is outside layout frame of size 2", id="frame"),
        pytest.param(0, (float("inf"), 0.0), "atom 0 has non-finite position", id="infinity"),
    ],
)
def test_reaction_layout_set_side_position_error(atom_id, position, message):
    side = MoleculeLayout([(0.0, 0.0), (1.0, 0.0)])
    layout = ReactionLayout(side, side, (-1.0, 0.0), (1.0, 0.0))

    with pytest.raises(ValueError, match=message):
        layout.set_lhs_position(atom_id, position)
    with pytest.raises(ValueError, match=message):
        layout.set_rhs_position(atom_id, position)

    assert layout == ReactionLayout(side, side, (-1.0, 0.0), (1.0, 0.0))


def test_reaction_layout_set_arrow():
    side = MoleculeLayout([(0.0, 0.0)])
    layout = ReactionLayout(side, side, (-1.0, 0.0), (1.0, 0.0))

    assert layout.set_arrow((0.0, -2.0), (0.0, 2.0)) is None

    assert layout.arrow_start == (0.0, -2.0)
    assert layout.arrow_end == (0.0, 2.0)


@pytest.mark.parametrize(
    ("start", "end", "message"),
    [
        pytest.param((3.0, 3.0), (3.0, 3.0), "reaction arrow starts and ends at", id="degenerate"),
        pytest.param((float("nan"), 0.0), (1.0, 0.0), "has a non-finite endpoint", id="nan"),
    ],
)
def test_reaction_layout_set_arrow_error(start, end, message):
    side = MoleculeLayout([(0.0, 0.0)])
    layout = ReactionLayout(side, side, (-1.0, 0.0), (1.0, 0.0))

    with pytest.raises(ValueError, match=message):
        layout.set_arrow(start, end)

    assert layout.arrow_start == (-1.0, 0.0)
    assert layout.arrow_end == (1.0, 0.0)


def test_reaction_layout():
    reaction = Reaction.parse(REACTION)

    layout = reaction.layout()

    assert isinstance(layout, ReactionLayout)
    assert layout == reaction.layout(algorithm=MoleculeLayoutAlgorithm.CoordGen())
    assert len(layout.lhs) == 2
    assert len(layout.rhs) == 3
    assert layout.arrow_start == (-0.75, 0.0)
    assert layout.arrow_end == (0.75, 0.0)


def test_reaction_layout_keyword_only_error():
    with pytest.raises(TypeError):
        Reaction.parse(REACTION).layout(MoleculeLayoutAlgorithm.CoordGen())


def test_reaction_verify_layout():
    reaction = Reaction.parse(REACTION)

    assert reaction.verify_layout(reaction.layout()) is None


def test_reaction_depict_layout():
    reaction = Reaction.parse(REACTION)

    assert reaction.depict_layout(reaction.layout()).render_svg() == reaction.depict().render_svg()


def test_reaction_depict_layout_frame_error():
    reaction = Reaction.parse(REACTION)
    layout = reaction.layout()
    mismatched = ReactionLayout(
        layout.lhs, MoleculeLayout([(2.0, 0.0)]), layout.arrow_start, layout.arrow_end
    )

    with pytest.raises(
        ValueError,
        match="^rhs depiction: layout frame: "
        "molecule atom count 3 does not match layout atom count 1$",
    ):
        reaction.depict_layout(mismatched)
    with pytest.raises(ValueError, match="^rhs depiction: layout frame"):
        reaction.verify_layout(mismatched)


def test_reaction_depict():
    depiction = Reaction.parse(REACTION).depict()
    text = depiction.render_svg()
    root = ET.fromstring(text)
    groups = list(root.findall("{http://www.w3.org/2000/svg}g"))
    mapping_groups = [
        group
        for group in groups
        if "correspondence-pair/" in group.attrib.get("data-umol-references", "")
    ]
    arrow = next(group for group in groups if group.attrib.get("data-umol-item") == "arrow")
    shaft, head = list(arrow)

    assert isinstance(depiction, Depiction)
    assert depiction._repr_svg_() == text
    assert [group.attrib["data-umol-references"] for group in mapping_groups] == [
        "reaction-lhs/atom/0 correspondence-pair/0",
        "reaction-lhs/atom/1 correspondence-pair/1",
        "reaction-rhs/atom/0 correspondence-pair/0",
        "reaction-rhs/atom/1 correspondence-pair/1",
    ]
    assert [list(group)[0].attrib["font-size"] for group in mapping_groups] == [
        "0.3825",
        "0.3825",
        "0.3825",
        "0.3825",
    ]
    assert shaft.attrib["class"] == "umol-arrow-shaft"
    assert shaft.attrib["x1"] == "-0.75"
    assert shaft.attrib["x2"] == "0.51"
    assert head.attrib == {
        "class": "umol-arrow-head",
        "points": "0.75,0 0.51,-0.11 0.51,0.11",
        "fill": "currentColor",
    }


def test_reaction_depict_with():
    reaction = Reaction.parse(
        """{:lhs {:atoms ["C" "O"] :bonds [[0 1 "1"]]}
             :deltas []}"""
    )

    assert (
        reaction.depict_with(DepictConfig()).render_svg()
        == reaction.depict().render_svg()
    )


@pytest.mark.parametrize(
    "operation",
    [
        pytest.param(lambda reaction: reaction.depict_with(DepictConfig()), id="depict_with"),
        pytest.param(lambda reaction: reaction.layout(), id="layout"),
    ],
)
def test_reaction_contradiction_error(operation):
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
        operation(Reaction(lhs, deltas))

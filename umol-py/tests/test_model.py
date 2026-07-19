import re

import pytest

from umol import AtomAst, AtomTypeRegistry, Element


def test_atom_type_registry_default():
    registry = AtomTypeRegistry.default()

    assert registry == AtomTypeRegistry.default()
    assert repr(registry) == "AtomTypeRegistry.default()"
    assert registry.content_hash_hex == f"{registry.content_hash:016x}"


def test_atom_type_registry_from_atoms():
    atom = AtomAst.parse("C#c0#v4")
    registry = AtomTypeRegistry.from_atoms(
        [atom, AtomAst.parse("C#c+#v3"), AtomAst.parse("O#c0#v2")]
    )

    assert registry == AtomTypeRegistry.from_atoms(
        [
            AtomAst.parse("C#c0#v4"),
            AtomAst.parse("C#c+#v3"),
            AtomAst.parse("O#c0#v2"),
        ]
    )
    assert registry != AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4")])

    atom.charge = 1
    assert registry.patterns_for_element_and_charge(Element("C"), 0) == [
        AtomAst.parse("C#c0#v4")
    ]


@pytest.mark.parametrize(
    ("atom", "message"),
    [
        (
            AtomAst.parse("*#c0"),
            "atom type registry entry 0 must have a literal element",
        ),
        (
            AtomAst.parse("C"),
            "atom type registry entry 0 must have a literal charge",
        ),
        (
            AtomAst.parse("C#c128"),
            "atom type registry entry 0 charge 128 is outside -128..=127",
        ),
    ],
)
def test_atom_type_registry_from_atoms_error(atom, message):
    with pytest.raises(ValueError, match=f"^{re.escape(message)}$"):
        AtomTypeRegistry.from_atoms([atom])


def test_atom_type_registry_from_toml():
    registry = AtomTypeRegistry.from_toml(
        '[C]\n0 = ["C#c0#v4"]\n1 = ["C#c+#v3"]\n[O]\n0 = ["O#c0#v2"]'
    )

    assert registry.patterns_for_element(Element("C")) == [
        AtomAst.parse("C#i=#c0#h0#n0#u0#s#v4#d0#t0#a!#m!"),
        AtomAst.parse("C#i=#c+#h0#n0#u0#s#v3#d0#t0#a!#m!"),
    ]
    assert registry.patterns_for_element_and_charge(Element("O"), 0) == [
        AtomAst.parse("O#i=#c0#h0#n0#u0#s#v2#d0#t0#a!#m!")
    ]


def test_atom_type_registry_from_toml_error():
    with pytest.raises(
        ValueError,
        match="^invalid atom type registry: unknown element: X$",
    ):
        AtomTypeRegistry.from_toml('[X]\n0 = ["X#c0"]')


def test_atom_type_registry_patterns_are_detached():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("C#c+#v3")]
    )
    patterns = registry.patterns_for_element(Element("C"))

    patterns[0].charge = -1
    patterns.pop()

    assert registry.patterns_for_element(Element("C")) == [
        AtomAst.parse("C#c0#v4"),
        AtomAst.parse("C#c+#v3"),
    ]


def test_atom_type_registry_repr():
    registry = AtomTypeRegistry.from_atoms(
        [AtomAst.parse("C#c0#v4"), AtomAst.parse("O#c0#v2")]
    )

    assert repr(registry) == (
        'AtomTypeRegistry.from_atoms([AtomAst.parse("C#c0#v4"), '
        'AtomAst.parse("O#c0#v2")])'
    )


def test_atom_type_registry_immutable():
    registry = AtomTypeRegistry.default()

    with pytest.raises(AttributeError):
        registry.content_hash = 0
    assert not hasattr(registry, "add")

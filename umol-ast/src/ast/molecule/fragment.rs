//! Fragment algebra — a `MoleculeAst` body with a positional, colour-typed **port interface**.
//! `attach` (operadic composition) wires two fragments through a pair of ports; `+` (monoidal
//! product) juxtaposes them without wiring; `finish` finalizes to the body (`finish_open` closes to a
//! pattern, capping each free port with a wildcard atom).

use std::ops::Add;

use super::super::atom::{AtomAst, ElementAst};
use super::super::bond::BondAst;
use super::super::correspondence::MoleculeCorrespondence;
use super::super::id::AtomId;
use super::super::traits::Lattice;
use super::MoleculeAst;
#[cfg(test)]
use super::MoleculeParts;

/// An attachment point on a fragment: which body `atom` exposes the free valence, the `bond` spec
/// (the port's *colour*) formed when it attaches, and an optional `name` to address it. Ports are
/// ordered; the name is a label, not part of the typing — compatibility on `attach` is the `meet` of
/// the two ports' bonds.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Port {
    pub atom: AtomId,
    pub bond: BondAst,
    pub name: Option<String>,
}

/// How a port is addressed on `attach` — by position or by name. `From<u32>`/`From<i32>` → index,
/// `From<&str>`/`From<String>` → name (mirroring `AtomArg`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PortArg {
    Index(u32),
    Name(String),
}

impl From<u32> for PortArg {
    fn from(index: u32) -> Self {
        Self::Index(index)
    }
}

impl From<i32> for PortArg {
    fn from(index: i32) -> Self {
        Self::Index(
            u32::try_from(index)
                .unwrap_or_else(|_| panic!("port index must be non-negative, got {index}")),
        )
    }
}

impl From<&str> for PortArg {
    fn from(name: &str) -> Self {
        Self::Name(name.to_string())
    }
}

impl From<String> for PortArg {
    fn from(name: String) -> Self {
        Self::Name(name)
    }
}

/// A subgraph with a port interface: a `MoleculeAst` body plus the ordered ports it may attach
/// through. Compose with `attach` / `+`, finalize with `finish` (or `finish_open` for a pattern).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fragment {
    body: MoleculeAst,
    ports: Vec<Port>,
}

impl Fragment {
    /// A fragment over `body` with no ports.
    pub fn new(body: MoleculeAst) -> Self {
        Self {
            body,
            ports: Vec::new(),
        }
    }

    /// Declare a named port on `atom`, forming `bond` when attached.
    pub fn with_port(
        mut self,
        name: impl Into<String>,
        atom: AtomId,
        bond: impl Into<BondAst>,
    ) -> Self {
        self.ports.push(Port {
            atom,
            bond: bond.into(),
            name: Some(name.into()),
        });
        self
    }

    /// Declare an unnamed port on `atom` (addressed by index only).
    pub fn with_unnamed_port(mut self, atom: AtomId, bond: impl Into<BondAst>) -> Self {
        self.ports.push(Port {
            atom,
            bond: bond.into(),
            name: None,
        });
        self
    }

    /// The body, for inspection.
    pub fn body(&self) -> &MoleculeAst {
        &self.body
    }

    /// The ports, in order.
    pub fn ports(&self) -> &[Port] {
        &self.ports
    }

    /// Finalize into the body. Every port must have been paired by `attach`; a remaining free port
    /// panics — a construction bug, like a bad port ref. For a pattern that leaves ports open, use
    /// [`finish_open`](Self::finish_open).
    pub fn finish(self) -> MoleculeAst {
        assert!(
            self.ports.is_empty(),
            "fragment has {} unpaired port(s); attach them or use `finish_open`",
            self.ports.len()
        );
        self.body
    }

    /// Finalize into a pattern: each free port becomes an undetermined wildcard atom bonded to the
    /// port's atom through the port's colour — "something attaches here".
    pub fn finish_open(self) -> MoleculeAst {
        let mut editor = self.body.edit();
        for port in self.ports {
            let wildcard = editor.add_atom(AtomAst::new(ElementAst::undetermined()));
            editor.add_bond(port.atom, wildcard, port.bond);
        }
        editor.build()
    }

    /// Attach `self`'s `self_port` to `other`'s `other_port` — the operadic composition (the only
    /// operation that wires). Joins the bodies, bonds the two port atoms with the `meet` of their
    /// bond specs, and drops the two consumed ports; remaining ports carry forward (`other`'s
    /// remapped). Panics on an unknown/ambiguous port ref or a ⊥-`meet` (incompatible ports).
    pub fn attach(
        self,
        self_port: impl Into<PortArg>,
        other: Fragment,
        other_port: impl Into<PortArg>,
    ) -> Fragment {
        let self_index = resolve_port(&self.ports, self_port.into());
        let other_index = resolve_port(&other.ports, other_port.into());

        let self_atom = self.ports[self_index].atom;
        let other_port_atom = other.ports[other_index].atom;
        let bond = self.ports[self_index]
            .bond
            .meet(&other.ports[other_index].bond)
            .unwrap_or_else(|| {
                panic!(
                    "incompatible ports: {:?} and {:?} have no common bond",
                    self.ports[self_index].bond, other.ports[other_index].bond
                )
            });

        let (body, correspondence) = self.body.combine(&other.body);
        let other_atom = correspondence
            .atoms()
            .right_of(other_port_atom)
            .expect("combine maps every atom of `other`");
        let mut editor = body.edit();
        editor.add_bond(self_atom, other_atom, bond);
        let body = editor.build();

        let mut ports: Vec<Port> = self
            .ports
            .into_iter()
            .enumerate()
            .filter(|(index, _)| *index != self_index)
            .map(|(_, port)| port)
            .collect();
        ports.extend(
            other
                .ports
                .into_iter()
                .enumerate()
                .filter(|(index, _)| *index != other_index)
                .map(|(_, port)| remap_port(port, &correspondence)),
        );
        Fragment { body, ports }
    }
}

/// Juxtapose two fragments — the monoidal product. Combines the bodies (no bond formed) and
/// concatenates the ports, `other`'s remapped through the combination correspondence.
impl Add<Fragment> for Fragment {
    type Output = Fragment;
    fn add(self, other: Fragment) -> Fragment {
        let (body, correspondence) = self.body.combine(&other.body);
        let mut ports = self.ports;
        ports.extend(
            other
                .ports
                .into_iter()
                .map(|port| remap_port(port, &correspondence)),
        );
        Fragment { body, ports }
    }
}

/// Move a port's atom from `other`'s id space into the combined body's, through the `other → union`
/// correspondence a `combine` returns (`self` is the prefix, so its ports are unchanged).
fn remap_port(port: Port, correspondence: &MoleculeCorrespondence) -> Port {
    let atom = correspondence
        .atoms()
        .right_of(port.atom)
        .expect("combine maps every atom of `other`");
    Port { atom, ..port }
}

/// Resolve a `PortArg` to an index into `ports`. Panics if the index is out of range, or a name is
/// unknown or ambiguous.
fn resolve_port(ports: &[Port], port_ref: PortArg) -> usize {
    match port_ref {
        PortArg::Index(index) => {
            let index = index as usize;
            assert!(
                index < ports.len(),
                "port index {index} out of range ({} ports)",
                ports.len()
            );
            index
        }
        PortArg::Name(name) => {
            let mut matches = ports
                .iter()
                .enumerate()
                .filter(|(_, port)| port.name.as_deref() == Some(name.as_str()));
            let (index, _) = matches
                .next()
                .unwrap_or_else(|| panic!("no port named {name:?}"));
            assert!(
                matches.next().is_none(),
                "port name {name:?} is ambiguous; use an index"
            );
            index
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::atom::{AtomAst, ElementAst};
    use crate::ast::id::BondId;

    #[rstest]
    fn test_fragment_with_port() {
        let body = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        });
        let fragment = Fragment::new(body)
            .with_port("open", AtomId(0), BondAst::from_order(1))
            .with_unnamed_port(AtomId(0), BondAst::from_order(2));

        assert_eq!(
            fragment.ports(),
            &[
                Port {
                    atom: AtomId(0),
                    bond: BondAst::from_order(1),
                    name: Some("open".to_string()),
                },
                Port {
                    atom: AtomId(0),
                    bond: BondAst::from_order(2),
                    name: None,
                },
            ]
        );
    }

    #[rstest]
    fn test_fragment_finish() {
        let body = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        });

        assert_eq!(Fragment::new(body.clone()).finish(), body);
    }

    #[rstest]
    #[should_panic(expected = "unpaired port")]
    fn test_fragment_finish_error() {
        let body = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        });
        Fragment::new(body)
            .with_port("open", AtomId(0), BondAst::from_order(1))
            .finish();
    }

    #[rstest]
    fn test_fragment_finish_open() {
        let body = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        });
        let pattern = Fragment::new(body)
            .with_port("r", AtomId(0), BondAst::from_order(2))
            .finish_open();

        // the free port became a wildcard atom double-bonded to the carbon
        assert_eq!(pattern.atoms().count(), 2);
        assert_eq!(
            pattern.atom(AtomId(1)).ast,
            &AtomAst::new(ElementAst::undetermined())
        );
        assert_eq!(pattern.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
        assert_eq!(pattern.bond(BondId(0)).ast, &BondAst::from_order(2));
    }

    #[rstest]
    #[case::index(PortArg::from(2_u32), PortArg::Index(2))]
    #[case::signed_index(PortArg::from(2_i32), PortArg::Index(2))]
    #[case::name(PortArg::from("left"), PortArg::Name("left".to_string()))]
    #[case::owned_name(PortArg::from("left".to_string()), PortArg::Name("left".to_string()))]
    fn test_port_arg_from(#[case] port_arg: PortArg, #[case] expected: PortArg) {
        assert_eq!(port_arg, expected);
    }

    #[rstest]
    fn test_fragment_add() {
        let left = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("left", AtomId(0), BondAst::from_order(1));
        let right = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::from_element(Element::O),
                AtomAst::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ..Default::default()
        }))
        .with_port("right", AtomId(1), BondAst::from_order(2));

        let combined = left + right;

        // Bodies are combined disjointly; no bond is formed between them.
        assert_eq!(combined.body().atoms().count(), 3);
        assert_eq!(combined.body().bonds().count(), 1);
        // left's port unchanged (prefix); right's port atom 1 remapped to 2
        assert_eq!(
            combined.ports(),
            &[
                Port {
                    atom: AtomId(0),
                    bond: BondAst::from_order(1),
                    name: Some("left".to_string()),
                },
                Port {
                    atom: AtomId(2),
                    bond: BondAst::from_order(2),
                    name: Some("right".to_string()),
                },
            ]
        );
    }

    #[rstest]
    fn test_fragment_attach() {
        let left = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("a", AtomId(0), BondAst::from_order(1));
        let right = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("b", AtomId(0), BondAst::from_order(1));

        let joined = left.attach("a", right, "b");

        assert_eq!(joined.body().atoms().count(), 2);
        assert_eq!(joined.body().bonds().count(), 1);
        assert_eq!(joined.body().bond(BondId(0)).ast, &BondAst::from_order(1));
        assert_eq!(
            joined.body().bond(BondId(0)).atom_ids(),
            [AtomId(0), AtomId(1)]
        );
        assert!(joined.ports().is_empty());
    }

    #[rstest]
    fn test_fragment_attach_meet() {
        // an undetermined (⊤) port absorbs the partner's bond spec through `meet`
        let left = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("a", AtomId(0), BondAst::default());
        let right = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("b", AtomId(0), BondAst::from_order(2));

        let body = left.attach("a", right, "b").finish();

        assert_eq!(body.bond(BondId(0)).ast, &BondAst::from_order(2));
    }

    #[rstest]
    #[should_panic(expected = "no port named")]
    fn test_fragment_attach_error_missing_port() {
        let left = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("a", AtomId(0), BondAst::from_order(1));
        let right = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("b", AtomId(0), BondAst::from_order(1));

        left.attach("a", right, "missing");
    }

    #[rstest]
    #[should_panic(expected = "incompatible ports")]
    fn test_fragment_attach_error_incompatible() {
        let left = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("a", AtomId(0), BondAst::from_order(1));
        let right = Fragment::new(MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::O)],
            bonds: Vec::new(),
            ..Default::default()
        }))
        .with_port("b", AtomId(0), BondAst::from_order(2));

        left.attach("a", right, "b");
    }
}

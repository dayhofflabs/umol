//! Python mirrors for resolved molecule deltas and their field-change payloads.

use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticSystemFieldChange as AstAromaticSystemFieldChange,
    AtomFieldChange as AstAtomFieldChange, BondFieldChange as AstBondFieldChange,
    DativeBondFieldChange as AstDativeBondFieldChange, ElectronCountsAst as AstElectronCountsAst,
    ElementAst as AstElementAst, IsotopeMassAst as AstIsotopeMassAst,
    MulticenterBondFieldChange as AstMulticenterBondFieldChange,
    NoncovalentBondFieldChange as AstNoncovalentBondFieldChange,
    NoncovalentBondKindAst as AstNoncovalentBondKindAst, SpinStateAst as AstSpinStateAst,
    StereoAtomFieldChange as AstStereoAtomFieldChange,
    StereoBondFieldChange as AstStereoBondFieldChange,
    StereoConfigurationAst as AstStereoConfigurationAst, ValueAst as AstValueAst,
};

use crate::atom::{ElementAst, IsotopeMassAst};
use crate::convert::into_py_variant;
use crate::electrons::ElectronCountsAst;
use crate::noncovalent::NoncovalentBondKindAst;
use crate::spin::SpinStateAst;
use crate::stereo::StereoConfigurationAst;
use crate::value::ValueAst;

/// Conversion shared by the bound value types used as old/new field payloads.
trait FieldValueMirror: Sized {
    type Ast;

    fn from_ast(py: Python<'_>, ast: &Self::Ast) -> PyResult<Self>;
    fn to_ast(&self, py: Python<'_>) -> Self::Ast;
}

impl FieldValueMirror for ElementAst {
    type Ast = AstElementAst;

    fn from_ast(_py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Ok(Self::from_ast(ast))
    }

    fn to_ast(&self, _py: Python<'_>) -> Self::Ast {
        self.to_ast()
    }
}

impl FieldValueMirror for IsotopeMassAst {
    type Ast = AstIsotopeMassAst;

    fn from_ast(_py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Ok(Self::from_ast(ast))
    }

    fn to_ast(&self, _py: Python<'_>) -> Self::Ast {
        self.to_ast()
    }
}

impl FieldValueMirror for ValueAst {
    type Ast = AstValueAst;

    fn from_ast(py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Self::from_ast(py, ast)
    }

    fn to_ast(&self, py: Python<'_>) -> Self::Ast {
        self.to_ast(py)
    }
}

impl FieldValueMirror for SpinStateAst {
    type Ast = AstSpinStateAst;

    fn from_ast(py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Self::from_ast(py, ast)
    }

    fn to_ast(&self, py: Python<'_>) -> Self::Ast {
        self.to_ast(py)
    }
}

impl FieldValueMirror for ElectronCountsAst {
    type Ast = AstElectronCountsAst;

    fn from_ast(_py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Ok(Self::from_ast(ast))
    }

    fn to_ast(&self, _py: Python<'_>) -> Self::Ast {
        self.to_ast()
    }
}

impl FieldValueMirror for NoncovalentBondKindAst {
    type Ast = AstNoncovalentBondKindAst;

    fn from_ast(_py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Ok(Self::from_ast(ast))
    }

    fn to_ast(&self, _py: Python<'_>) -> Self::Ast {
        self.to_ast()
    }
}

impl FieldValueMirror for StereoConfigurationAst {
    type Ast = AstStereoConfigurationAst;

    fn from_ast(py: Python<'_>, ast: &Self::Ast) -> PyResult<Self> {
        Self::from_ast(py, ast)
    }

    fn to_ast(&self, py: Python<'_>) -> Self::Ast {
        self.to_ast(py)
    }
}

/// Render a named old/new complex-enum variant using the child objects' reprs.
fn field_change_repr(obj: &Bound<'_, PyAny>, type_name: &str, variant: &str) -> PyResult<String> {
    let old = obj.getattr("old")?.repr()?.extract::<String>()?;
    let new = obj.getattr("new")?.repr()?.extract::<String>()?;
    Ok(format!("{type_name}.{variant}(old={old}, new={new})"))
}

macro_rules! field_change {
    (
        $(#[$meta:meta])*
        $name:ident => $ast:ident {
            $(
                $variant:ident($value:ty)
            ),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[pyclass]
        pub enum $name {
            $(
                $variant {
                    old: Py<$value>,
                    new: Py<$value>,
                },
            )+
        }

        #[pymethods]
        impl $name {
            fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
                self.to_ast(py) == other.to_ast(py)
            }

            fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
                let variant = match &*slf.bind(py).borrow() {
                    $(Self::$variant { .. } => stringify!($variant),)+
                };
                field_change_repr(
                    slf.bind(py).as_any(),
                    stringify!($name),
                    variant,
                )
            }

            /// Return the same field change with its old and new values exchanged.
            fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>> {
                into_py_variant(py, Self::from_ast(py, &self.to_ast(py).inverse())?)
            }
        }

        impl $name {
            pub(crate) fn from_ast(py: Python<'_>, change: &$ast) -> PyResult<Self> {
                Ok(match change {
                    $(
                        $ast::$variant { old, new } => Self::$variant {
                            old: into_py_variant(
                                py,
                                <$value as FieldValueMirror>::from_ast(py, old)?,
                            )?,
                            new: into_py_variant(
                                py,
                                <$value as FieldValueMirror>::from_ast(py, new)?,
                            )?,
                        },
                    )+
                })
            }

            pub(crate) fn to_ast(&self, py: Python<'_>) -> $ast {
                match self {
                    $(
                        Self::$variant { old, new } => $ast::$variant {
                            old: <$value as FieldValueMirror>::to_ast(
                                &old.bind(py).borrow(),
                                py,
                            ),
                            new: <$value as FieldValueMirror>::to_ast(
                                &new.bind(py).borrow(),
                                py,
                            ),
                        },
                    )+
                }
            }
        }
    };
}

field_change! {
    /// An atom attribute change carrying the field's old and new AST values.
    AtomFieldChange => AstAtomFieldChange {
        Element(ElementAst),
        IsotopeMass(IsotopeMassAst),
        Charge(ValueAst),
        ImplicitHydrogens(ValueAst),
        LonePairs(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A covalent-bond attribute change carrying the field's old and new AST values.
    BondFieldChange => AstBondFieldChange {
        Order(ValueAst),
        Charge(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A dative-bond attribute change carrying the field's old and new AST values.
    DativeBondFieldChange => AstDativeBondFieldChange {
        Order(ValueAst),
    }
}

field_change! {
    /// An aromatic-system attribute change carrying the field's old and new AST values.
    AromaticSystemFieldChange => AstAromaticSystemFieldChange {
        Electrons(ElectronCountsAst),
        Charge(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A multicenter-bond attribute change carrying the field's old and new AST values.
    MulticenterBondFieldChange => AstMulticenterBondFieldChange {
        Electrons(ElectronCountsAst),
        Charge(ValueAst),
        Spin(SpinStateAst),
    }
}

field_change! {
    /// A noncovalent-bond kind change carrying the field's old and new AST values.
    NoncovalentBondFieldChange => AstNoncovalentBondFieldChange {
        Kind(NoncovalentBondKindAst),
    }
}

field_change! {
    /// A stereo-atom configuration change carrying the field's old and new AST values.
    StereoAtomFieldChange => AstStereoAtomFieldChange {
        Configuration(StereoConfigurationAst),
    }
}

field_change! {
    /// A stereo-bond configuration change carrying the field's old and new AST values.
    StereoBondFieldChange => AstStereoBondFieldChange {
        Configuration(StereoConfigurationAst),
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        NoncovalentBondKind as AstNoncovalentBondKind, SpinStateAst as AstSpinStateAst,
        StereoCosetAst as AstStereoCosetAst, StereoKind as AstStereoKind, ValueAst as AstValueAst,
    };
    use umol_chem::element::Element;

    use super::*;

    #[rstest]
    #[case::element(AstAtomFieldChange::Element {
        old: AstElementAst::Lit(Element::C),
        new: AstElementAst::Lit(Element::N),
    })]
    #[case::isotope_mass(AstAtomFieldChange::IsotopeMass {
        old: AstIsotopeMassAst::Lit(12),
        new: AstIsotopeMassAst::Lit(13),
    })]
    #[case::charge(AstAtomFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::implicit_hydrogens(AstAtomFieldChange::ImplicitHydrogens {
        old: AstValueAst::Lit(3),
        new: AstValueAst::Lit(2),
    })]
    #[case::lone_pairs(AstAtomFieldChange::LonePairs {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::spin(AstAtomFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_atom_field_change_roundtrip(#[case] change: AstAtomFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                AtomFieldChange::from_ast(py, &change).unwrap().to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        true,
    )]
    #[case::different(
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(1),
        },
        false,
    )]
    fn test_atom_field_change_eq(
        #[case] lhs: AstAtomFieldChange,
        #[case] rhs: AstAtomFieldChange,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = AtomFieldChange::from_ast(py, &lhs).unwrap();
            let rhs = AtomFieldChange::from_ast(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::element(
        AstAtomFieldChange::Element {
            old: AstElementAst::Lit(Element::C),
            new: AstElementAst::Lit(Element::N),
        },
        "ElementAst.Lit(Element('C'))",
        "ElementAst.Lit(Element('N'))",
        "AtomFieldChange.Element(old=ElementAst.Lit(Element('C')), new=ElementAst.Lit(Element('N')))"
    )]
    #[case::isotope_mass(
        AstAtomFieldChange::IsotopeMass {
            old: AstIsotopeMassAst::Lit(12),
            new: AstIsotopeMassAst::Lit(13),
        },
        "IsotopeMassAst.Lit(12)",
        "IsotopeMassAst.Lit(13)",
        "AtomFieldChange.IsotopeMass(old=IsotopeMassAst.Lit(12), new=IsotopeMassAst.Lit(13))"
    )]
    #[case::charge(
        AstAtomFieldChange::Charge {
            old: AstValueAst::Lit(0),
            new: AstValueAst::Lit(-1),
        },
        "ValueAst.Lit(0)",
        "ValueAst.Lit(-1)",
        "AtomFieldChange.Charge(old=ValueAst.Lit(0), new=ValueAst.Lit(-1))"
    )]
    #[case::implicit_hydrogens(
        AstAtomFieldChange::ImplicitHydrogens {
            old: AstValueAst::Lit(3),
            new: AstValueAst::Lit(2),
        },
        "ValueAst.Lit(3)",
        "ValueAst.Lit(2)",
        "AtomFieldChange.ImplicitHydrogens(old=ValueAst.Lit(3), new=ValueAst.Lit(2))"
    )]
    #[case::lone_pairs(
        AstAtomFieldChange::LonePairs {
            old: AstValueAst::Lit(1),
            new: AstValueAst::Lit(2),
        },
        "ValueAst.Lit(1)",
        "ValueAst.Lit(2)",
        "AtomFieldChange.LonePairs(old=ValueAst.Lit(1), new=ValueAst.Lit(2))"
    )]
    #[case::spin(
        AstAtomFieldChange::Spin {
            old: AstSpinStateAst {
                unpaired: AstValueAst::Lit(0),
                multiplicity: AstValueAst::Lit(1),
            },
            new: AstSpinStateAst {
                unpaired: AstValueAst::Lit(1),
                multiplicity: AstValueAst::Lit(2),
            },
        },
        "SpinStateAst(ValueAst.Lit(0), ValueAst.Lit(1))",
        "SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2))",
        "AtomFieldChange.Spin(old=SpinStateAst(ValueAst.Lit(0), ValueAst.Lit(1)), new=SpinStateAst(ValueAst.Lit(1), ValueAst.Lit(2)))"
    )]
    fn test_atom_field_change_repr(
        #[case] change: AstAtomFieldChange,
        #[case] old: &str,
        #[case] new: &str,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let change =
                into_py_variant(py, AtomFieldChange::from_ast(py, &change).unwrap()).unwrap();
            let bound = change.bind(py).as_any();
            assert_eq!(
                bound
                    .getattr("old")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                old
            );
            assert_eq!(
                bound
                    .getattr("new")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                new
            );
            assert_eq!(bound.repr().unwrap().extract::<String>().unwrap(), expected);
        });
    }

    #[rstest]
    #[case::element(AstAtomFieldChange::Element {
        old: AstElementAst::Lit(Element::C),
        new: AstElementAst::Lit(Element::N),
    })]
    #[case::isotope_mass(AstAtomFieldChange::IsotopeMass {
        old: AstIsotopeMassAst::Lit(12),
        new: AstIsotopeMassAst::Lit(13),
    })]
    #[case::charge(AstAtomFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::implicit_hydrogens(AstAtomFieldChange::ImplicitHydrogens {
        old: AstValueAst::Lit(3),
        new: AstValueAst::Lit(2),
    })]
    #[case::lone_pairs(AstAtomFieldChange::LonePairs {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::spin(AstAtomFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_atom_field_change_inverse(#[case] change: AstAtomFieldChange) {
        Python::attach(|py| {
            let mirror = AtomFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::order(AstBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::charge(AstBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_bond_field_change_roundtrip(#[case] change: AstBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                BondFieldChange::from_ast(py, &change).unwrap().to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::order(AstBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    #[case::charge(AstBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_bond_field_change_inverse(#[case] change: AstBondFieldChange) {
        Python::attach(|py| {
            let mirror = BondFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::order(AstDativeBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    fn test_dative_bond_field_change_roundtrip(#[case] change: AstDativeBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                DativeBondFieldChange::from_ast(py, &change)
                    .unwrap()
                    .to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::order(AstDativeBondFieldChange::Order {
        old: AstValueAst::Lit(1),
        new: AstValueAst::Lit(2),
    })]
    fn test_dative_bond_field_change_inverse(#[case] change: AstDativeBondFieldChange) {
        Python::attach(|py| {
            let mirror = DativeBondFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::electrons(AstAromaticSystemFieldChange::Electrons {
        old: AstElectronCountsAst::Undetermined,
        new: AstElectronCountsAst::Lit(vec![1, 1, 1]),
    })]
    #[case::charge(AstAromaticSystemFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::spin(AstAromaticSystemFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_aromatic_system_field_change_roundtrip(#[case] change: AstAromaticSystemFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                AromaticSystemFieldChange::from_ast(py, &change)
                    .unwrap()
                    .to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::electrons(AstAromaticSystemFieldChange::Electrons {
        old: AstElectronCountsAst::Undetermined,
        new: AstElectronCountsAst::Lit(vec![1, 1, 1]),
    })]
    #[case::charge(AstAromaticSystemFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(-1),
    })]
    #[case::spin(AstAromaticSystemFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(1),
            multiplicity: AstValueAst::Lit(2),
        },
    })]
    fn test_aromatic_system_field_change_inverse(#[case] change: AstAromaticSystemFieldChange) {
        Python::attach(|py| {
            let mirror = AromaticSystemFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::electrons(AstMulticenterBondFieldChange::Electrons {
        old: AstElectronCountsAst::Lit(vec![1, 0, 1]),
        new: AstElectronCountsAst::Lit(vec![2, 0, 1]),
    })]
    #[case::charge(AstMulticenterBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstMulticenterBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(2),
            multiplicity: AstValueAst::Lit(3),
        },
    })]
    fn test_multicenter_bond_field_change_roundtrip(#[case] change: AstMulticenterBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterBondFieldChange::from_ast(py, &change)
                    .unwrap()
                    .to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::electrons(AstMulticenterBondFieldChange::Electrons {
        old: AstElectronCountsAst::Lit(vec![1, 0, 1]),
        new: AstElectronCountsAst::Lit(vec![2, 0, 1]),
    })]
    #[case::charge(AstMulticenterBondFieldChange::Charge {
        old: AstValueAst::Lit(0),
        new: AstValueAst::Lit(1),
    })]
    #[case::spin(AstMulticenterBondFieldChange::Spin {
        old: AstSpinStateAst {
            unpaired: AstValueAst::Lit(0),
            multiplicity: AstValueAst::Lit(1),
        },
        new: AstSpinStateAst {
            unpaired: AstValueAst::Lit(2),
            multiplicity: AstValueAst::Lit(3),
        },
    })]
    fn test_multicenter_bond_field_change_inverse(#[case] change: AstMulticenterBondFieldChange) {
        Python::attach(|py| {
            let mirror = MulticenterBondFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::kind(AstNoncovalentBondFieldChange::Kind {
        old: AstNoncovalentBondKindAst::Undetermined,
        new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
    })]
    fn test_noncovalent_bond_field_change_roundtrip(#[case] change: AstNoncovalentBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                NoncovalentBondFieldChange::from_ast(py, &change)
                    .unwrap()
                    .to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::kind(AstNoncovalentBondFieldChange::Kind {
        old: AstNoncovalentBondKindAst::Undetermined,
        new: AstNoncovalentBondKindAst::Lit(AstNoncovalentBondKind::HydrogenBond),
    })]
    fn test_noncovalent_bond_field_change_inverse(#[case] change: AstNoncovalentBondFieldChange) {
        Python::attach(|py| {
            let mirror = NoncovalentBondFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::geometry_unknown(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_atom_field_change_roundtrip(#[case] change: AstStereoAtomFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                StereoAtomFieldChange::from_ast(py, &change)
                    .unwrap()
                    .to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        true,
    )]
    #[case::different(
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(1),
            ),
        },
        false,
    )]
    fn test_stereo_atom_field_change_eq(
        #[case] lhs: AstStereoAtomFieldChange,
        #[case] rhs: AstStereoAtomFieldChange,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = StereoAtomFieldChange::from_ast(py, &lhs).unwrap();
            let rhs = StereoAtomFieldChange::from_ast(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
        },
        "StereoConfigurationAst.Undetermined()",
        "StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined())",
        "StereoAtomFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined()))",
    )]
    #[case::coset_resolved(
        AstStereoAtomFieldChange::Configuration {
            old: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Undetermined,
            ),
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::Tetrahedral,
                AstStereoCosetAst::Lit(1),
            ),
        },
        "StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined())",
        "StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Lit(1))",
        "StereoAtomFieldChange.Configuration(old=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Undetermined()), new=StereoConfigurationAst.Kinded(StereoKind.Tetrahedral, StereoCosetAst.Lit(1)))",
    )]
    fn test_stereo_atom_field_change_repr(
        #[case] change: AstStereoAtomFieldChange,
        #[case] old: &str,
        #[case] new: &str,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let change =
                into_py_variant(py, StereoAtomFieldChange::from_ast(py, &change).unwrap()).unwrap();
            let bound = change.bind(py).as_any();
            assert_eq!(
                bound
                    .getattr("old")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                old
            );
            assert_eq!(
                bound
                    .getattr("new")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                new
            );
            assert_eq!(bound.repr().unwrap().extract::<String>().unwrap(), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoAtomFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::Tetrahedral,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_atom_field_change_inverse(#[case] change: AstStereoAtomFieldChange) {
        Python::attach(|py| {
            let mirror = StereoAtomFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }

    #[rstest]
    #[case::geometry_unknown(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_bond_field_change_roundtrip(#[case] change: AstStereoBondFieldChange) {
        Python::attach(|py| {
            assert_eq!(
                StereoBondFieldChange::from_ast(py, &change)
                    .unwrap()
                    .to_ast(py),
                change
            );
        });
    }

    #[rstest]
    #[case::equal(
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        true,
    )]
    #[case::different(
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(1),
            ),
        },
        false,
    )]
    fn test_stereo_bond_field_change_eq(
        #[case] lhs: AstStereoBondFieldChange,
        #[case] rhs: AstStereoBondFieldChange,
        #[case] expected: bool,
    ) {
        Python::attach(|py| {
            let lhs = StereoBondFieldChange::from_ast(py, &lhs).unwrap();
            let rhs = StereoBondFieldChange::from_ast(py, &rhs).unwrap();
            assert_eq!(lhs.__eq__(&rhs, py), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Undetermined,
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
        },
        "StereoConfigurationAst.Undetermined()",
        "StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined())",
        "StereoBondFieldChange.Configuration(old=StereoConfigurationAst.Undetermined(), new=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined()))",
    )]
    #[case::coset_resolved(
        AstStereoBondFieldChange::Configuration {
            old: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Undetermined,
            ),
            new: AstStereoConfigurationAst::Kinded(
                AstStereoKind::CisTrans,
                AstStereoCosetAst::Lit(1),
            ),
        },
        "StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined())",
        "StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Lit(1))",
        "StereoBondFieldChange.Configuration(old=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Undetermined()), new=StereoConfigurationAst.Kinded(StereoKind.CisTrans, StereoCosetAst.Lit(1)))",
    )]
    fn test_stereo_bond_field_change_repr(
        #[case] change: AstStereoBondFieldChange,
        #[case] old: &str,
        #[case] new: &str,
        #[case] expected: &str,
    ) {
        Python::attach(|py| {
            let change =
                into_py_variant(py, StereoBondFieldChange::from_ast(py, &change).unwrap()).unwrap();
            let bound = change.bind(py).as_any();
            assert_eq!(
                bound
                    .getattr("old")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                old
            );
            assert_eq!(
                bound
                    .getattr("new")
                    .unwrap()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                new
            );
            assert_eq!(bound.repr().unwrap().extract::<String>().unwrap(), expected);
        });
    }

    #[rstest]
    #[case::geometry_unknown(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Undetermined,
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
    })]
    #[case::coset_resolved(AstStereoBondFieldChange::Configuration {
        old: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Undetermined,
        ),
        new: AstStereoConfigurationAst::Kinded(
            AstStereoKind::CisTrans,
            AstStereoCosetAst::Lit(1),
        ),
    })]
    fn test_stereo_bond_field_change_inverse(#[case] change: AstStereoBondFieldChange) {
        Python::attach(|py| {
            let mirror = StereoBondFieldChange::from_ast(py, &change).unwrap();
            let inverse = mirror.inverse(py).unwrap();
            assert_eq!(
                inverse.bind(py).borrow().to_ast(py),
                change.clone().inverse()
            );
            let roundtrip = inverse.bind(py).borrow().inverse(py).unwrap();
            assert_eq!(roundtrip.bind(py).borrow().to_ast(py), change);
        });
    }
}

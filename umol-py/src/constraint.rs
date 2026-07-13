//! Shared constraint value/scope leaves mirroring `umol_ast::ast::constraint`: the
//! aromatic/multicenter valence states, ring scope, and ring membership.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_ast::ast::{
    AromaticValenceAst as AstAromaticValenceAst, MulticenterValenceAst as AstMulticenterValenceAst,
    RingMembershipAst as AstRingMembershipAst, RingScope as AstRingScope,
    TetrahedralStereoAst as AstTetrahedralStereoAst,
};

use crate::convert::{hash_ast, into_py_variant, variant_repr};
use crate::stereo::{TetrahedralStereo, TetrahedralStereoAst};
use crate::value::{ValueArg, ValueAst};

/// Aromatic-valence state: undetermined, explicitly not aromatic, or aromatic with
/// an aromatic-valence count. `Aromatic` coerces `int | ValueAst` on construction.
#[pyclass]
pub enum AromaticValenceAst {
    Undetermined(),
    NotAromatic(),
    Aromatic(ValueArg),
}

#[pymethods]
impl AromaticValenceAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            AromaticValenceAst::Undetermined() => ("Undetermined", 0),
            AromaticValenceAst::NotAromatic() => ("NotAromatic", 0),
            AromaticValenceAst::Aromatic(_) => ("Aromatic", 1),
        };
        variant_repr(slf.bind(py).as_any(), "AromaticValenceAst", variant, arity)
    }
}

impl AromaticValenceAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstAromaticValenceAst) -> PyResult<Self> {
        Ok(match ast {
            AstAromaticValenceAst::Undetermined => Self::Undetermined(),
            AstAromaticValenceAst::NotAromatic => Self::NotAromatic(),
            AstAromaticValenceAst::Aromatic(v) => Self::Aromatic(ValueArg::Ast(into_py_variant(
                py,
                ValueAst::from_ast(py, v)?,
            )?)),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstAromaticValenceAst {
        match self {
            Self::Undetermined() => AstAromaticValenceAst::Undetermined,
            Self::NotAromatic() => AstAromaticValenceAst::NotAromatic,
            Self::Aromatic(v) => AstAromaticValenceAst::Aromatic(v.to_ast(py)),
        }
    }
}

/// Setter coercion for `aromatic_valence`: `False` → not aromatic, an `int`/`ValueAst`
/// → aromatic with that valence, or an `AromaticValenceAst` passthrough.
#[derive(FromPyObject)]
pub(crate) enum AromaticValenceArg {
    Flag(bool),
    Value(ValueArg),
    Ast(Py<AromaticValenceAst>),
}

impl AromaticValenceArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> PyResult<AstAromaticValenceAst> {
        Ok(match self {
            AromaticValenceArg::Flag(false) => AstAromaticValenceAst::NotAromatic,
            AromaticValenceArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "aromatic_valence = True is not meaningful; use an int count or False",
                ))
            }
            AromaticValenceArg::Value(v) => AstAromaticValenceAst::Aromatic(v.to_ast(py)),
            AromaticValenceArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// Multicenter-valence state: undetermined, explicitly not multicenter, or
/// multicenter with a multicenter-valence count. `Multicenter` coerces
/// `int | ValueAst` on construction.
#[pyclass]
pub enum MulticenterValenceAst {
    Undetermined(),
    NotMulticenter(),
    Multicenter(ValueArg),
}

#[pymethods]
impl MulticenterValenceAst {
    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            MulticenterValenceAst::Undetermined() => ("Undetermined", 0),
            MulticenterValenceAst::NotMulticenter() => ("NotMulticenter", 0),
            MulticenterValenceAst::Multicenter(_) => ("Multicenter", 1),
        };
        variant_repr(
            slf.bind(py).as_any(),
            "MulticenterValenceAst",
            variant,
            arity,
        )
    }
}

impl MulticenterValenceAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstMulticenterValenceAst) -> PyResult<Self> {
        Ok(match ast {
            AstMulticenterValenceAst::Undetermined => Self::Undetermined(),
            AstMulticenterValenceAst::NotMulticenter => Self::NotMulticenter(),
            AstMulticenterValenceAst::Multicenter(v) => Self::Multicenter(ValueArg::Ast(
                into_py_variant(py, ValueAst::from_ast(py, v)?)?,
            )),
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstMulticenterValenceAst {
        match self {
            Self::Undetermined() => AstMulticenterValenceAst::Undetermined,
            Self::NotMulticenter() => AstMulticenterValenceAst::NotMulticenter,
            Self::Multicenter(v) => AstMulticenterValenceAst::Multicenter(v.to_ast(py)),
        }
    }
}

/// Setter coercion for `multicenter_valence`: `False` → not multicenter, an
/// `int`/`ValueAst` → multicenter with that valence, or a `MulticenterValenceAst`
/// passthrough.
#[derive(FromPyObject)]
pub(crate) enum MulticenterValenceArg {
    Flag(bool),
    Value(ValueArg),
    Ast(Py<MulticenterValenceAst>),
}

impl MulticenterValenceArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> PyResult<AstMulticenterValenceAst> {
        Ok(match self {
            MulticenterValenceArg::Flag(false) => AstMulticenterValenceAst::NotMulticenter,
            MulticenterValenceArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "multicenter_valence = True is not meaningful; use an int count or False",
                ))
            }
            MulticenterValenceArg::Value(v) => AstMulticenterValenceAst::Multicenter(v.to_ast(py)),
            MulticenterValenceArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// Setter coercion for `tetrahedral_stereo`: `False` → not stereogenic, a
/// `TetrahedralStereo` (`Ccw`/`Cw`) → that coset, or a `TetrahedralStereoAst`
/// passthrough.
#[derive(FromPyObject)]
pub(crate) enum TetrahedralStereoArg {
    Flag(bool),
    Config(TetrahedralStereo),
    Ast(Py<TetrahedralStereoAst>),
}

impl TetrahedralStereoArg {
    pub(crate) fn to_ast(&self, py: Python<'_>) -> PyResult<AstTetrahedralStereoAst> {
        Ok(match self {
            TetrahedralStereoArg::Flag(false) => AstTetrahedralStereoAst::NotStereo,
            TetrahedralStereoArg::Flag(true) => {
                return Err(PyValueError::new_err(
                    "tetrahedral_stereo = True is not meaningful; use TetrahedralStereo.Ccw/Cw or False",
                ))
            }
            TetrahedralStereoArg::Config(ts) => ts.to_ast(),
            TetrahedralStereoArg::Ast(a) => a.bind(py).borrow().to_ast(py),
        })
    }
}

/// Ring scope: all rings, or rings of a given size.
#[pyclass]
pub enum RingScope {
    All(),
    Size(u8),
}

#[pymethods]
impl RingScope {
    fn __eq__(&self, other: &Self) -> bool {
        self.to_ast() == other.to_ast()
    }

    fn __hash__(&self) -> u64 {
        hash_ast(&self.to_ast())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            RingScope::All() => ("All", 0),
            RingScope::Size(_) => ("Size", 1),
        };
        variant_repr(slf.bind(py).as_any(), "RingScope", variant, arity)
    }
}

impl RingScope {
    pub(crate) fn from_ast(ast: &AstRingScope) -> Self {
        match ast {
            AstRingScope::All => Self::All(),
            AstRingScope::Size(size) => Self::Size(*size),
        }
    }

    pub(crate) fn to_ast(&self) -> AstRingScope {
        match self {
            Self::All() => AstRingScope::All,
            Self::Size(size) => AstRingScope::Size(*size),
        }
    }
}

/// Ring-membership fact: a ring scope and a membership count.
#[pyclass]
pub struct RingMembershipAst {
    #[pyo3(get)]
    scope: Py<RingScope>,
    #[pyo3(get)]
    count: Py<ValueAst>,
}

#[pymethods]
impl RingMembershipAst {
    #[new]
    fn new(py: Python<'_>, scope: Py<RingScope>, count: ValueArg) -> PyResult<Self> {
        Ok(RingMembershipAst {
            scope,
            count: count.to_py(py)?,
        })
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_ast(py) == other.to_ast(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_ast(&self.to_ast(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "RingMembershipAst({}, {})",
            self.scope.bind(py).as_any().repr()?.extract::<String>()?,
            self.count.bind(py).as_any().repr()?.extract::<String>()?,
        ))
    }
}

impl RingMembershipAst {
    pub(crate) fn from_ast(py: Python<'_>, ast: &AstRingMembershipAst) -> PyResult<Self> {
        Ok(RingMembershipAst {
            scope: into_py_variant(py, RingScope::from_ast(&ast.scope))?,
            count: into_py_variant(py, ValueAst::from_ast(py, &ast.count)?)?,
        })
    }

    pub(crate) fn to_ast(&self, py: Python<'_>) -> AstRingMembershipAst {
        AstRingMembershipAst::new(
            self.scope.bind(py).borrow().to_ast(),
            self.count.bind(py).borrow().to_ast(py),
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    #[rstest]
    #[case(AstAromaticValenceAst::Undetermined)]
    #[case(AstAromaticValenceAst::NotAromatic)]
    #[case(AstAromaticValenceAst::aromatic(1))]
    fn test_aromatic_valence_ast_roundtrip(#[case] ast: AstAromaticValenceAst) {
        Python::attach(|py| {
            assert_eq!(
                AromaticValenceAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstMulticenterValenceAst::Undetermined)]
    #[case(AstMulticenterValenceAst::NotMulticenter)]
    #[case(AstMulticenterValenceAst::multicenter(2))]
    fn test_multicenter_valence_ast_roundtrip(#[case] ast: AstMulticenterValenceAst) {
        Python::attach(|py| {
            assert_eq!(
                MulticenterValenceAst::from_ast(py, &ast)
                    .unwrap()
                    .to_ast(py),
                ast
            );
        });
    }

    #[rstest]
    #[case(AstRingScope::All)]
    #[case(AstRingScope::Size(6))]
    fn test_ring_scope_roundtrip(#[case] ast: AstRingScope) {
        assert_eq!(RingScope::from_ast(&ast).to_ast(), ast);
    }

    #[rstest]
    #[case(AstRingMembershipAst::new(AstRingScope::All, 2))]
    #[case(AstRingMembershipAst::new(AstRingScope::Size(6), 1))]
    fn test_ring_membership_ast_roundtrip(#[case] ast: AstRingMembershipAst) {
        Python::attach(|py| {
            assert_eq!(
                RingMembershipAst::from_ast(py, &ast).unwrap().to_ast(py),
                ast
            );
        });
    }
}

//! Python values for reaction fingerprint results.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::prelude::*;
use umol_graph::fingerprint::Side as GraphReactionSide;

/// Side of a reaction from which a fingerprint feature originates.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReactionSide {
    Reactant,
    Product,
}

impl ReactionSide {
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Rust-to-Python conversion is used by reaction fingerprint operations"
        )
    )]
    pub(crate) fn from_rust(side: GraphReactionSide) -> Self {
        match side {
            GraphReactionSide::Reactant => Self::Reactant,
            GraphReactionSide::Product => Self::Product,
        }
    }

    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "Python-to-Rust conversion is used by reaction fingerprint operations"
        )
    )]
    pub(crate) fn to_rust(self) -> GraphReactionSide {
        match self {
            Self::Reactant => GraphReactionSide::Reactant,
            Self::Product => GraphReactionSide::Product,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::reactant(GraphReactionSide::Reactant, ReactionSide::Reactant)]
    #[case::product(GraphReactionSide::Product, ReactionSide::Product)]
    fn test_reaction_side_from_rust(
        #[case] side: GraphReactionSide,
        #[case] expected: ReactionSide,
    ) {
        assert_eq!(ReactionSide::from_rust(side), expected);
    }

    #[rstest]
    #[case::reactant(ReactionSide::Reactant, GraphReactionSide::Reactant)]
    #[case::product(ReactionSide::Product, GraphReactionSide::Product)]
    fn test_reaction_side_to_rust(#[case] side: ReactionSide, #[case] expected: GraphReactionSide) {
        assert_eq!(side.to_rust(), expected);
    }

    #[rstest]
    #[case::reactant(
        ReactionSide::Reactant,
        ReactionSide::Reactant,
        ReactionSide::Product,
        "ReactionSide.Reactant"
    )]
    #[case::product(
        ReactionSide::Product,
        ReactionSide::Product,
        ReactionSide::Reactant,
        "ReactionSide.Product"
    )]
    fn test_reaction_side_python_value(
        #[case] side: ReactionSide,
        #[case] equal: ReactionSide,
        #[case] unequal: ReactionSide,
        #[case] expected_repr: &str,
    ) {
        Python::attach(|py| {
            let side = Py::new(py, side).unwrap();
            let equal = Py::new(py, equal).unwrap();
            let unequal = Py::new(py, unequal).unwrap();

            assert!(side.bind(py).as_any().eq(equal.bind(py).as_any()).unwrap());
            assert!(!side
                .bind(py)
                .as_any()
                .eq(unequal.bind(py).as_any())
                .unwrap());
            assert_eq!(
                side.bind(py).as_any().hash().unwrap(),
                equal.bind(py).as_any().hash().unwrap()
            );
            assert_eq!(
                side.bind(py)
                    .as_any()
                    .repr()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                expected_repr
            );
        });
    }
}

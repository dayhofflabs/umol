//! Model algebra.
//!
//! Algebraic operations on models:
//! - Ensembles of chemical models (e.g., conformers, resonance structures,
//!   positive fractional coefficients)
//! - Aggregates of chemical models (e.g., mixtures, fixed-stoichiometry
//!   complexes, positive integer coefficients)
//! - Processes (e.g., reactions, positive and negative integer coefficients)

use crate::Model;

/// Represents an ensemble of models with weights (e.g., conformers, resonance
/// structures, positive fractional coefficients)
pub trait Ensemble<M: Model> {
    /// Returns a slice of model-weight pairs
    fn components(&self) -> &[(M, Option<f64>)];

    /// Returns a vector of weights
    fn weights(&self) -> Option<Vec<f64>> {
        self.components().iter().map(|(_, w)| *w).collect()
    }

    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}

/// Represents an aggregate of models with coefficients (e.g., mixtures,
/// fixed-stoichiometry complexes, positive integer coefficients)
pub trait Aggregate<M: Model> {
    /// Returns a slice of model-coefficient pairs
    fn components(&self) -> &[(M, u32)];

    /// Returns a vector of coefficients
    fn coefficients(&self) -> Vec<u32> {
        self.components().iter().map(|(_, c)| *c).collect()
    }

    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}

/// Represents a process of chemical models (e.g., reactions, positive and
/// negative integer coefficients)
pub trait Process<M: Model> {
    /// Returns a slice of model-coefficient pairs
    fn components(&self) -> &[(M, i32)];

    /// Returns a vector of coefficients
    fn coefficients(&self) -> Vec<i32> {
        self.components().iter().map(|(_, c)| *c).collect()
    }

    /// Returns a vector of references to models
    fn models(&self) -> Vec<&M> {
        self.components().iter().map(|(m, _)| m).collect()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use serde::{Deserialize, Serialize};
    use serde_json;

    use super::*;
    use crate::Capability;

    /// A simple model that counts C, H, and O atoms
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct ElementCount {
        data: ElementCountData,
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct ElementCountData {
        c: usize,
        h: usize,
        o: usize,
    }

    impl Model for ElementCount {
        type Data = ElementCountData;

        fn data(&self) -> &Self::Data {
            &self.data
        }

        fn capabilities(&self) -> HashSet<Capability> {
            let mut caps = HashSet::new();
            caps.insert(Capability::local("element_count", 1));
            caps
        }
    }

    impl ElementCount {
        pub fn new(c: usize, h: usize, o: usize) -> Self {
            Self {
                data: ElementCountData { c, h, o },
            }
        }
    }

    // Test implementations of the algebra traits
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct SimpleEnsemble<M> {
        components: Vec<(M, Option<f64>)>,
    }

    impl<M> Ensemble<M> for SimpleEnsemble<M>
    where
        M: Model + Serialize + for<'de> Deserialize<'de>,
    {
        fn components(&self) -> &[(M, Option<f64>)] {
            &self.components
        }
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct SimpleAggregate<M> {
        components: Vec<(M, u32)>,
    }

    impl<M> Aggregate<M> for SimpleAggregate<M>
    where
        M: Model + Serialize + for<'de> Deserialize<'de>,
    {
        fn components(&self) -> &[(M, u32)] {
            &self.components
        }
    }

    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct SimpleProcess<M> {
        components: Vec<(M, i32)>,
    }

    impl<M> Process<M> for SimpleProcess<M>
    where
        M: Model + Serialize + for<'de> Deserialize<'de>,
    {
        fn components(&self) -> &[(M, i32)] {
            &self.components
        }
    }

    #[test]
    fn test_element_count_serialization() {
        let model = ElementCount::new(1, 4, 0); // CH4

        let json = serde_json::to_string(&model).unwrap();
        assert_eq!(json, r#"{"data":{"c":1,"h":4,"o":0}}"#);

        let deserialized: ElementCount = serde_json::from_str(&json).unwrap();
        assert_eq!(deserialized.data.c, 1);
        assert_eq!(deserialized.data.h, 4);
        assert_eq!(deserialized.data.o, 0);
    }

    #[test]
    fn test_ensemble_serialization() {
        let ch4 = ElementCount::new(1, 4, 0);
        let h2o = ElementCount::new(0, 2, 1);

        let ensemble = SimpleEnsemble {
            components: vec![(ch4, Some(0.5)), (h2o, Some(0.5))],
        };

        let json = serde_json::to_string(&ensemble).unwrap();
        let deserialized: SimpleEnsemble<ElementCount> = serde_json::from_str(&json).unwrap();

        assert_eq!(deserialized.components.len(), 2);
        assert_eq!(deserialized.components[0].1, Some(0.5));
        assert_eq!(deserialized.components[1].1, Some(0.5));
    }

    #[test]
    fn test_aggregate_serialization() {
        let ch4 = ElementCount::new(1, 4, 0);
        let h2o = ElementCount::new(0, 2, 1);

        let aggregate = SimpleAggregate {
            components: vec![(ch4, 2), (h2o, 3)],
        };

        let json = serde_json::to_string(&aggregate).unwrap();
        let deserialized: SimpleAggregate<ElementCount> = serde_json::from_str(&json).unwrap();

        assert_eq!(deserialized.components.len(), 2);
        assert_eq!(deserialized.components[0].1, 2);
        assert_eq!(deserialized.components[1].1, 3);
    }

    #[test]
    fn test_process_serialization() {
        let ch4 = ElementCount::new(1, 4, 0);
        let o2 = ElementCount::new(0, 0, 2);
        let co2 = ElementCount::new(1, 0, 2);
        let h2o = ElementCount::new(0, 2, 1);

        // CH4 + 2O2 -> CO2 + 2H2O
        let process = SimpleProcess {
            components: vec![
                (ch4, -1), // reactant
                (o2, -2),  // reactant
                (co2, 1),  // product
                (h2o, 2),  // product
            ],
        };

        let json = serde_json::to_string(&process).unwrap();
        let deserialized: SimpleProcess<ElementCount> = serde_json::from_str(&json).unwrap();

        assert_eq!(deserialized.components.len(), 4);
        assert_eq!(deserialized.components[0].1, -1);
        assert_eq!(deserialized.components[1].1, -2);
        assert_eq!(deserialized.components[2].1, 1);
        assert_eq!(deserialized.components[3].1, 2);
    }
}

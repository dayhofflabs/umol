// Traits and structs for pattern matching

use crate::graph::{AtomIndex, BondIndex, GraphMolecule};

pub trait PatternLanguage: Sized {
    type Error;
    type Match: PatternMatch;

    fn parse(pattern_string: &str) -> Result<Self, Self::Error>;
}

pub trait PatternMatch {
    fn matched_atoms(&self) -> &[AtomIndex];
    fn matched_bonds(&self) -> &[BondIndex];
}

pub struct Pattern<L: PatternLanguage> {
    language: L,
}

impl<L: PatternLanguage> Pattern<L> {
    pub fn new(pattern_string: &str) -> Result<Self, L::Error> {
        Ok(Self {
            language: L::parse(pattern_string)?,
        })
    }

    pub fn find_match<'a>(&self, molecule: &'a GraphMolecule) -> Option<L::Match> {
        // Implementation
        todo!()
    }

    pub fn find_matches(&self, molecule: &GraphMolecule) -> Vec<L::Match> {
        // Implementation
        todo!()
    }
}

pub struct SmartsLanguage;
pub struct SmartsMatch;

impl PatternMatch for SmartsMatch {
    fn matched_atoms(&self) -> &[AtomIndex] {
        // Implementation
        todo!()
    }

    fn matched_bonds(&self) -> &[BondIndex] {
        // Implementation
        todo!()
    }
}

pub struct SmartsError;

pub struct ReactionSmartsLanguage;
pub struct ReactionSmartsMatch;

impl PatternMatch for ReactionSmartsMatch {
    fn matched_atoms(&self) -> &[AtomIndex] {
        // Implementation
        todo!()
    }

    fn matched_bonds(&self) -> &[BondIndex] {
        // Implementation
        todo!()
    }
}

pub struct ReactionSmartsError;

impl PatternLanguage for SmartsLanguage {
    type Error = SmartsError;
    type Match = SmartsMatch;

    fn parse(pattern_string: &str) -> Result<Self, Self::Error> {
        // Implementation
        todo!()
    }
}


impl PatternLanguage for ReactionSmartsLanguage {
    type Error = ReactionSmartsError;
    type Match = ReactionSmartsMatch;

    fn parse(pattern_string: &str) -> Result<Self, Self::Error> {
        // Implementation
        todo!()
    }
}

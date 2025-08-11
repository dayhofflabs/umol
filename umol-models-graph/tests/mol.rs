//! MOL format compliance suite

mod common;

#[path = "mol/real_world.rs"]
pub mod real_world;

#[path = "mol/properties.rs"] 
pub mod properties;

#[path = "mol/edge_cases.rs"]
pub mod edge_cases;

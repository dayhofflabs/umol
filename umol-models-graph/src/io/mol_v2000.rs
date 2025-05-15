//! MOL v2000 file parser and writer.
//!
//! # Components of the MOL v2000 Format (based on the Connection Table Structure):
//!
//! * [Counts line](read_mol_v2000.rs)
//! * [Atom block](atom.rs)
//! * [Bond block](bond.rs)
//! // * Atom list block
//! // * Stext block
//! * [Properties block](property.rs) incl. [3D features block](conformer.rs)
//!
//! # References
//! * MOL v2000: https://en.wikipedia.org/wiki/Chemical_table_file
//! * RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules
//! * ChemAxon: https://docs.chemaxon.com/display/docs/formats_mdl-molfiles-rgfiles-sdfiles-rxnfiles-rdfiles-formats.md
//! * CDK: https://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/io/MDLV2000Reader.html

mod atom;
mod bond;
mod conformer;
mod property;
mod reader;
// mod writer;

#[cfg(test)]
mod tests;

pub use reader::read_mol_v2000;
// pub use writer::write_mol_v2000;
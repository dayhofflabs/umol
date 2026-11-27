//! Pure-Rust trail-backed edge-colored individualization--refinement engine.

mod components;
mod partition;
mod preprocess;
mod refine;
mod search;
#[cfg(test)]
mod tests;

pub(super) use search::TrailIrWorkspace;

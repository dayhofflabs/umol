//! Graph connectivity.
//!
//! Current operations provide breadth-first connected components and Tarjan
//! biconnected components. The biconnected result contains vertex sets for
//! blocks of at least three vertices: bridge blocks are omitted and
//! articulation points are not returned separately. See
//! [Tarjan, *Depth-First Search and Linear Graph Algorithms*
//! (1972)](https://doi.org/10.1137/0201010).

mod biconnected;
mod components;

pub use biconnected::BiconnectedComponentsAlgorithm;
pub use components::ConnectedComponentsAlgorithm;

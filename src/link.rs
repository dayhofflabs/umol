// Core AtomLink trait

use super::atom::AtomSite;

pub trait AtomLink<A: AtomSite + Sized> {
    type SiteRef;
}

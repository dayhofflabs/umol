// Core AtomLink trait

use super::atom::AtomSite;

pub trait AtomLink<A: AtomSite + Sized> {
    type SiteRef;
    fn between(&self) -> (Self::SiteRef, Self::SiteRef);
}

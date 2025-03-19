use crate::element::Element;

pub trait AtomSite {
    fn element(&self) -> Option<Element>;
}

pub trait AtomLink {
    type Site: AtomSite;
    type SiteRef;
}

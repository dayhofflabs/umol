use rstest::rstest;
use umol_graph_ir::frag;
use umol_graph_ir::ir::{AtomAst, AtomId, BondAst, BondId, ElementForm, Fragment, Port};

#[rstest]
fn test_frag_port() {
    // a `^name` marker declares a port on the bonded real atom; colour read off the bond op
    let fragment: Fragment = frag!((c: C) - ^x);

    assert_eq!(fragment.body().atoms().count(), 1);
    assert_eq!(fragment.body().bonds().count(), 0);
    assert_eq!(
        fragment.ports(),
        &[Port {
            atom: AtomId(0),
            bond: BondAst::from_order(1),
            name: Some("x".to_string()),
        }]
    );
}

#[rstest]
fn test_frag_multiple_ports() {
    // a real atom flanked by two ports of different colours
    let fragment: Fragment = frag!(^a - (c: C) = ^b);

    assert_eq!(fragment.body().atoms().count(), 1);
    assert_eq!(
        fragment.ports(),
        &[
            Port {
                atom: AtomId(0),
                bond: BondAst::from_order(1),
                name: Some("a".to_string()),
            },
            Port {
                atom: AtomId(0),
                bond: BondAst::from_order(2),
                name: Some("b".to_string()),
            },
        ]
    );
}

#[rstest]
fn test_frag_spec_port() {
    // a port colour carried by a rich DSL bond spec
    let fragment: Fragment = frag!((c: C) -[ "1#a" ]- ^r);

    assert_eq!(
        fragment.ports(),
        &[Port {
            atom: AtomId(0),
            bond: "1#a".parse::<BondAst>().unwrap(),
            name: Some("r".to_string()),
        }]
    );
}

#[rstest]
fn test_frag_finish_open() {
    // closing a fragment with a free port caps it with a wildcard atom via the port colour
    let pattern = frag!((c: C) - ^x).finish_open();

    assert_eq!(pattern.atoms().count(), 2);
    assert_eq!(
        pattern.atom(AtomId(1)).ast,
        &AtomAst::new(ElementForm::undetermined())
    );
    assert_eq!(pattern.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    assert_eq!(pattern.bond(BondId(0)).ast, &BondAst::from_order(1));
}

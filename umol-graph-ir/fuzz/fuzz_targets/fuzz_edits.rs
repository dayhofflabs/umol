#![no_main]

use libfuzzer_sys::fuzz_target;
use umol_edn::FromEdn;
use umol_graph_ir::dsl::{EditsDsl, MoleculeDefaults};
use umol_graph_ir::ir::{FromIr, IntoIr};

fuzz_target!(|data: &str| {
    if let Ok(dsl) = EditsDsl::from_edn_str(data) {
        let defaults = MoleculeDefaults::new();
        let edits = dsl.into_ir(&defaults);
        let rendered = EditsDsl::from_ir(&edits, &defaults).to_string();
        EditsDsl::from_edn_str(&rendered).expect("rendered edits must parse");
    }
});

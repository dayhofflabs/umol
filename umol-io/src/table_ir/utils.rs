//! Utility functions for TableIR.

use std::collections::BTreeMap;

use umol_chem::element::Element;

/// Format sum formula according to Hill notation
pub(super) fn format_sum_formula(
    c_count: usize,
    h_count: usize,
    atom_counts: BTreeMap<[u8; 2], (Element, usize)>,
    charge: i32,
) -> String {
    let mut sum_formula = String::new();

    // Carbon first
    if c_count > 1 {
        sum_formula.push_str(&format!("C{}", c_count));
    } else if c_count == 1 {
        sum_formula.push('C');
    }

    // Hydrogen second
    if h_count > 1 {
        sum_formula.push_str(&format!("H{}", h_count));
    } else if h_count == 1 {
        sum_formula.push('H');
    }

    // Other elements alphabetically by symbol
    for (_, (element, count)) in atom_counts {
        if count > 1 {
            sum_formula.push_str(&format!("{}{}", element, count));
        } else {
            sum_formula.push_str(&element.to_string());
        }
    }

    // Charge at the end
    if charge != 0 {
        if charge == 1 {
            sum_formula.push('+');
        } else if charge == -1 {
            sum_formula.push('-');
        } else {
            sum_formula.push_str(&format!("{:+}", charge));
        }
    }

    sum_formula
}

/// Convert element symbol to [u8; 2] key for alphabetical sorting
pub(super) fn element_symbol_key(element: Element) -> [u8; 2] {
    let symbol = element.symbol();
    let bytes = symbol.as_bytes();
    [
        bytes.first().copied().unwrap_or(0),
        bytes.get(1).copied().unwrap_or(0),
    ]
}

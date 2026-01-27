Full valence pattern table (snapshot)

Copy the contents below back into `ValencePatternTable::organic_strict()` when ready. This is the version with H, B, C, N, O, F, P, S, Cl, Br, I that was previously in `checker.rs`.

```rust
// Patterns for the organic elements and hydrogen. Uses umol_data::e!(X) macro.
Self { patterns: vec![
    // H
    ValencePattern { element: Some(e!(H)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // HR
    ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // H^
    ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // H+
    ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // H-

    // B
    ValencePattern { element: Some(e!(B)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BR3, B(OH)3, B2O3
    ValencePattern { element: Some(e!(B)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HBR2
    ValencePattern { element: Some(e!(B)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2BR
    ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3B
    ValencePattern { element: Some(e!(B)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // BR2^
    ValencePattern { element: Some(e!(B)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(1) }, // HBR^
    ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(1) }, // H2B^
    ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(3) }, // B^^^
    ValencePattern { element: Some(e!(B)), bond_sum: Some(4), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BR4-
    ValencePattern { element: Some(e!(B)), bond_sum: Some(3), charge: Some(-1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HBR3-
    ValencePattern { element: Some(e!(B)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2BR2-
    ValencePattern { element: Some(e!(B)), bond_sum: Some(1), charge: Some(-1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3BR-
    ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(4), unpaired: Some(0) }, // H4B-

    // C
    ValencePattern { element: Some(e!(C)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // CR4, CO2, HOC(=O)R
    ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HCR3, HC(=O)R, HC(=O)OH
    ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2CR2, H2C=O
    ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3CR
    ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(4), unpaired: Some(0) }, // H4C
    ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // CR3^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(1) }, // HCR2^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(1) }, // H2CR^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(3), unpaired: Some(1) }, // H3C^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(2) }, // CR2^^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(2) }, // HCR^^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(2) }, // H2C^^
    ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(2) }, // C^^, CO
    ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // CR3+
    ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HCR2+
    ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2CR+
    ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3C+
    ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // CR3-
    ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HCR2-
    ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(-1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2CR-
    ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3C-

    // N
    ValencePattern { element: Some(e!(N)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // NR3, HONO
    ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HNR2
    ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2NR
    ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3N
    ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // NR2^, NO
    ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(1) }, // HNR^
    ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(1) }, // H2N^
    ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(2) }, // NR^^
    ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(2) }, // HN^^
    ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(3) }, // N^^^
    ValencePattern { element: Some(e!(N)), bond_sum: Some(4), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // NR4+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(3), charge: Some(1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HNR3+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2NR2+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3NR+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(1), implicit_hydrogens: Some(4), unpaired: Some(0) }, // H4N+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // NO+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // NR2-
    ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(-1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HNR-
    ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2N-
    ValencePattern { element: Some(e!(N)), bond_sum: Some(5), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // NRO2, HONO2
    ValencePattern { element: Some(e!(N)), bond_sum: Some(4), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // NO2+
    ValencePattern { element: Some(e!(N)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // NO2

    // O
    ValencePattern { element: Some(e!(O)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // OR2
    ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HOR
    ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2O
    ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // OR^
    ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // HO^
    ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(2) }, // O^^
    ValencePattern { element: Some(e!(O)), bond_sum: Some(3), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // OR3+
    ValencePattern { element: Some(e!(O)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HOR2+
    ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2OR+
    ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3O+
    ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // OR-
    ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HO-

    // F
    ValencePattern { element: Some(e!(F)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // FR
    ValencePattern { element: Some(e!(F)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HF
    ValencePattern { element: Some(e!(F)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // F^
    ValencePattern { element: Some(e!(F)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // F-

    // P
    ValencePattern { element: Some(e!(P)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // PR3, P(OH)3, P4O6
    ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HPR2
    ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2PR
    ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3P
    ValencePattern { element: Some(e!(P)), bond_sum: Some(4), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // PR4+
    ValencePattern { element: Some(e!(P)), bond_sum: Some(3), charge: Some(1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HPR3+
    ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2PR2+
    ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3PR+
    ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(1), implicit_hydrogens: Some(4), unpaired: Some(0) }, // H4P+
    ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // PR2^
    ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(1) }, // HPR^
    ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(1) }, // H2P^
    ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // PR2-
    ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(-1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HPR-
    ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2P-
    ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // P^^^
    ValencePattern { element: Some(e!(P)), bond_sum: Some(5), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // PR3(=O), PO(OH)3, P4O10

    // S
    ValencePattern { element: Some(e!(S)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // SR2
    ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HSR
    ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2S
    ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // SR^
    ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(1) }, // SH^
    ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(2) }, // S^^
    ValencePattern { element: Some(e!(S)), bond_sum: Some(3), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // SR3+
    ValencePattern { element: Some(e!(S)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HSR2+
    ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(1), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2SR+
    ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(1), implicit_hydrogens: Some(3), unpaired: Some(0) }, // H3S+
    ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // SR-
    ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HS-
    ValencePattern { element: Some(e!(S)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // S(=O)R2, SO2, SO(OH)2
    ValencePattern { element: Some(e!(S)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HS(=O)R
    ValencePattern { element: Some(e!(S)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2S=O
    ValencePattern { element: Some(e!(S)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // S(=O)(=O)R2, SO3, SO(OH)2
    ValencePattern { element: Some(e!(S)), bond_sum: Some(5), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HS(=O)(=O)R
    ValencePattern { element: Some(e!(S)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(2), unpaired: Some(0) }, // H2S(=O)=O

    // Cl
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClR
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HCl, HOCl, Cl2O
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // Cl^
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Cl-
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // ClO
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClO+
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClO-
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Cl(=O)OH, RClO
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HClO
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // ClO2
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClO2+
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClO2-
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(5), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Cl(=O)(=O)OH, RClO2
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HCl(=O)(=O)
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // ClO3
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClO3+
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // ClO3-
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(7), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Cl(=O)(=O)(=O)OH, RClO3, ClO4-
    ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(1) }, // HClO3

    // Br
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrR
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HBr, HOBr
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // Br^
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Br-
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // BrO
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrO+
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrO-
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Br(=O)OH, RBrO
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HBrO
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // BrO2
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrO2+
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrO2-
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(5), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Br(=O)(=O)OH, RBrO2
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HBrO2
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // BrO3
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrO3+
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // BrO3-
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(7), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // Br(=O)(=O)(=O)OH, RBrO3, BrO4-
    ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HBrO3

    // I
    ValencePattern { element: Some(e!(I)), bond_sum: Some(1), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IR
    ValencePattern { element: Some(e!(I)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HI, HOI
    ValencePattern { element: Some(e!(I)), bond_sum: Some(0), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // I^
    ValencePattern { element: Some(e!(I)), bond_sum: Some(0), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // I-
    ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // IO
    ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IO+
    ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IO-
    ValencePattern { element: Some(e!(I)), bond_sum: Some(3), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // I(=O)OH, RIO
    ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HIO
    ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // IO2
    ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IO2+
    ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IO2-
    ValencePattern { element: Some(e!(I)), bond_sum: Some(5), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // I(=O)(=O)OH, RIO2
    ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HIO
    ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HIO2
    ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(1) }, // IO3
    ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IO3+
    ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(-1), implicit_hydrogens: Some(0), unpaired: Some(0) }, // IO3-
    ValencePattern { element: Some(e!(I)), bond_sum: Some(7), charge: Some(0), implicit_hydrogens: Some(0), unpaired: Some(0) }, // I(=O)(=O)(=O)OH, RIO3, IO4-
    ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(0), implicit_hydrogens: Some(1), unpaired: Some(0) }, // HIO3
]} 
```



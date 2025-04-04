# Encoding of atomic properties

* Element: 1 byte
* Valence: 1 byte (charge + electronic conf.)
* Bonding: 1 byte
* 0b1111_1111 indicates information unavailable (or outside coding domain)

```rust
use map_macro::hash_map;
use element::Element;

static ELEMENT_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // period: 3 bits, group: 5 bits using Gray coding
    // (based on 32-group layout instead of traditional 18-group)
    Element::H  => 0b000_00000u8, // 1
    Element::He => 0b000_10000u8, // 32
    Element::Li => 0b001_00000u8, // 1
    Element::Be => 0b001_00001u8, // 2
    Element::B  => 0b001_10111u8, // 27
    Element::C  => 0b001_10110u8, // 28
    Element::N  => 0b001_10010u8, // 29
    Element::O  => 0b001_10011u8, // 30
    Element::F  => 0b001_10001u8, // 31
    Element::Ne => 0b001_10000u8, // 32
    Element::Na => 0b010_00000u8, // 1
    Element::Mg => 0b010_00001u8, // 2
    Element::Al => 0b010_10111u8, // 27
    Element::Si => 0b010_10110u8, // 28
    Element::P  => 0b010_10010u8, // 29
    Element::S  => 0b010_10011u8, // 30
    Element::Cl => 0b010_10001u8, // 31
    Element::Ar => 0b010_10000u8, // 32
    Element::K  => 0b011_00000u8, // 1
    Element::Ca => 0b011_00001u8, // 2
    Element::Sc => 0b011_11000u8, // 17
    Element::Ti => 0b011_11001u8, // 18
    Element::V  => 0b011_11011u8, // 19
    Element::Cr => 0b011_11010u8, // 20
    Element::Mn => 0b011_11110u8, // 21
    Element::Fe => 0b011_11111u8, // 22
    Element::Co => 0b011_11101u8, // 23
    Element::Ni => 0b011_11100u8, // 24
    Element::Cu => 0b011_10100u8, // 25
    Element::Zn => 0b011_10101u8, // 26
    Element::Ga => 0b011_10111u8, // 27
    Element::Ge => 0b011_10110u8, // 28
    Element::As => 0b011_10010u8, // 29
    Element::Se => 0b011_10011u8, // 30
    Element::Br => 0b011_10001u8, // 31
    Element::Kr => 0b011_10000u8, // 32
    Element::Rb => 0b100_00000u8, // 1
    Element::Sr => 0b100_00001u8, // 2
    Element::Y  => 0b100_11000u8, // 17
    Element::Zr => 0b100_11001u8, // 18
    Element::Nb => 0b100_11011u8, // 19
    Element::Mo => 0b100_11010u8, // 20
    Element::Tc => 0b100_11110u8, // 21
    Element::Ru => 0b100_11111u8, // 22
    Element::Rh => 0b100_11101u8, // 23
    Element::Pd => 0b100_11100u8, // 24
    Element::Ag => 0b100_10100u8, // 25
    Element::Cd => 0b100_10101u8, // 26
    Element::In => 0b100_10111u8, // 27
    Element::Sn => 0b100_10110u8, // 28
    Element::Sb => 0b100_10010u8, // 29
    Element::Te => 0b100_10011u8, // 30
    Element::I  => 0b100_10001u8, // 31
    Element::Xe => 0b100_10000u8, // 32
    Element::Cs => 0b101_00000u8, // 1
    Element::Ba => 0b101_00001u8, // 2
    Element::La => 0b101_00011u8, // 3
    Element::Ce => 0b101_00010u8, // 4
    Element::Pr => 0b101_00110u8, // 5
    Element::Nd => 0b101_00111u8, // 6
    Element::Pm => 0b101_00101u8, // 7
    Element::Sm => 0b101_00100u8, // 8
    Element::Eu => 0b101_01100u8, // 9
    Element::Gd => 0b101_01101u8, // 10
    Element::Tb => 0b101_01111u8, // 11
    Element::Dy => 0b101_01110u8, // 12
    Element::Ho => 0b101_01010u8, // 13
    Element::Er => 0b101_01011u8, // 14
    Element::Tm => 0b101_01001u8, // 15
    Element::Yb => 0b101_01000u8, // 16
    Element::Lu => 0b101_11000u8, // 17
    Element::Hf => 0b101_11001u8, // 18
    Element::Ta => 0b101_11011u8, // 19
    Element::W  => 0b101_11010u8, // 20
    Element::Re => 0b101_11110u8, // 21
    Element::Os => 0b101_11111u8, // 22
    Element::Ir => 0b101_11101u8, // 23
    Element::Pt => 0b101_11100u8, // 24
    Element::Au => 0b101_10100u8, // 25
    Element::Hg => 0b101_10101u8, // 26
    Element::Tl => 0b101_10111u8, // 27
    Element::Pb => 0b101_10110u8, // 28
    Element::Bi => 0b101_10010u8, // 29
    Element::Po => 0b101_10011u8, // 30
    Element::At => 0b101_10001u8, // 31
    Element::Rn => 0b101_10000u8, // 32
    Element::Fr => 0b110_00000u8, // 1
    Element::Ra => 0b110_00001u8, // 2
    Element::Ac => 0b110_00011u8, // 3
    Element::Th => 0b110_00010u8, // 4
    Element::Pa => 0b110_00110u8, // 5
    Element::U  => 0b110_00111u8, // 6
    Element::Np => 0b110_00101u8, // 7
    Element::Pu => 0b110_00100u8, // 8
    Element::Am => 0b110_01100u8, // 9
    Element::Cm => 0b110_01101u8, // 10
    Element::Bk => 0b110_01111u8, // 11
    Element::Cf => 0b110_01110u8, // 12
    Element::Es => 0b110_01010u8, // 13
    Element::Fm => 0b110_01011u8, // 14
    Element::Md => 0b110_01001u8, // 15
    Element::No => 0b110_01000u8, // 16
    Element::Lr => 0b110_11000u8, // 17
    Element::Rf => 0b110_11001u8, // 18
    Element::Db => 0b110_11011u8, // 19
    Element::Sg => 0b110_11010u8, // 20
    Element::Bh => 0b110_11110u8, // 21
    Element::Hs => 0b110_11111u8, // 22
    Element::Mt => 0b110_11101u8, // 23
    Element::Ds => 0b110_11100u8, // 24
    Element::Rg => 0b110_10100u8, // 25
    Element::Cn => 0b110_10101u8, // 26
    Element::Nh => 0b110_10111u8, // 27
    Element::Fl => 0b110_10110u8, // 28
    Element::Mc => 0b110_10010u8, // 29
    Element::Lv => 0b110_10011u8, // 30
    Element::Ts => 0b110_10001u8, // 31
    Element::Og => 0b110_10000u8, // 32
    }
});

static CHARGE_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // sign: 1 bit, absolute charge: 3 bits using Gray coding
    // sign: 1 bit, argument: 3 bits
    // negative charges are encoded by setting the sign bit, not as 1's or 2's complement
    Charge::Zero       => 0b0_000u8,
    Charge::PlusOne    => 0b0_001u8,
    Charge::PlusTwo    => 0b0_011u8,
    Charge::PlusThree  => 0b0_010u8,
    Charge::PlusFour   => 0b0_110u8,
    Charge::MinusOne   => 0b1_001u8,
    Charge::MinusTwo   => 0b1_011u8,
    Charge::MinusThree => 0b1_010u8,
    Charge::MinusFour  => 0b1_110u8,
}});

static EL_CONF_CODES = Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // 4 bits, last bit is always 0 for singlets (also some doublets)
    Conf::Zero         => 0b000_0u8,
    Conf::OneDoublet   => 0b000_1u8,
    Conf::TwoTriplet   => 0b001_1u8,
    Conf::TwoSinglet   => 0b001_0u8,
    Conf::ThreeQuartet => 0b010_1u8,
    Conf::ThreeDoublet => 0b010_0u8,
    Conf::FourTriplet  => 0b011_1u8,
    Conf::FourSinglet  => 0b011_0u8,
    Conf::FiveQuartet  => 0b100_1u8,
    Conf::FiveDoublet  => 0b100_0u8,
    Conf::SixTriplet   => 0b101_1u8,
    Conf::SixSinglet   => 0b101_0u8,
    Conf::SevenQuartet => 0b110_1u8,
    Conf::SevenDoublet => 0b110_0u8,
    Conf::EightTriplet => 0b111_1u8,
    Conf::EightSinglet => 0b111_0u8,
}});

static BONDING_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // >= quadruple: 2 bits, double and triple: 2 bits, single: 4 bits using Gray coding
    Bonding::Zero                    => 0b00_00_0000u8,
    Bonding::Single                  => 0b00_00_0001u8,
    Bonding::Double                  => 0b00_01_0000u8,
    Bonding::TwoSingle               => 0b00_00_0011u8,
    Bonding::Triple                  => 0b00_10_0000u8,
    Bonding::DoubleAndSingle         => 0b00_01_0001u8,
    Bonding::ThreeSingle             => 0b00_00_0010u8,
    Bonding::Quadruple               => 0b01_00_0000u8,
    Bonding::TripleAndSingle         => 0b00_10_0001u8,
    Bonding::TwoDouble               => 0b00_11_0000u8,
    Bonding::DoubleAndTwoSingle      => 0b00_01_0011u8,
    Bonding::FourSingle              => 0b00_00_0110u8,
    Bonding::Quintuple               => 0b10_00_0000u8,
    Bonding::QuadrupleAndSingle      => 0b01_00_0001u8,
    Bonding::TripleAndTwoSingle      => 0b00_10_0011u8,
    Bonding::DoubleAndThreeSingle    => 0b00_01_0010u8,
    Bonding::FiveSingle              => 0b00_00_0111u8,
    Bonding::Sextuple                => 0b11_00_0000u8,
    Bonding::QuintupleAndSingle      => 0b10_00_0001u8,
    Bonding::QuadrupleAndTwoSingle   => 0b01_00_0011u8,
    Bonding::TripleAndThreeSingle    => 0b00_10_0010u8,
    Bonding::DoubleAndFourSingle     => 0b00_01_0110u8,
    Bonding::SixSingle               => 0b00_00_0101u8,
    Bonding::QuadrupleAndThreeSingle => 0b01_00_0010u8,
    Bonding::TripleAndFourSingle     => 0b00_10_0110u8,
    Bonding::DoubleAndFiveSingle     => 0b00_01_0111u8,
    Bonding::SevenSingle             => 0b00_00_0100u8,
    Bonding::QuadrupleAndFourSingle  => 0b01_00_0110u8,
    Bonding::TripleAndFiveSingle     => 0b00_10_0111u8,
    Bonding::DoubleAndSixSingle      => 0b00_01_0101u8,
    Bonding::EightSingle             => 0b00_00_1100u8,
 }});

static RINGS_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // # of rings: 3, ring sizes: 5 bits
    // using Gray coding for number of rings
    Rings::Zero          => 0b000_00000u8,
    Rings::Three         => 0b001_00000u8,
    Rings::Four          => 0b001_00001u8,
    Rings::Five          => 0b001_00011u8,
    Rings::Six           => 0b001_00010u8,
    Rings::Seven         => 0b001_00110u8,
    Rings::Eight         => 0b001_00111u8,
    Rings::Nine          => 0b001_00101u8,
    Rings::Ten           => 0b001_00100u8,
    Rings::Eleven        => 0b001_01100u8,
    Rings::Twelve        => 0b001_01101u8,
    Rings::Thirteen      => 0b001_01111u8,
    Rings::Fourteen      => 0b001_01110u8,
    Rings::Fifteen       => 0b001_01010u8,
    Rings::Sixteen       => 0b001_01011u8,

    // Ordering by sum of ring sizes, then lexicographic, using Gray coding
    Rings::TwoThree       => 0b011_00000u8, // 6
    Rings::FourAndThree   => 0b011_00001u8, // 7
    Rings::TwoFour        => 0b011_00011u8, // 8
    Rings::FiveAndThree   => 0b011_00010u8,
    Rings::FiveAndFour    => 0b011_00110u8, // 9
    Rings::SixAndThree    => 0b011_00111u8, 
    Rings::TwoFive        => 0b011_00101u8, // 10
    Rings::SixAndFour     => 0b011_00100u8,
    Rings::SevenAndThree  => 0b011_01100u8,
    Rings::SixAndFive     => 0b011_01101u8, // 11
    Rings::SevenAndFour   => 0b011_01111u8,
    Rings::EightAndThree  => 0b011_01110u8,
    Rings::TwoSix         => 0b011_01010u8, // 12
    Rings::SevenAndFive   => 0b011_01011u8,
    Rings::EightAndFour   => 0b011_01001u8,
    Rings::SevenAndSix    => 0b011_01000u8, // 13
    Rings::EightAndFive   => 0b011_11000u8,
    Rings::TwoSeven       => 0b011_11001u8, // 14
    Rings::EightAndSix    => 0b011_11011u8,
    Rings::NineAndFive    => 0b011_11010u8,
    Rings::EightAndSeven  => 0b011_11110u8, // 15
    Rings::NineAndSix     => 0b011_11111u8,
    Rings::TenAndFive     => 0b011_11101u8,
    Rings::TwoEight       => 0b011_11100u8, // 16
    Rings::NineAndSeven   => 0b011_10100u8,

    // Only six-membered rings, use Gray coding for number of rings
    Rings::ThreeSix      => 0b010_00000u8,
    Rings::FourSix       => 0b010_00001u8,
    Rings::FiveSix       => 0b010_00011u8,
    Rings::SixSix        => 0b010_00010u8,
}});

```
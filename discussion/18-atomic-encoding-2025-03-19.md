# Encoding of atomic properties

* Element: 1 byte
* Valence: 1 byte (charge + electronic conf.)
* Bonding: 1 byte
* 0b1111_1111 indicates information unavailable (or outside coding domain)

```rust
use map_marco::hash_map;
use element::Element;

static ELEMENT_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // period: 3 bits, group: 5 bits (using 32-group layout instead of traditional 18-group)
    Element::H  => 0b000_00000u8,
    Element::He => 0b000_11111u8,
    Element::Li => 0b001_00000u8,
    Element::Be => 0b001_00001u8,
    Element::B  => 0b001_11010u8,
    Element::C  => 0b001_11011u8,
    Element::N  => 0b001_11100u8,
    Element::O  => 0b001_11101u8,
    Element::F  => 0b001_11110u8,
    Element::Ne => 0b001_11111u8,
    Element::Na => 0b010_00000u8,
    Element::Mg => 0b010_00001u8,
    Element::Al => 0b010_11010u8,
    Element::Si => 0b010_11011u8,
    Element::P  => 0b010_11100u8,
    Element::S  => 0b010_11101u8,
    Element::Cl => 0b010_11110u8,
    Element::Ar => 0b010_11111u8,
    Element::K  => 0b011_00000u8,
    Element::Ca => 0b011_00001u8,
    Element::Sc => 0b011_10000u8
    Element::Ti => 0b011_10001u8,
    Element::V  => 0b011_10010u8,
    Element::Cr => 0b011_10011u8,
    Element::Mn => 0b011_10100u8,
    Element::Fe => 0b011_10101u8,
    Element::Co => 0b011_10110u8,
    Element::Ni => 0b011_10111u8
    Element::Cu => 0b011_11000u8
    Element::Zn => 0b011_11001u8,
    Element::Ga => 0b011_11010u8,
    Element::Ge => 0b011_11011u8,
    Element::As => 0b011_11100u8,
    Element::Se => 0b011_11101u8,
    Element::Br => 0b011_11110u8,
    Element::Kr => 0b011_11111u8,
    Element::Rb => 0b100_00000u8,
    Element::Sr => 0b100_00001u8,
    Element::Y  => 0b100_10000u8,
    Element::Zr => 0b100_10001u8,
    Element::Nb => 0b100_10010u8,
    Element::Mo => 0b100_10011u8,
    Element::Tc => 0b100_10100u8,
    Element::Ru => 0b100_10101u8,
    Element::Rh => 0b100_10110u8,
    Element::Pd => 0b100_10111u8,
    Element::Ag => 0b100_11000u8,
    Element::Cd => 0b100_11001u8,
    Element::In => 0b100_11010u8,
    Element::Sn => 0b100_11011u8,
    Element::Sb => 0b100_11100u8,
    Element::Te => 0b100_11101u8,
    Element::I  => 0b100_11110u8,
    Element::Xe => 0b100_11111u8,
    Element::Cs => 0b101_00000u8,
    Element::Ba => 0b101_00001u8,
    Element::La => 0b101_00010u8,
    Element::Ce => 0b101_00011u8,
    Element::Pr => 0b101_00100u8,
    Element::Nd => 0b101_00101u8,
    Element::Pm => 0b101_00110u8,
    Element::Sm => 0b101_00111u8,
    Element::Eu => 0b101_01000u8,
    Element::Gd => 0b101_01001u8,
    Element::Tb => 0b101_01010u8,
    Element::Dy => 0b101_01011u8,
    Element::Ho => 0b101_01100u8,
    Element::Er => 0b101_01101u8,
    Element::Tm => 0b101_01110u8,
    Element::Yb => 0b101_01111u8,
    Element::Lu => 0b101_10000u8,
    Element::Hf => 0b101_10001u8,
    Element::Ta => 0b101_10010u8,
    Element::W  => 0b101_10011u8,
    Element::Re => 0b101_10100u8,
    Element::Os => 0b101_10101u8,
    Element::Ir => 0b101_10110u8,
    Element::Pt => 0b101_10111u8,
    Element::Au => 0b101_11000u8,
    Element::Hg => 0b101_11001u8,
    Element::Tl => 0b101_11010u8,
    Element::Pb => 0b101_11011u8,
    Element::Bi => 0b101_11100u8,
    Element::Po => 0b101_11101u8,
    Element::At => 0b101_11110u8,
    Element::Rn => 0b101_11111u8,
    Element::Fr => 0b110_00000u8,
    Element::Ra => 0b110_00001u8,
    Element::Ac => 0b110_00010u8,
    Element::Th => 0b110_00011u8,
    Element::Pa => 0b110_00100u8,
    Element::U  => 0b110_00101u8,
    Element::Np => 0b110_00110u8,
    Element::Pu => 0b110_00111u8,
    Element::Am => 0b110_01000u8,
    Element::Cm => 0b110_01001u8,
    Element::Bk => 0b110_01010u8,
    Element::Cf => 0b110_01011u8,
    Element::Es => 0b110_01100u8,
    Element::Fm => 0b110_01101u8,
    Element::Md => 0b110_01110u8,
    Element::No => 0b110_01111u8,
    Element::Lr => 0b110_10000u8,
    Element::Rf => 0b110_10001u8,
    Element::Db => 0b110_10010u8,
    Element::Sg => 0b110_10011u8,
    Element::Bh => 0b110_10100u8,
    Element::Hs => 0b110_10101u8,
    Element::Mt => 0b110_10110u8,
    Element::Ds => 0b110_10111u8,
    Element::Rg => 0b110_11000u8,
    Element::Cn => 0b110_11001u8,
    Element::Nh => 0b110_11010u8,
    Element::Fl => 0b110_11011u8,
    Element::Mc => 0b110_11100u8,
    Element::Lv => 0b110_11101u8,
    Element::Ts => 0b110_11110u8,
    Element::Og => 0b110_11111u8,
    }
});


static CHARGE_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // sign: 1 bit, absolute charge: 3 bits
    // Uses 1's complement for negative charges
    // TODO: Consider Gray coding
    Charge::Zero       => 0b0_000,
    Charge::PlusOne    => 0b0_001,
    Charge::PlusTwo    => 0b0_010,
    Charge::PlusThree  => 0b0_011,
    Charge::PlusFour   => 0b0_100,
    Charge::MinusOne   => 0b1_110,
    Charge::MinusTwo   => 0b1_101,
    Charge::MinusThree => 0b1_100,
    Charge::MinusFour  => 0b1_011
}});

static EL_CONF_CODES = Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // 4 bits, last bit is always 0 for singlets (also some doublets)
    Conf::Zero         => 0b000_0,
    Conf::OneDoublet   => 0b000_1,
    Conf::TwoTriplet   => 0b001_1,
    Conf::TwoSinglet   => 0b001_0,
    Conf::ThreeQuartet => 0b010_1,
    Conf::ThreeDoublet => 0b010_0,
    Conf::FourTriplet  => 0b011_1,
    Conf::FourSinglet  => 0b011_0,
    Conf::FiveQuartet  => 0b100_1,
    Conf::FiveDoublet  => 0b100_0,
    Conf::SixTriplet   => 0b101_1,
    Conf::SixSinglet   => 0b101_0,
    Conf::SevenQuartet => 0b110_1,
    Conf::SevenDoublet => 0b110_0,
    Conf::EightTriplet => 0b111_1,
    Conf::EightSinglet => 0b111_0,
}});

static BONDING_CODES: Lazy<HashMap<Element, u8>> = Lazy({||hash_map! {
    // >= quadruple: 2 bits, double and triple: 2 bits, single: 4 bits
    Bonding::Zero                    => 0b00_00_0000,
    Bonding::Single                  => 0b00_00_0001,
    Bonding::Double                  => 0b00_01_0000,
    Bonding::TwoSingle               => 0b00_00_0010,
    Bonding::Triple                  => 0b00_10_0000,
    Bonding::DoubleAndSingle         => 0b00_01_0001,
    Bonding::ThreeSingle             => 0b00_00_0011,
    Bonding::Quadruple               => 0b01_00_0000,
    Bonding::TripleAndSingle         => 0b00_10_0001,
    Bonding::TwoDouble               => 0b00_11_0000,
    Bonding::DoubleAndTwoSingle      => 0b00_01_0010,
    Bonding::FourSingle              => 0b00_00_0100,
    Bonding::Quintuple               => 0b10_00_0000,
    Bonding::QuadrupleAndSingle      => 0b01_00_0001,
    Bonding::TripleAndTwoSingle      => 0b00_10_0010,
    Bonding::DoubleAndThreeSingle    => 0b00_01_0011,
    Bonding::FiveSingle              => 0b00_00_0101,
    Bonding::Sextuple                => 0b11_00_0000,
    Bonding::QuintupleAndSingle      => 0b10_00_0001,
    Bonding::QuadrupleAndTwoSingle   => 0b01_00_0010,
    Bonding::TripleAndThreeSingle    => 0b00_10_0011,
    Bonding::DoubleAndFourSingle     => 0b00_01_0100,
    Bonding::SixSingle               => 0b00_00_0110,
    Bonding::QuadrupleAndThreeSingle => 0b01_00_0011,
    Bonding::TripleAndFourSingle     => 0b00_10_0100,
    Bonding::DoubleAndFiveSingle     => 0b00_01_0101,
    Bonding::SevenSingle             => 0b00_00_0111,
    Bonding::QuadrupleAndFourSingle  => 0b01_00_0100,
    Bonding::TripleAndFiveSingle     => 0b00_10_0101,
    Bonding::DoubleAndSixSingle      => 0b00_01_0110,
    Bonding::EightSingle             => 0b00_00_1000,
 }});

static RINGS_CODES: Lazy<HashMap<Element, u8>> = Lazy{||hash_map! 
    // # of rings: 3 bits, ring sizes: 5 bits
    // TODO: Consider using Gray coding
    Rings::Zero          => 0b000_00000,
    Rings::Three         => 0b001_00000,
    Rings::Four          => 0b001_00001,
    Rings::Five          => 0b001_00010,
    Rings::Six           => 0b001_00011,
    Rings::Seven         => 0b001_00100,
    Rings::Eight         => 0b001_00101,
    Rings::Nine          => 0b001_00110,
    Rings::Ten           => 0b001_00111,
    Rings::Eleven        => 0b001_01000,
    Rings::Twelve        => 0b001_01001,
    Rings::Thirteen      => 0b001_01010,
    Rings::Fourteen      => 0b001_01011,
    Rings::Fifteen       => 0b001_01100,
    Rings::Sixteen       => 0b001_01101,

    // Ordering by sum of ring sizes, then lexicographic
    Rings::TwoThree       => 0b010_00000, // 6
    Rings::FourAndThree   => 0b010_00001, // 7
    Rings::TwoFour        => 0b010_00010, // 8
    Rings::FiveAndThree   => 0b010_00011,
    Rings::FiveAndFour    => 0b010_00100, // 9
    Rings::SixAndThree    => 0b010_00101, 
    Rings::TwoFive        => 0b010_00110, // 10
    Rings::SixAndFour     => 0b010_00111,
    Rings::SevenAndThree  => 0b010_01000,
    Rings::SixAndFive     => 0b010_01001, // 11
    Rings::SevenAndFour   => 0b010_01010,
    Rings::EightAndThree  => 0b010_01011,
    Rings::TwoSix         => 0b010_01100, // 12
    Rings::SevenAndFive   => 0b010_01101,
    Rings::EightAndFour   => 0b010_01110,
    Rings::SevenAndSix    => 0b010_01111, // 13
    Rings::EightAndFive   => 0b010_10000,
    Rings::TwoSeven       => 0b010_10001, // 14
    Rings::EightAndSix    => 0b010_10010,
    Rings::NineAndFive    => 0b010_10011,
    Rings::EightAndSeven  => 0b010_10100, // 15
    Rings::NineAndSix     => 0b010_10101,
    Rings::TenAndFive     => 0b010_10110,
    Rings::TwoEight       => 0b010_10111, // 16
    Rings::NineAndSeven   => 0b010_11000,

    // Only six-membered rings
    Rings::ThreeSix      => 0b011_00000,
    Rings::FourSix       => 0b100_00000,
    Rings::FiveSix       => 0b101_00000,
    Rings::SixSix        => 0b110_00000,
});

```
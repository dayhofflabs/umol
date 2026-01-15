# Comparison of MOL parsing benchmarks

## Revision 6e85fb01 (2026-01-14)

     Running benches/mol_parsing.rs (/Users/dr/.cargo-target/release/deps/mol_parsing-783530498285640c)
mol_parsing/counts/valid
                        time:   [65.762 ns 66.104 ns 66.659 ns]
Found 8 outliers among 100 measurements (8.00%)
  3 (3.00%) high mild
  5 (5.00%) high severe
mol_parsing/counts/invalid
                        time:   [44.629 ns 44.730 ns 44.849 ns]
Found 9 outliers among 100 measurements (9.00%)
  3 (3.00%) high mild
  6 (6.00%) high severe

mol_parsing/atom/len69  time:   [199.45 ns 199.58 ns 199.73 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/atom/len69_invalid
                        time:   [194.70 ns 194.88 ns 195.08 ns]
Found 8 outliers among 100 measurements (8.00%)
  5 (5.00%) high mild
  3 (3.00%) high severe
mol_parsing/atom/len51  time:   [167.10 ns 167.31 ns 167.51 ns]
Found 6 outliers among 100 measurements (6.00%)
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/atom/len51_invalid
                        time:   [161.47 ns 161.76 ns 162.16 ns]
Found 10 outliers among 100 measurements (10.00%)
  5 (5.00%) high mild
  5 (5.00%) high severe
mol_parsing/atom/len42  time:   [142.54 ns 143.62 ns 145.20 ns]
Found 8 outliers among 100 measurements (8.00%)
  3 (3.00%) high mild
  5 (5.00%) high severe
mol_parsing/atom/len42_invalid
                        time:   [128.86 ns 128.99 ns 129.14 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/atom/len39  time:   [123.27 ns 123.47 ns 123.71 ns]
Found 7 outliers among 100 measurements (7.00%)
  3 (3.00%) high mild
  4 (4.00%) high severe
mol_parsing/atom/len39_invalid
                        time:   [120.54 ns 122.47 ns 124.60 ns]
Found 20 outliers among 100 measurements (20.00%)
  9 (9.00%) high mild
  11 (11.00%) high severe
mol_parsing/atom/len36  time:   [118.92 ns 119.17 ns 119.47 ns]
Found 8 outliers among 100 measurements (8.00%)
  5 (5.00%) high mild
  3 (3.00%) high severe
mol_parsing/atom/len36_invalid
                        time:   [115.63 ns 115.89 ns 116.16 ns]
Found 11 outliers among 100 measurements (11.00%)
  5 (5.00%) high mild
  6 (6.00%) high severe
mol_parsing/atom/len34  time:   [101.10 ns 101.19 ns 101.31 ns]
Found 9 outliers among 100 measurements (9.00%)
  7 (7.00%) high mild
  2 (2.00%) high severe
mol_parsing/atom/len34_invalid
                        time:   [97.473 ns 97.539 ns 97.619 ns]
Found 9 outliers among 100 measurements (9.00%)
  6 (6.00%) high mild
  3 (3.00%) high severe

mol_parsing/extended_atom/len69_extended
                        time:   [161.83 ns 161.98 ns 162.13 ns]
Found 11 outliers among 100 measurements (11.00%)
  4 (4.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_atom/len69_extended_invalid
                        time:   [169.26 ns 169.42 ns 169.63 ns]
Found 8 outliers among 100 measurements (8.00%)
  2 (2.00%) high mild
  6 (6.00%) high severe
mol_parsing/extended_atom/len51_extended
                        time:   [134.10 ns 134.25 ns 134.42 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_atom/len51_extended_invalid
                        time:   [124.45 ns 124.62 ns 124.81 ns]
Found 7 outliers among 100 measurements (7.00%)
  5 (5.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_atom/len42_extended
                        time:   [108.50 ns 108.64 ns 108.80 ns]
Found 8 outliers among 100 measurements (8.00%)
  5 (5.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_atom/len42_extended_invalid
                        time:   [107.59 ns 108.03 ns 108.75 ns]
Found 10 outliers among 100 measurements (10.00%)
  4 (4.00%) high mild
  6 (6.00%) high severe
mol_parsing/extended_atom/len39_extended
                        time:   [99.937 ns 100.81 ns 102.44 ns]
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_atom/len39_extended_invalid
                        time:   [89.137 ns 89.284 ns 89.443 ns]
mol_parsing/extended_atom/len36_extended
                        time:   [95.907 ns 96.008 ns 96.123 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_atom/len36_extended_invalid
                        time:   [82.356 ns 82.562 ns 82.863 ns]
Found 10 outliers among 100 measurements (10.00%)
  4 (4.00%) high mild
  6 (6.00%) high severe
mol_parsing/extended_atom/len34_extended
                        time:   [93.175 ns 93.245 ns 93.327 ns]
Found 9 outliers among 100 measurements (9.00%)
  4 (4.00%) high mild
  5 (5.00%) high severe
mol_parsing/extended_atom/len34_extended_invalid
                        time:   [93.938 ns 94.092 ns 94.254 ns]
Found 2 outliers among 100 measurements (2.00%)
  1 (1.00%) high mild
  1 (1.00%) high severe

mol_parsing/bond/len21  time:   [63.055 ns 63.154 ns 63.254 ns]
Found 3 outliers among 100 measurements (3.00%)
  3 (3.00%) high mild
mol_parsing/bond/len21_invalid
                        time:   [42.023 ns 42.063 ns 42.109 ns]
Found 5 outliers among 100 measurements (5.00%)
  4 (4.00%) high mild
  1 (1.00%) high severe
mol_parsing/bond/len12  time:   [45.107 ns 45.138 ns 45.172 ns]
Found 11 outliers among 100 measurements (11.00%)
  3 (3.00%) high mild
  8 (8.00%) high severe
mol_parsing/bond/len12_invalid
                        time:   [41.984 ns 42.021 ns 42.063 ns]
Found 8 outliers among 100 measurements (8.00%)
  2 (2.00%) low mild
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/bond/len9   time:   [32.693 ns 32.726 ns 32.769 ns]
Found 12 outliers among 100 measurements (12.00%)
  4 (4.00%) high mild
  8 (8.00%) high severe
mol_parsing/bond/len9_invalid
                        time:   [30.875 ns 30.900 ns 30.931 ns]
Found 13 outliers among 100 measurements (13.00%)
  2 (2.00%) low mild
  3 (3.00%) high mild
  8 (8.00%) high severe

mol_parsing/extended_bond/len21_extended
                        time:   [62.877 ns 62.971 ns 63.120 ns]
Found 11 outliers among 100 measurements (11.00%)
  4 (4.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_bond/len21_extended_invalid
                        time:   [44.345 ns 44.391 ns 44.448 ns]
Found 8 outliers among 100 measurements (8.00%)
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_bond/len18_ring_extended
                        time:   [51.932 ns 51.975 ns 52.029 ns]
Found 10 outliers among 100 measurements (10.00%)
  4 (4.00%) high mild
  6 (6.00%) high severe
mol_parsing/extended_bond/len18_ring_extended_invalid
                        time:   [46.197 ns 46.463 ns 46.944 ns]
Found 6 outliers among 100 measurements (6.00%)
  1 (1.00%) high mild
  5 (5.00%) high severe
mol_parsing/extended_bond/len12
                        time:   [40.401 ns 40.459 ns 40.513 ns]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high severe
mol_parsing/extended_bond/len12_extended_invalid
                        time:   [32.779 ns 32.885 ns 33.058 ns]
Found 8 outliers among 100 measurements (8.00%)
  2 (2.00%) high mild
  6 (6.00%) high severe
mol_parsing/extended_bond/len9_extended
                        time:   [32.191 ns 32.222 ns 32.255 ns]
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_bond/len9_extended_invalid
                        time:   [27.804 ns 27.831 ns 27.861 ns]
Found 6 outliers among 100 measurements (6.00%)
  5 (5.00%) high mild
  1 (1.00%) high severe

mol_parsing/legacy_atom_list/no_exclusion
                        time:   [4.3177 ns 4.3291 ns 4.3401 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) low mild
  1 (1.00%) high severe
mol_parsing/legacy_atom_list/exclusion
                        time:   [4.3090 ns 4.3196 ns 4.3302 ns]
Found 2 outliers among 100 measurements (2.00%)
  1 (1.00%) low mild
  1 (1.00%) high severe


mol_parsing/properties/chg
                        time:   [54.743 ns 55.172 ns 55.937 ns]
Found 4 outliers among 100 measurements (4.00%)
  1 (1.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/chg6
                        time:   [165.85 ns 166.25 ns 166.67 ns]
Found 2 outliers among 100 measurements (2.00%)
  2 (2.00%) high mild
mol_parsing/properties/chg_invalid
                        time:   [54.726 ns 54.865 ns 55.005 ns]
Found 2 outliers among 100 measurements (2.00%)
  2 (2.00%) high mild
mol_parsing/properties/rad1
                        time:   [51.521 ns 51.938 ns 52.574 ns]
Found 14 outliers among 100 measurements (14.00%)
  12 (12.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/rad6
                        time:   [148.20 ns 148.41 ns 148.66 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/rad_invalid
                        time:   [51.646 ns 51.977 ns 52.488 ns]
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/iso1
                        time:   [49.972 ns 50.085 ns 50.207 ns]
mol_parsing/properties/iso8
                        time:   [178.32 ns 178.68 ns 179.12 ns]
Found 3 outliers among 100 measurements (3.00%)
  1 (1.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/iso_invalid
                        time:   [49.749 ns 49.809 ns 49.877 ns]
Found 8 outliers among 100 measurements (8.00%)
  1 (1.00%) high mild
  7 (7.00%) high severe
mol_parsing/properties/sty1
                        time:   [4.8231 ns 4.8680 ns 4.9538 ns]
Found 9 outliers among 100 measurements (9.00%)
  1 (1.00%) low mild
  3 (3.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/sty2
                        time:   [4.7802 ns 4.7871 ns 4.7945 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/slb1
                        time:   [4.8009 ns 4.8271 ns 4.8814 ns]
Found 6 outliers among 100 measurements (6.00%)
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/slb2
                        time:   [4.7905 ns 4.7966 ns 4.8029 ns]
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) low mild
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/sal_simple
                        time:   [4.8160 ns 4.8604 ns 4.9486 ns]
Found 3 outliers among 100 measurements (3.00%)
  1 (1.00%) low mild
  1 (1.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/sal_multi
                        time:   [4.8311 ns 4.8392 ns 4.8470 ns]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high mild
mol_parsing/properties/sbl_simple
                        time:   [4.7875 ns 4.7958 ns 4.8041 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/sbl_multi
                        time:   [4.7979 ns 4.8517 ns 4.9603 ns]
Found 5 outliers among 100 measurements (5.00%)
  2 (2.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/alias_simple
                        time:   [4.7082 ns 4.7404 ns 4.7926 ns]
Found 2 outliers among 100 measurements (2.00%)
  2 (2.00%) high severe
mol_parsing/properties/alias_long
                        time:   [4.7136 ns 4.7699 ns 4.8764 ns]
Found 2 outliers among 100 measurements (2.00%)
  2 (2.00%) high severe
mol_parsing/properties/value_simple
                        time:   [37.608 ns 38.026 ns 38.775 ns]
Found 3 outliers among 100 measurements (3.00%)
  3 (3.00%) high severe
mol_parsing/properties/value_long
                        time:   [39.430 ns 40.469 ns 41.727 ns]
Found 9 outliers among 100 measurements (9.00%)
  3 (3.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/sst1
                        time:   [4.7938 ns 4.8421 ns 4.9326 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/sst2
                        time:   [4.8384 ns 4.9073 ns 4.9981 ns]
Found 9 outliers among 100 measurements (9.00%)
  4 (4.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/smt_simple
                        time:   [4.8351 ns 4.8858 ns 4.9830 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/smt_long
                        time:   [4.8385 ns 4.8490 ns 4.8601 ns]
Found 10 outliers among 100 measurements (10.00%)
  2 (2.00%) high mild
  8 (8.00%) high severe
mol_parsing/properties/zbo1
                        time:   [5.0968 ns 5.1015 ns 5.1071 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/zbo4
                        time:   [5.1247 ns 5.1753 ns 5.2688 ns]
Found 8 outliers among 100 measurements (8.00%)
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/zch1
                        time:   [5.1526 ns 5.1590 ns 5.1657 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/zch4
                        time:   [5.1567 ns 5.1633 ns 5.1706 ns]
Found 7 outliers among 100 measurements (7.00%)
  6 (6.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/hyd1
                        time:   [4.9980 ns 5.0504 ns 5.1523 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/hyd6
                        time:   [5.0112 ns 5.0274 ns 5.0475 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe

mol_parsing/extended_properties/chg1_extended
                        time:   [54.297 ns 54.389 ns 54.480 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/chg8_extended
                        time:   [206.79 ns 209.20 ns 214.01 ns]
Found 7 outliers among 100 measurements (7.00%)
  1 (1.00%) low mild
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/rad1_extended
                        time:   [52.514 ns 52.659 ns 52.868 ns]
Found 9 outliers among 100 measurements (9.00%)
  4 (4.00%) high mild
  5 (5.00%) high severe
mol_parsing/extended_properties/iso1_extended
                        time:   [50.210 ns 50.426 ns 50.715 ns]
Found 9 outliers among 100 measurements (9.00%)
  4 (4.00%) high mild
  5 (5.00%) high severe
mol_parsing/extended_properties/iso8_extended
                        time:   [177.82 ns 178.98 ns 181.13 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/sty1_extended
                        time:   [47.902 ns 48.152 ns 48.626 ns]
Found 11 outliers among 100 measurements (11.00%)
  6 (6.00%) high mild
  5 (5.00%) high severe
mol_parsing/extended_properties/sty8_extended
                        time:   [146.30 ns 148.46 ns 152.28 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_properties/sst1_extended
                        time:   [46.107 ns 46.442 ns 47.171 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_properties/smt_multiplier_extended
                        time:   [61.416 ns 61.480 ns 61.550 ns]
Found 3 outliers among 100 measurements (3.00%)
  3 (3.00%) high severe
mol_parsing/extended_properties/smt_subscript_extended
                        time:   [55.346 ns 55.443 ns 55.552 ns]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high mild
mol_parsing/extended_properties/slb1_extended
                        time:   [47.518 ns 47.952 ns 48.793 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_properties/sal5_extended
                        time:   [76.500 ns 76.569 ns 76.656 ns]
Found 10 outliers among 100 measurements (10.00%)
  3 (3.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/sbl4_extended
                        time:   [68.382 ns 68.438 ns 68.506 ns]
Found 12 outliers among 100 measurements (12.00%)
  5 (5.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/alias_extended
                        time:   [4.9997 ns 5.0090 ns 5.0190 ns]
Found 3 outliers among 100 measurements (3.00%)
  1 (1.00%) low mild
  1 (1.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/value_extended
                        time:   [37.335 ns 37.359 ns 37.388 ns]
Found 15 outliers among 100 measurements (15.00%)
  8 (8.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/als1_extended
                        time:   [19.738 ns 19.759 ns 19.788 ns]
Found 11 outliers among 100 measurements (11.00%)
  3 (3.00%) high mild
  8 (8.00%) high severe
mol_parsing/extended_properties/als4_extended
                        time:   [19.721 ns 19.750 ns 19.790 ns]
Found 14 outliers among 100 measurements (14.00%)
  8 (8.00%) high mild
  6 (6.00%) high severe
mol_parsing/extended_properties/apo1_extended
                        time:   [52.152 ns 52.219 ns 52.321 ns]
Found 10 outliers among 100 measurements (10.00%)
  6 (6.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/aal1_extended
                        time:   [52.043 ns 52.088 ns 52.147 ns]
Found 13 outliers among 100 measurements (13.00%)
  4 (4.00%) high mild
  9 (9.00%) high severe
mol_parsing/extended_properties/aal2_extended
                        time:   [65.661 ns 65.750 ns 65.853 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_properties/rbc1_extended
                        time:   [56.331 ns 56.968 ns 58.130 ns]
Found 11 outliers among 100 measurements (11.00%)
  6 (6.00%) high mild
  5 (5.00%) high severe
mol_parsing/extended_properties/sub1_extended
                        time:   [55.246 ns 55.322 ns 55.421 ns]
Found 5 outliers among 100 measurements (5.00%)
  4 (4.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/uns1_extended
                        time:   [52.202 ns 52.235 ns 52.280 ns]
Found 9 outliers among 100 measurements (9.00%)
  1 (1.00%) low mild
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/lin1_extended
                        time:   [71.733 ns 71.849 ns 71.978 ns]
Found 7 outliers among 100 measurements (7.00%)
  3 (3.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/zbo1_extended
                        time:   [5.1743 ns 5.1827 ns 5.1918 ns]
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_properties/zbo4_extended
                        time:   [5.1914 ns 5.2174 ns 5.2552 ns]
Found 6 outliers among 100 measurements (6.00%)
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/zch1_extended
                        time:   [5.2186 ns 5.2234 ns 5.2291 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_properties/zch4_extended
                        time:   [5.2181 ns 5.2244 ns 5.2326 ns]
Found 10 outliers among 100 measurements (10.00%)
  6 (6.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/hyd1_extended
                        time:   [5.1649 ns 5.1708 ns 5.1777 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/hyd6_extended
                        time:   [5.1600 ns 5.1644 ns 5.1694 ns]
Found 4 outliers among 100 measurements (4.00%)
  1 (1.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_properties/scn1_extended
                        time:   [47.310 ns 47.346 ns 47.390 ns]
Found 8 outliers among 100 measurements (8.00%)
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/scn2_extended
                        time:   [57.645 ns 57.952 ns 58.512 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_properties/sds_extended
                        time:   [46.257 ns 46.874 ns 47.956 ns]
Found 11 outliers among 100 measurements (11.00%)
  1 (1.00%) high mild
  10 (10.00%) high severe
mol_parsing/extended_properties/spa12_extended
                        time:   [130.56 ns 130.70 ns 130.84 ns]
Found 7 outliers among 100 measurements (7.00%)
  5 (5.00%) high mild
  2 (2.00%) high severe
mol_parsing/extended_properties/crs_extended
                        time:   [72.566 ns 72.606 ns 72.649 ns]
Found 13 outliers among 100 measurements (13.00%)
  3 (3.00%) high mild
  10 (10.00%) high severe
mol_parsing/extended_properties/sdi1_extended
                        time:   [61.061 ns 61.102 ns 61.153 ns]
Found 14 outliers among 100 measurements (14.00%)
  7 (7.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/sdi2_extended
                        time:   [61.127 ns 61.208 ns 61.308 ns]
Found 13 outliers among 100 measurements (13.00%)
  9 (9.00%) high mild
  4 (4.00%) high severe
mol_parsing/extended_properties/sbv_extended
                        time:   [52.257 ns 52.600 ns 53.002 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/extended_properties/sdt_extended
                        time:   [5.0176 ns 5.0288 ns 5.0408 ns]
Found 6 outliers among 100 measurements (6.00%)
  1 (1.00%) low mild
  4 (4.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/sdd_extended
                        time:   [103.47 ns 103.67 ns 103.90 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/extended_properties/scd_extended
                        time:   [42.872 ns 42.955 ns 43.056 ns]
Found 14 outliers among 100 measurements (14.00%)
  5 (5.00%) high mild
  9 (9.00%) high severe
mol_parsing/extended_properties/sed_extended
                        time:   [45.093 ns 45.202 ns 45.326 ns]
Found 11 outliers among 100 measurements (11.00%)
  3 (3.00%) high mild
  8 (8.00%) high severe
mol_parsing/extended_properties/spl1_extended
                        time:   [51.228 ns 51.291 ns 51.366 ns]
Found 11 outliers among 100 measurements (11.00%)
  4 (4.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/snc1_extended
                        time:   [51.710 ns 51.816 ns 51.942 ns]
Found 12 outliers among 100 measurements (12.00%)
  5 (5.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/rgp1_extended
                        time:   [50.988 ns 51.179 ns 51.502 ns]
Found 13 outliers among 100 measurements (13.00%)
  6 (6.00%) high mild
  7 (7.00%) high severe
mol_parsing/extended_properties/log1_extended
                        time:   [175.31 ns 176.47 ns 178.05 ns]
Found 5 outliers among 100 measurements (5.00%)
  4 (4.00%) high mild
  1 (1.00%) high severe

## Revision 649d6120 (2025-09-05)

     Running benches/parsing_bench.rs (/Users/dr/.cargo-target/release/deps/parsing_bench-81be5e344550426e)
mol_parsing/counts/valid
                        time:   [71.883 ns 72.111 ns 72.344 ns]
                        change: [+8.3348% +8.9695% +9.5100%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) low mild
  3 (3.00%) high mild
mol_parsing/counts/invalid
                        time:   [71.237 ns 71.379 ns 71.518 ns]
                        change: [+59.300% +59.864% +60.441%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 5 outliers among 100 measurements (5.00%)
  2 (2.00%) high mild
  3 (3.00%) high severe

mol_parsing/atom/len69_basic
                        time:   [134.37 ns 135.34 ns 136.83 ns]
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe
mol_parsing/atom/len69_basic_invalid
                        time:   [133.18 ns 133.39 ns 133.63 ns]
Found 15 outliers among 100 measurements (15.00%)
  2 (2.00%) low mild
  9 (9.00%) high mild
  4 (4.00%) high severe
mol_parsing/atom/len51_basic
                        time:   [120.52 ns 120.78 ns 121.07 ns]
Found 11 outliers among 100 measurements (11.00%)
  4 (4.00%) low mild
  7 (7.00%) high mild
mol_parsing/atom/len51_basic_invalid
                        time:   [117.47 ns 117.72 ns 117.99 ns]
Found 2 outliers among 100 measurements (2.00%)
  2 (2.00%) high mild
mol_parsing/atom/len42_basic
                        time:   [108.40 ns 108.63 ns 108.87 ns]
Found 9 outliers among 100 measurements (9.00%)
  8 (8.00%) high mild
  1 (1.00%) high severe
mol_parsing/atom/len42_basic_invalid
                        time:   [109.80 ns 110.10 ns 110.42 ns]
Found 7 outliers among 100 measurements (7.00%)
  7 (7.00%) high mild
mol_parsing/atom/len39_basic
                        time:   [97.723 ns 97.990 ns 98.331 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/atom/len39_basic_invalid
                        time:   [87.622 ns 88.919 ns 91.173 ns]
Found 5 outliers among 100 measurements (5.00%)
  2 (2.00%) high mild
  3 (3.00%) high severe
mol_parsing/atom/len36_basic
                        time:   [88.886 ns 89.081 ns 89.280 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/atom/len36_basic_invalid
                        time:   [79.763 ns 80.090 ns 80.572 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/atom/len34_basic
                        time:   [77.321 ns 77.460 ns 77.609 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/atom/len34_basic_invalid
                        time:   [72.547 ns 72.689 ns 72.844 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe

mol_parsing/atomlike/len69
                        time:   [174.88 ns 175.24 ns 175.66 ns]
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
mol_parsing/atomlike/len69_invalid
                        time:   [168.21 ns 170.18 ns 174.05 ns]
Found 11 outliers among 100 measurements (11.00%)
  6 (6.00%) high mild
  5 (5.00%) high severe
mol_parsing/atomlike/len51
                        time:   [145.17 ns 145.44 ns 145.75 ns]
Found 5 outliers among 100 measurements (5.00%)
  4 (4.00%) high mild
  1 (1.00%) high severe
mol_parsing/atomlike/len51_invalid
                        time:   [127.30 ns 127.41 ns 127.54 ns]
Found 8 outliers among 100 measurements (8.00%)
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/atomlike/len42
                        time:   [115.81 ns 116.09 ns 116.39 ns]
Found 8 outliers among 100 measurements (8.00%)
  7 (7.00%) high mild
  1 (1.00%) high severe
mol_parsing/atomlike/len42_invalid
                        time:   [119.62 ns 124.63 ns 130.66 ns]
Found 6 outliers among 100 measurements (6.00%)
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/atomlike/len39
                        time:   [130.26 ns 135.36 ns 141.25 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/atomlike/len39_invalid
                        time:   [89.787 ns 89.938 ns 90.118 ns]
Found 14 outliers among 100 measurements (14.00%)
  1 (1.00%) low severe
  4 (4.00%) high mild
  9 (9.00%) high severe
mol_parsing/atomlike/len36
                        time:   [100.91 ns 100.96 ns 101.02 ns]
Found 7 outliers among 100 measurements (7.00%)
  1 (1.00%) low mild
  2 (2.00%) high mild
  4 (4.00%) high severe
mol_parsing/atomlike/len36_invalid
                        time:   [82.239 ns 82.315 ns 82.392 ns]
Found 7 outliers among 100 measurements (7.00%)
  2 (2.00%) low mild
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/atomlike/len34
                        time:   [99.408 ns 99.507 ns 99.607 ns]
Found 9 outliers among 100 measurements (9.00%)
  2 (2.00%) low severe
  3 (3.00%) low mild
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/atomlike/len34_invalid
                        time:   [99.487 ns 99.615 ns 99.778 ns]
Found 12 outliers among 100 measurements (12.00%)
  3 (3.00%) low mild
  5 (5.00%) high mild
  4 (4.00%) high severe

mol_parsing/bond/len21_basic
                        time:   [44.521 ns 44.572 ns 44.627 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/bond/len21_basic_invalid
                        time:   [36.458 ns 36.543 ns 36.625 ns]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high mild
mol_parsing/bond/len12_basic
                        time:   [44.000 ns 44.226 ns 44.610 ns]
Found 4 outliers among 100 measurements (4.00%)
  1 (1.00%) high mild
  3 (3.00%) high severe
mol_parsing/bond/len12_basic_invalid
                        time:   [37.357 ns 37.427 ns 37.497 ns]
Found 8 outliers among 100 measurements (8.00%)
  1 (1.00%) low severe
  3 (3.00%) low mild
  4 (4.00%) high mild
mol_parsing/bond/len9_basic
                        time:   [32.692 ns 32.722 ns 32.755 ns]
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe
mol_parsing/bond/len9_basic_invalid
                        time:   [26.533 ns 26.567 ns 26.608 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe

mol_parsing/bondlike/len21
                        time:   [68.539 ns 68.627 ns 68.722 ns]
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
mol_parsing/bondlike/len21_invalid
                        time:   [86.722 ns 86.816 ns 86.916 ns]
Found 9 outliers among 100 measurements (9.00%)
  7 (7.00%) high mild
  2 (2.00%) high severe
mol_parsing/bondlike/len18_ring
                        time:   [55.290 ns 55.354 ns 55.436 ns]
Found 10 outliers among 100 measurements (10.00%)
  3 (3.00%) high mild
  7 (7.00%) high severe
mol_parsing/bondlike/len18_invalid
                        time:   [55.698 ns 55.744 ns 55.795 ns]
Found 6 outliers among 100 measurements (6.00%)
  1 (1.00%) high mild
  5 (5.00%) high severe
mol_parsing/bondlike/len12
                        time:   [42.511 ns 42.551 ns 42.593 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/bondlike/len12_invalid
                        time:   [34.302 ns 34.349 ns 34.409 ns]
Found 12 outliers among 100 measurements (12.00%)
  6 (6.00%) high mild
  6 (6.00%) high severe
mol_parsing/bondlike/len9
                        time:   [34.258 ns 34.495 ns 34.974 ns]
Found 6 outliers among 100 measurements (6.00%)
  1 (1.00%) high mild
  5 (5.00%) high severe
mol_parsing/bondlike/len9_invalid
                        time:   [35.328 ns 35.356 ns 35.388 ns]
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe

mol_parsing/legacy_atom_list/no_exclusion
                        time:   [77.575 ns 77.663 ns 77.760 ns]
                        change: [+1686.9% +1692.9% +1698.7%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
mol_parsing/legacy_atom_list/exclusion
                        time:   [77.541 ns 77.618 ns 77.698 ns]
                        change: [+1690.7% +1696.5% +1702.3%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe

mol_parsing/basic_properties/chg_basic
                        time:   [59.115 ns 59.162 ns 59.214 ns]
Found 9 outliers among 100 measurements (9.00%)
  3 (3.00%) high mild
  6 (6.00%) high severe
mol_parsing/basic_properties/chg6_basic
                        time:   [164.99 ns 165.23 ns 165.51 ns]
Found 8 outliers among 100 measurements (8.00%)
  5 (5.00%) high mild
  3 (3.00%) high severe
mol_parsing/basic_properties/chg_basic_invalid
                        time:   [59.515 ns 59.574 ns 59.633 ns]
Found 3 outliers among 100 measurements (3.00%)
  1 (1.00%) high mild
  2 (2.00%) high severe
mol_parsing/basic_properties/rad1_basic
                        time:   [55.676 ns 55.899 ns 56.347 ns]
Found 7 outliers among 100 measurements (7.00%)
  2 (2.00%) high mild
  5 (5.00%) high severe
mol_parsing/basic_properties/rad6_basic
                        time:   [153.92 ns 154.04 ns 154.19 ns]
Found 9 outliers among 100 measurements (9.00%)
  5 (5.00%) high mild
  4 (4.00%) high severe
mol_parsing/basic_properties/rad_basic_invalid
                        time:   [53.580 ns 54.114 ns 55.205 ns]
Found 12 outliers among 100 measurements (12.00%)
  4 (4.00%) high mild
  8 (8.00%) high severe
mol_parsing/basic_properties/iso1_basic
                        time:   [53.295 ns 53.825 ns 54.916 ns]
Found 11 outliers among 100 measurements (11.00%)
  6 (6.00%) high mild
  5 (5.00%) high severe
mol_parsing/basic_properties/iso8_basic
                        time:   [168.76 ns 170.49 ns 173.88 ns]
Found 8 outliers among 100 measurements (8.00%)
  5 (5.00%) high mild
  3 (3.00%) high severe
mol_parsing/basic_properties/iso_basic_invalid
                        time:   [53.259 ns 53.317 ns 53.385 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/basic_properties/sty1_basic
                        time:   [51.220 ns 51.384 ns 51.613 ns]
Found 9 outliers among 100 measurements (9.00%)
  3 (3.00%) high mild
  6 (6.00%) high severe
mol_parsing/basic_properties/sty2_basic
                        time:   [67.218 ns 67.996 ns 69.588 ns]
Found 16 outliers among 100 measurements (16.00%)
  2 (2.00%) high mild
  14 (14.00%) high severe
mol_parsing/basic_properties/slb1_basic
                        time:   [57.241 ns 57.826 ns 58.651 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/basic_properties/slb2_basic
                        time:   [73.431 ns 73.815 ns 74.583 ns]
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/basic_properties/sal_simple_basic
                        time:   [50.982 ns 51.253 ns 51.734 ns]
Found 9 outliers among 100 measurements (9.00%)
  4 (4.00%) high mild
  5 (5.00%) high severe
mol_parsing/basic_properties/sal_multi_basic
                        time:   [66.892 ns 66.964 ns 67.045 ns]
Found 7 outliers among 100 measurements (7.00%)
  1 (1.00%) low mild
  1 (1.00%) high mild
  5 (5.00%) high severe
mol_parsing/basic_properties/sbl_simple_basic
                        time:   [51.020 ns 51.566 ns 52.736 ns]
Found 5 outliers among 100 measurements (5.00%)
  2 (2.00%) high mild
  3 (3.00%) high severe
mol_parsing/basic_properties/sbl_multi_basic
                        time:   [59.220 ns 59.871 ns 61.119 ns]
Found 1 outliers among 100 measurements (1.00%)
  1 (1.00%) high severe
mol_parsing/basic_properties/alias_simple_basic
                        time:   [17.830 ns 18.130 ns 18.528 ns]
Found 8 outliers among 100 measurements (8.00%)
  5 (5.00%) high mild
  3 (3.00%) high severe
mol_parsing/basic_properties/alias_long_basic
                        time:   [16.665 ns 16.915 ns 17.262 ns]
Found 8 outliers among 100 measurements (8.00%)
  1 (1.00%) high mild
  7 (7.00%) high severe
mol_parsing/basic_properties/value_simple_basic
                        time:   [42.019 ns 42.430 ns 43.285 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/basic_properties/value_long_basic
                        time:   [45.369 ns 45.740 ns 46.338 ns]
Found 7 outliers among 100 measurements (7.00%)
  3 (3.00%) high mild
  4 (4.00%) high severe
mol_parsing/basic_properties/sst1_basic
                        time:   [52.474 ns 52.962 ns 53.987 ns]
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) high mild
  4 (4.00%) high severe
mol_parsing/basic_properties/sst2_basic
                        time:   [65.465 ns 66.189 ns 67.632 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/basic_properties/smt_simple_basic
                        time:   [88.751 ns 89.129 ns 89.786 ns]
Found 12 outliers among 100 measurements (12.00%)
  8 (8.00%) high mild
  4 (4.00%) high severe
mol_parsing/basic_properties/smt_long_basic
                        time:   [81.281 ns 82.160 ns 83.728 ns]
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe
mol_parsing/basic_properties/zbo1_basic
                        time:   [54.988 ns 55.315 ns 55.979 ns]
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) high mild
  4 (4.00%) high severe
mol_parsing/basic_properties/zbo4_basic
                        time:   [114.70 ns 115.74 ns 117.59 ns]
Found 5 outliers among 100 measurements (5.00%)
  1 (1.00%) high mild
  4 (4.00%) high severe
mol_parsing/basic_properties/zch1_basic
                        time:   [57.044 ns 57.323 ns 57.704 ns]
Found 13 outliers among 100 measurements (13.00%)
  3 (3.00%) high mild
  10 (10.00%) high severe
mol_parsing/basic_properties/zch4_basic
                        time:   [113.42 ns 114.39 ns 116.47 ns]
Found 10 outliers among 100 measurements (10.00%)
  2 (2.00%) high mild
  8 (8.00%) high severe
mol_parsing/basic_properties/hyd1_basic
                        time:   [54.635 ns 55.363 ns 56.886 ns]
Found 19 outliers among 100 measurements (19.00%)
  19 (19.00%) high severe
mol_parsing/basic_properties/hyd6_basic
                        time:   [148.01 ns 148.82 ns 150.30 ns]
Found 10 outliers among 100 measurements (10.00%)
  3 (3.00%) high mild
  7 (7.00%) high severe

mol_parsing/properties/chg1
                        time:   [59.356 ns 59.652 ns 60.239 ns]
Found 12 outliers among 100 measurements (12.00%)
  1 (1.00%) low mild
  6 (6.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/chg8
                        time:   [206.04 ns 208.21 ns 212.40 ns]
Found 11 outliers among 100 measurements (11.00%)
  8 (8.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/rad1
                        time:   [56.431 ns 56.506 ns 56.604 ns]
                        change: [+8.7195% +12.372% +18.164%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 11 outliers among 100 measurements (11.00%)
  4 (4.00%) high mild
  7 (7.00%) high severe
mol_parsing/properties/iso1
                        time:   [55.491 ns 56.156 ns 57.279 ns]
                        change: [+10.567% +11.464% +12.762%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 12 outliers among 100 measurements (12.00%)
  6 (6.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/iso8
                        time:   [169.27 ns 171.03 ns 173.66 ns]
                        change: [−5.0883% −4.6060% −3.9804%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 15 outliers among 100 measurements (15.00%)
  4 (4.00%) high mild
  11 (11.00%) high severe
mol_parsing/properties/sty1
                        time:   [52.381 ns 52.801 ns 53.663 ns]
                        change: [+965.46% +981.48% +995.52%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 12 outliers among 100 measurements (12.00%)
  6 (6.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/sty8
                        time:   [150.70 ns 151.83 ns 154.23 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/sst1
                        time:   [53.204 ns 53.711 ns 54.763 ns]
                        change: [+996.21% +1009.1% +1023.2%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 8 outliers among 100 measurements (8.00%)
  3 (3.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/smt_multiplier
                        time:   [90.514 ns 91.009 ns 92.085 ns]
Found 6 outliers among 100 measurements (6.00%)
  1 (1.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/smt_subscript
                        time:   [81.940 ns 82.497 ns 83.558 ns]
Found 11 outliers among 100 measurements (11.00%)
  3 (3.00%) high mild
  8 (8.00%) high severe
mol_parsing/properties/slb1
                        time:   [49.849 ns 50.254 ns 51.074 ns]
                        change: [+912.46% +933.04% +948.20%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 10 outliers among 100 measurements (10.00%)
  3 (3.00%) high mild
  7 (7.00%) high severe
mol_parsing/properties/sal5
                        time:   [83.040 ns 83.785 ns 85.300 ns]
Found 14 outliers among 100 measurements (14.00%)
  6 (6.00%) high mild
  8 (8.00%) high severe
mol_parsing/properties/sbl4
                        time:   [74.732 ns 74.807 ns 74.895 ns]
Found 10 outliers among 100 measurements (10.00%)
  5 (5.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/alias
                        time:   [17.753 ns 17.769 ns 17.786 ns]
Found 12 outliers among 100 measurements (12.00%)
  5 (5.00%) high mild
  7 (7.00%) high severe
mol_parsing/properties/value
                        time:   [42.722 ns 42.774 ns 42.832 ns]
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/als1
                        time:   [25.276 ns 25.517 ns 25.972 ns]
Found 8 outliers among 100 measurements (8.00%)
  6 (6.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/als4
                        time:   [25.016 ns 25.042 ns 25.074 ns]
Found 10 outliers among 100 measurements (10.00%)
  5 (5.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/apo1
                        time:   [56.574 ns 56.632 ns 56.705 ns]
Found 10 outliers among 100 measurements (10.00%)
  3 (3.00%) high mild
  7 (7.00%) high severe
mol_parsing/properties/aal1
                        time:   [57.526 ns 57.624 ns 57.736 ns]
Found 11 outliers among 100 measurements (11.00%)
  2 (2.00%) high mild
  9 (9.00%) high severe
mol_parsing/properties/aal2
                        time:   [72.046 ns 72.177 ns 72.320 ns]
Found 11 outliers among 100 measurements (11.00%)
  7 (7.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/rbc1
                        time:   [61.258 ns 61.383 ns 61.526 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/sub1
                        time:   [59.457 ns 59.806 ns 60.479 ns]
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/uns1
                        time:   [56.441 ns 56.972 ns 58.053 ns]
Found 8 outliers among 100 measurements (8.00%)
  6 (6.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/lin1
                        time:   [78.219 ns 78.347 ns 78.489 ns]
Found 8 outliers among 100 measurements (8.00%)
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/zbo1
                        time:   [56.067 ns 56.143 ns 56.231 ns]
                        change: [+1001.2% +1003.8% +1006.9%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/zbo4
                        time:   [115.33 ns 115.57 ns 115.82 ns]
                        change: [+2121.4% +2147.0% +2164.1%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 2 outliers among 100 measurements (2.00%)
  1 (1.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/zch1
                        time:   [58.786 ns 58.863 ns 58.989 ns]
                        change: [+1040.8% +1043.1% +1045.6%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 14 outliers among 100 measurements (14.00%)
  5 (5.00%) high mild
  9 (9.00%) high severe
mol_parsing/properties/zch4
                        time:   [115.18 ns 115.45 ns 115.75 ns]
                        change: [+2129.5% +2135.8% +2142.5%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 4 outliers among 100 measurements (4.00%)
  3 (3.00%) high mild
  1 (1.00%) high severe
mol_parsing/properties/hyd1
                        time:   [55.853 ns 55.897 ns 55.950 ns]
                        change: [+1007.2% +1017.5% +1024.2%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 10 outliers among 100 measurements (10.00%)
  5 (5.00%) high mild
  5 (5.00%) high severe
mol_parsing/properties/hyd6
                        time:   [149.29 ns 149.40 ns 149.52 ns]
                        change: [+2870.8% +2878.4% +2885.0%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 10 outliers among 100 measurements (10.00%)
  6 (6.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/scn1
                        time:   [53.224 ns 53.332 ns 53.455 ns]
Found 8 outliers among 100 measurements (8.00%)
  2 (2.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/scn2
                        time:   [60.179 ns 60.246 ns 60.316 ns]
Found 8 outliers among 100 measurements (8.00%)
  6 (6.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/sds1
                        time:   [52.981 ns 53.013 ns 53.048 ns]
Found 8 outliers among 100 measurements (8.00%)
  8 (8.00%) high severe
mol_parsing/properties/spa12
                        time:   [148.93 ns 149.10 ns 149.29 ns]
Found 6 outliers among 100 measurements (6.00%)
  3 (3.00%) high mild
  3 (3.00%) high severe
mol_parsing/properties/crs
                        time:   [80.363 ns 80.719 ns 81.395 ns]
Found 7 outliers among 100 measurements (7.00%)
  3 (3.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/sdi1
                        time:   [69.298 ns 69.372 ns 69.459 ns]
Found 11 outliers among 100 measurements (11.00%)
  5 (5.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/sdi2
                        time:   [69.271 ns 69.352 ns 69.444 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/sbv
                        time:   [56.247 ns 56.299 ns 56.354 ns]
Found 9 outliers among 100 measurements (9.00%)
  3 (3.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/sdt
                        time:   [7.7789 ns 7.7850 ns 7.7918 ns]
Found 8 outliers among 100 measurements (8.00%)
  4 (4.00%) high mild
  4 (4.00%) high severe
mol_parsing/properties/sdd
                        time:   [115.39 ns 115.49 ns 115.59 ns]
Found 5 outliers among 100 measurements (5.00%)
  5 (5.00%) high mild
mol_parsing/properties/scd
                        time:   [51.822 ns 51.860 ns 51.904 ns]
Found 4 outliers among 100 measurements (4.00%)
  4 (4.00%) high mild
mol_parsing/properties/sed
                        time:   [54.714 ns 55.145 ns 55.979 ns]
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/spl1
                        time:   [58.310 ns 58.942 ns 60.258 ns]
Found 6 outliers among 100 measurements (6.00%)
  4 (4.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/snc1
                        time:   [56.743 ns 56.830 ns 56.938 ns]
Found 8 outliers among 100 measurements (8.00%)
  2 (2.00%) high mild
  6 (6.00%) high severe
mol_parsing/properties/rgp1
                        time:   [58.193 ns 58.696 ns 59.723 ns]
Found 5 outliers among 100 measurements (5.00%)
  3 (3.00%) high mild
  2 (2.00%) high severe
mol_parsing/properties/log1
                        time:   [183.67 ns 185.20 ns 186.74 ns]

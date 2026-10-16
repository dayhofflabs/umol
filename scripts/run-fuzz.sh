#!/bin/bash

EDN_TARGETS=(fuzz_deser fuzz_dyn_edn fuzz_parse fuzz_parse_all fuzz_roundtrip)
GRAPH_IR_TARGETS=(fuzz_edits fuzz_entity_strings fuzz_molecule fuzz_reaction fuzz_reaction_span)
IO_TARGETS=(fuzz_parse_opensmiles)
MAX_TOTAL_TIME=600

cd umol-edn
for target in "${EDN_TARGETS[@]}"; do
	cargo +nightly fuzz run "$target" -- -max_total_time=$MAX_TOTAL_TIME
	if [ $? -ne 0 ]; then
		echo "Target $target failed"
		exit 1
	fi
done
cd ..

cd umol-graph-ir
for target in "${GRAPH_IR_TARGETS[@]}"; do
	mkdir -p "fuzz/corpus/$target"
	cargo +nightly fuzz run "$target" "fuzz/corpus/$target" "fuzz/seeds/$target" -- -max_total_time=$MAX_TOTAL_TIME
	if [ $? -ne 0 ]; then
		echo "Target $target failed"
		exit 1
	fi
done
cd ..

cd umol-io
for target in "${IO_TARGETS[@]}"; do
	mkdir -p "fuzz/corpus/$target"
	cargo +nightly fuzz run "$target" "fuzz/corpus/$target" "fuzz/seeds/$target" -- \
		-dict=fuzz/dictionaries/smiles.dict -max_total_time=$MAX_TOTAL_TIME
	if [ $? -ne 0 ]; then
		echo "Target $target failed"
		exit 1
	fi
done
cd ..

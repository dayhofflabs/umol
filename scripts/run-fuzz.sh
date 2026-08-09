#!/bin/bash

EDN_TARGETS=(fuzz_deser fuzz_dyn_edn fuzz_parse fuzz_parse_all fuzz_roundtrip)
GRAPH_IR_TARGETS=(fuzz_constraints fuzz_entity_strings fuzz_molecule fuzz_reaction fuzz_reaction_span)
IO_TARGETS=(fuzz_parse_opensmiles)
MAX_TOTAL_TIME=10

cd umol-edn
for target in ${EDN_TARGETS[@]}; do
	cargo fuzz run $target -- -max_total_time=$MAX_TOTAL_TIME
	if [ $? -ne 0 ]; then
		echo "Target $target failed"
		exit 1
	fi
done
cd ..

cd umol-graph-ir
for target in ${GRAPH_IR_TARGETS[@]}; do
	cargo fuzz run $target -- -max_total_time=$MAX_TOTAL_TIME
	if [ $? -ne 0 ]; then
		echo "Target $target failed"
		exit 1
	fi
done
cd ..

cd umol-io
for target in ${IO_TARGETS[@]}; do
	cargo fuzz run $target -- -max_total_time=$MAX_TOTAL_TIME
	if [ $? -ne 0 ]; then
		echo "Target $target failed"
		exit 1
	fi
done
cd ..

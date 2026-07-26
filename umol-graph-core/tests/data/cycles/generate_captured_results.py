#!/usr/bin/env python3
"""Generate captured external MCB, relevant-cycle, and URF results."""

import argparse
import ctypes
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path

import igraph as ig
import networkx as nx


JAVA_SOURCE = Path(__file__).with_name("CdkCycleFamilies.java")


def parse_graph6(source):
    data = source.encode()
    node_count = data[0] - 63
    edges = []
    bit = 0
    for second in range(1, node_count):
        for first in range(second):
            value = data[1 + bit // 6] - 63
            if value & (1 << (5 - bit % 6)):
                edges.append((first, second))
            bit += 1
    return node_count, edges


def encode_edges(edges):
    return ";".join(f"{first},{second}" for first, second in edges)


def encode_sequences(sequences, separator=";"):
    if not sequences:
        return "-"
    return separator.join(",".join(map(str, sequence)) for sequence in sequences)


def cycle_order(cycle):
    return len(cycle), tuple(cycle)


def normalize_cycles(cycles):
    return sorted((tuple(sorted(cycle)) for cycle in cycles), key=cycle_order)


def git_revision(source):
    return subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def cdk_jars(source):
    interfaces = list((source / "base/interfaces/target").glob("cdk-interfaces-*.jar"))
    core = list((source / "base/core/target").glob("cdk-core-*.jar"))
    if len(interfaces) != 1 or len(core) != 1:
        raise RuntimeError("CDK interfaces and core must be built before validation")
    return interfaces[0], core[0]


def run_cdk(source, records):
    interfaces, core = cdk_jars(source)
    dependency_path = os.pathsep.join((str(interfaces), str(core)))
    with tempfile.TemporaryDirectory() as classes:
        subprocess.run(
            ["javac", "-cp", dependency_path, "-d", classes, str(JAVA_SOURCE)],
            check=True,
        )
        result = subprocess.run(
            [
                "java",
                "-cp",
                os.pathsep.join((classes, dependency_path)),
                "CdkCycleFamilies",
            ],
            input="".join(
                f"{graph6}\t{node_count}\t{encode_edges(edges)}\n"
                for graph6, node_count, edges in records
            ),
            check=True,
            capture_output=True,
            text=True,
        )

    answers = {}
    for line in result.stdout.splitlines():
        graph6, dimension, total_length, relevant = line.split("\t")
        answers[graph6] = (
            int(dimension),
            int(total_length),
            decode_sequences(relevant),
        )
    if len(answers) != len(records):
        raise AssertionError("CDK did not return exactly one row per input graph")
    return answers


def decode_sequences(source, separator=";"):
    if source == "-":
        return []
    return [
        tuple(map(int, sequence.split(",")))
        for sequence in source.split(separator)
    ]


class RingDecomposer:
    Edge = ctypes.c_uint * 2

    class Cycle(ctypes.Structure):
        pass

    Cycle._fields_ = [
        ("edges", ctypes.POINTER(Edge)),
        ("weight", ctypes.c_uint),
        ("urf", ctypes.c_uint),
        ("rcf", ctypes.c_uint),
    ]

    def __init__(self, library):
        self.library = ctypes.CDLL(str(library))
        self.libc = ctypes.CDLL(None)
        self.libc.free.argtypes = [ctypes.c_void_p]

        cycle_pointer = ctypes.POINTER(self.Cycle)
        cycle_array = ctypes.POINTER(cycle_pointer)
        edge_pointer = ctypes.POINTER(self.Edge)
        node_pointer = ctypes.POINTER(ctypes.c_uint)

        self.library.RDL_initNewGraph.argtypes = [ctypes.c_uint]
        self.library.RDL_initNewGraph.restype = ctypes.c_void_p
        self.library.RDL_addUEdge.argtypes = [
            ctypes.c_void_p,
            ctypes.c_uint,
            ctypes.c_uint,
        ]
        self.library.RDL_addUEdge.restype = ctypes.c_uint
        self.library.RDL_calculate.argtypes = [ctypes.c_void_p]
        self.library.RDL_calculate.restype = ctypes.c_void_p
        self.library.RDL_deleteData.argtypes = [ctypes.c_void_p]
        self.library.RDL_getNofURF.argtypes = [ctypes.c_void_p]
        self.library.RDL_getNofURF.restype = ctypes.c_uint
        self.library.RDL_getWeightForURF.argtypes = [ctypes.c_void_p, ctypes.c_uint]
        self.library.RDL_getWeightForURF.restype = ctypes.c_uint
        self.library.RDL_getNofRCForURF.argtypes = [ctypes.c_void_p, ctypes.c_uint]
        self.library.RDL_getNofRCForURF.restype = ctypes.c_double
        self.library.RDL_getNofRC.argtypes = [ctypes.c_void_p]
        self.library.RDL_getNofRC.restype = ctypes.c_double
        self.library.RDL_getNodesForURF.argtypes = [
            ctypes.c_void_p,
            ctypes.c_uint,
            ctypes.POINTER(node_pointer),
        ]
        self.library.RDL_getNodesForURF.restype = ctypes.c_uint
        self.library.RDL_getEdgesForURF.argtypes = [
            ctypes.c_void_p,
            ctypes.c_uint,
            ctypes.POINTER(edge_pointer),
        ]
        self.library.RDL_getEdgesForURF.restype = ctypes.c_uint
        self.library.RDL_getRCyclesForURFIterator.argtypes = [
            ctypes.c_void_p,
            ctypes.c_uint,
        ]
        self.library.RDL_getRCyclesForURFIterator.restype = ctypes.c_void_p
        self.library.RDL_cycleIteratorAtEnd.argtypes = [ctypes.c_void_p]
        self.library.RDL_cycleIteratorAtEnd.restype = ctypes.c_int
        self.library.RDL_cycleIteratorGetCycle.argtypes = [ctypes.c_void_p]
        self.library.RDL_cycleIteratorGetCycle.restype = cycle_pointer
        self.library.RDL_cycleIteratorNext.argtypes = [ctypes.c_void_p]
        self.library.RDL_cycleIteratorNext.restype = ctypes.c_void_p
        self.library.RDL_deleteCycleIterator.argtypes = [ctypes.c_void_p]
        self.library.RDL_deleteCycle.argtypes = [cycle_pointer]
        self.library.RDL_getRCycles.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(cycle_array),
        ]
        self.library.RDL_getRCycles.restype = ctypes.c_uint
        self.library.RDL_getSSSR.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(cycle_array),
        ]
        self.library.RDL_getSSSR.restype = ctypes.c_uint
        self.library.RDL_deleteCycles.argtypes = [cycle_array, ctypes.c_uint]
        self.library.RDL_validateRingFamilies.argtypes = [
            ctypes.POINTER(ctypes.c_char_p),
            ctypes.POINTER(ctypes.c_uint),
            ctypes.c_uint,
            ctypes.c_uint,
            ctypes.c_void_p,
        ]
        self.library.RDL_validateRingFamilies.restype = ctypes.c_int

    def calculate(self, node_count, edges, validate):
        graph = self.library.RDL_initNewGraph(node_count)
        if not graph:
            raise RuntimeError("RingDecomposerLib could not allocate a graph")
        for edge_id, (first, second) in enumerate(edges):
            actual = self.library.RDL_addUEdge(graph, first, second)
            if actual != edge_id:
                raise RuntimeError("RingDecomposerLib rejected an input edge")
        data = self.library.RDL_calculate(graph)
        if not data:
            raise RuntimeError("RingDecomposerLib calculation failed")
        try:
            edge_ids = {tuple(sorted(edge)): index for index, edge in enumerate(edges)}
            edges_by_id = dict(enumerate(edges))
            dimension, total_length = self._minimum_cycle_basis(data)
            families, attributed_cycles = self._families(data, edge_ids, edges_by_id)
            relevant = normalize_cycles(cycle for family in families for cycle in family["cycles"])
            eager = self._eager_relevant_cycles(data, edge_ids)
            if eager != relevant:
                raise AssertionError("RingDecomposerLib eager and iterator results differ")
            if int(self.library.RDL_getNofRC(data)) != len(relevant):
                raise AssertionError("RingDecomposerLib global relevant-cycle count differs")

            validation = "-"
            if relevant and validate:
                self._validate_families(attributed_cycles, len(edges))
                validation = "1"
            elif relevant:
                validation = "0"

            return dimension, total_length, relevant, families, validation
        finally:
            self.library.RDL_deleteData(data)

    def _minimum_cycle_basis(self, data):
        cycles = ctypes.POINTER(ctypes.POINTER(self.Cycle))()
        dimension = self.library.RDL_getSSSR(data, ctypes.byref(cycles))
        total_length = sum(cycles[index].contents.weight for index in range(dimension))
        self.library.RDL_deleteCycles(cycles, dimension)
        return dimension, total_length

    def _families(self, data, edge_ids, edges_by_id):
        families = []
        attributed_cycles = []
        for urf in range(self.library.RDL_getNofURF(data)):
            nodes = self._family_nodes(data, urf)
            edges = self._family_edges(data, urf, edge_ids)
            cycles = self._family_cycles(data, urf, edge_ids)
            count = int(self.library.RDL_getNofRCForURF(data, urf))
            weight = self.library.RDL_getWeightForURF(data, urf)
            if count != len(cycles):
                raise AssertionError("RingDecomposerLib URF count differs from iterator")
            if any(len(cycle) != weight for cycle in cycles):
                raise AssertionError("RingDecomposerLib URF weight differs from cycles")
            derived_nodes = sorted(
                {
                    node
                    for cycle in cycles
                    for edge in cycle
                    for node in edges_by_id[edge]
                }
            )
            derived_edges = sorted({edge for cycle in cycles for edge in cycle})
            if nodes != derived_nodes or edges != derived_edges:
                raise AssertionError("RingDecomposerLib URF metadata differs from cycles")
            families.append(
                {
                    "weight": weight,
                    "count": count,
                    "nodes": nodes,
                    "edges": edges,
                    "cycles": cycles,
                }
            )
            attributed_cycles.extend((cycle, urf) for cycle in cycles)
        families.sort(
            key=lambda family: (
                family["weight"],
                family["cycles"],
                family["nodes"],
                family["edges"],
            )
        )
        return families, attributed_cycles

    def _family_nodes(self, data, urf):
        nodes = ctypes.POINTER(ctypes.c_uint)()
        count = self.library.RDL_getNodesForURF(data, urf, ctypes.byref(nodes))
        result = sorted(nodes[index] for index in range(count))
        self.libc.free(nodes)
        return result

    def _family_edges(self, data, urf, edge_ids):
        edges = ctypes.POINTER(self.Edge)()
        count = self.library.RDL_getEdgesForURF(data, urf, ctypes.byref(edges))
        result = sorted(
            edge_ids[tuple(sorted(edges[index]))]
            for index in range(count)
        )
        self.libc.free(edges)
        return result

    def _family_cycles(self, data, urf, edge_ids):
        iterator = self.library.RDL_getRCyclesForURFIterator(data, urf)
        if not iterator:
            raise RuntimeError("RingDecomposerLib could not allocate a cycle iterator")
        cycles = []
        try:
            while not self.library.RDL_cycleIteratorAtEnd(iterator):
                cycle = self.library.RDL_cycleIteratorGetCycle(iterator)
                if not cycle:
                    raise RuntimeError("RingDecomposerLib cycle iteration failed")
                try:
                    cycles.append(self._copy_cycle(cycle.contents, edge_ids))
                finally:
                    self.library.RDL_deleteCycle(cycle)
                if not self.library.RDL_cycleIteratorNext(iterator):
                    raise RuntimeError("RingDecomposerLib cycle iterator advance failed")
        finally:
            self.library.RDL_deleteCycleIterator(iterator)
        return normalize_cycles(cycles)

    def _eager_relevant_cycles(self, data, edge_ids):
        cycles = ctypes.POINTER(ctypes.POINTER(self.Cycle))()
        count = self.library.RDL_getRCycles(data, ctypes.byref(cycles))
        result = normalize_cycles(
            self._copy_cycle(cycles[index].contents, edge_ids)
            for index in range(count)
        )
        self.library.RDL_deleteCycles(cycles, count)
        return result

    @staticmethod
    def _copy_cycle(cycle, edge_ids):
        return tuple(
            sorted(
                edge_ids[tuple(sorted(cycle.edges[index]))]
                for index in range(cycle.weight)
            )
        )

    def _validate_families(self, attributed_cycles, edge_count):
        buffers = []
        urfs = []
        for cycle, urf in attributed_cycles:
            buffer = ctypes.create_string_buffer(edge_count)
            for edge in cycle:
                buffer[edge] = b"\1"
            buffers.append(buffer)
            urfs.append(urf)
        cycle_array = (ctypes.c_char_p * len(buffers))(
            *(ctypes.cast(buffer, ctypes.c_char_p) for buffer in buffers)
        )
        urf_array = (ctypes.c_uint * len(urfs))(*urfs)
        stop = ctypes.cast(self.library.RDL_no_stop_fun, ctypes.c_void_p)
        result = self.library.RDL_validateRingFamilies(
            cycle_array,
            urf_array,
            len(buffers),
            edge_count,
            stop,
        )
        if result != 0:
            raise AssertionError(
                f"RingDecomposerLib exponential validation failed with {result}"
            )

def encode_families(families):
    if not families:
        return "-"
    encoded = []
    for family in families:
        encoded.append(
            ":".join(
                (
                    str(family["weight"]),
                    str(family["count"]),
                    encode_sequences([family["nodes"]]),
                    encode_sequences([family["edges"]]),
                    encode_sequences(family["cycles"], separator="."),
                )
            )
        )
    return "/".join(encoded)


def external_minimum_cycle_bases(node_count, edges):
    networkx_graph = nx.Graph()
    networkx_graph.add_nodes_from(range(node_count))
    networkx_graph.add_edges_from(edges)
    networkx_basis = nx.minimum_cycle_basis(networkx_graph)

    igraph_graph = ig.Graph(n=node_count, edges=edges, directed=False)
    igraph_basis = igraph_graph.minimum_cycle_basis(use_cycle_order=False)

    return (
        (len(networkx_basis), sum(map(len, networkx_basis))),
        (len(igraph_basis), sum(map(len, igraph_basis))),
    )


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("corpus", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cdk-source", type=Path, required=True)
    parser.add_argument("--rdl-source", type=Path, required=True)
    parser.add_argument("--rdl-library", type=Path, required=True)
    parser.add_argument("--validate-through-order", type=int, default=5)
    args = parser.parse_args()

    records = [
        (source, *parse_graph6(source))
        for source in args.corpus.read_text().splitlines()
    ]
    cdk = run_cdk(args.cdk_source, records)
    rdl = RingDecomposer(args.rdl_library)

    graph_count = 0
    relevant_cycle_count = 0
    urf_count = 0
    validated_count = 0
    with args.output.open("w") as output:
        for source, node_count, edges in records:
            networkx_mcb, igraph_mcb = external_minimum_cycle_bases(node_count, edges)
            cdk_dimension, cdk_weight, cdk_relevant = cdk[source]
            (
                rdl_dimension,
                rdl_weight,
                rdl_relevant,
                families,
                validation,
            ) = rdl.calculate(
                node_count,
                edges,
                node_count <= args.validate_through_order,
            )

            minimum_cycle_bases = {
                "networkx": networkx_mcb,
                "igraph": igraph_mcb,
                "cdk": (cdk_dimension, cdk_weight),
                "rdl": (rdl_dimension, rdl_weight),
            }
            if len(set(minimum_cycle_bases.values())) != 1:
                raise AssertionError(f"{source}: MCB mismatch {minimum_cycle_bases}")
            if cdk_relevant != rdl_relevant:
                raise AssertionError(f"{source}: relevant-cycle mismatch")

            output.write(
                "\t".join(
                    (
                        source,
                        str(cdk_dimension),
                        str(cdk_weight),
                        encode_sequences(cdk_relevant),
                        encode_families(families),
                        validation,
                    )
                )
                + "\n"
            )
            graph_count += 1
            relevant_cycle_count += len(cdk_relevant)
            urf_count += len(families)
            validated_count += validation == "1"

    print(f"networkx={nx.__version__}")
    print(f"python_igraph={ig.__version__}")
    print(f"igraph={ig.__igraph_version__}")
    print(f"cdk_revision={git_revision(args.cdk_source)}")
    print(f"ringdecomposerlib_revision={git_revision(args.rdl_source)}")
    print(f"validate_through_order={args.validate_through_order}")
    print(f"graphs={graph_count}")
    print(f"relevant_cycles={relevant_cycle_count}")
    print(f"unique_ring_families={urf_count}")
    print(f"exponentially_validated_cyclic_graphs={validated_count}")
    print("mcb_comparison=dimension_and_total_weight")
    print("relevant_cycle_comparison=normalized_edge_sets")
    print("urf_comparison=normalized_partition_metadata_and_iterator_emission")
    print("semantic_scope=finite_simple_unweighted_undirected_graphs")
    print("normalized_failures=0")
    print(f"corpus_sha256={digest(args.corpus)}")
    print(f"results_sha256={digest(args.output)}")


if __name__ == "__main__":
    main()

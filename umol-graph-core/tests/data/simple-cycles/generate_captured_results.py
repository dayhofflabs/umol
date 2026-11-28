#!/usr/bin/env python3
"""Generate captured external simple-cycle results for graph-core tests."""

import argparse
import hashlib
from pathlib import Path

import igraph as ig
import networkx as nx


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


def normalize_node_cycle(cycle):
    cycle = tuple(cycle)
    candidates = []
    for direction in (cycle, tuple(reversed(cycle))):
        for offset in range(len(cycle)):
            candidates.append(direction[offset:] + direction[:offset])
    return min(candidates)


def node_cycles_networkx(node_count, edges):
    graph = nx.Graph()
    graph.add_nodes_from(range(node_count))
    graph.add_edges_from(edges)
    cycles = [normalize_node_cycle(cycle) for cycle in nx.simple_cycles(graph)]
    return sorted(cycles, key=lambda cycle: (len(cycle), cycle))


def node_cycles_igraph(node_count, edges):
    graph = ig.Graph(n=node_count, edges=edges, directed=False)
    cycles = [normalize_node_cycle(cycle) for cycle in graph.simple_cycles(output="vpath")]
    return sorted(cycles, key=lambda cycle: (len(cycle), cycle))


def encode_sequences(sequences):
    if not sequences:
        return "-"
    return ";".join(",".join(map(str, sequence)) for sequence in sequences)


def write_simple_results(corpus_path, output_path):
    graph_count = 0
    cycle_count = 0
    with corpus_path.open() as corpus, output_path.open("w") as output:
        for line in corpus:
            source = line.rstrip("\n")
            node_count, edges = parse_graph6(source)
            networkx_cycles = node_cycles_networkx(node_count, edges)
            igraph_cycles = node_cycles_igraph(node_count, edges)
            if networkx_cycles != igraph_cycles:
                raise AssertionError(f"external libraries disagree for {source}")
            output.write(f"{source}\t{encode_sequences(networkx_cycles)}\n")
            graph_count += 1
            cycle_count += len(networkx_cycles)
    return graph_count, cycle_count


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("corpus", type=Path)
    parser.add_argument("simple_output", type=Path)
    args = parser.parse_args()

    simple_graphs, simple_cycles = write_simple_results(
        args.corpus, args.simple_output
    )

    print(f"networkx={nx.__version__}")
    print(f"python_igraph={ig.__version__}")
    print(f"igraph={ig.__igraph_version__}")
    print(f"simple_graphs={simple_graphs}")
    print(f"simple_cycles={simple_cycles}")
    print(f"simple_sha256={digest(args.simple_output)}")


if __name__ == "__main__":
    main()

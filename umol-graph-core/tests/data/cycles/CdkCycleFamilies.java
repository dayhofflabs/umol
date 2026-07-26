import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.graph.MinimumCycleBasis;
import org.openscience.cdk.graph.RelevantCycles;

public final class CdkCycleFamilies {
    private static final Comparator<int[]> CYCLE_ORDER = (left, right) -> {
        int length = Integer.compare(left.length, right.length);
        if (length != 0) {
            return length;
        }
        for (int index = 0; index < left.length; index++) {
            int element = Integer.compare(left[index], right[index]);
            if (element != 0) {
                return element;
            }
        }
        return 0;
    };

    private CdkCycleFamilies() {
    }

    public static void main(String[] args) throws Exception {
        try (BufferedReader input = new BufferedReader(new InputStreamReader(System.in))) {
            String line;
            while ((line = input.readLine()) != null) {
                String[] fields = line.split("\\t", -1);
                String source = fields[0];
                int nodeCount = Integer.parseInt(fields[1]);
                int[][] edges = parseEdges(fields[2]);
                int[][] graph = adjacencyList(nodeCount, edges);

                MinimumCycleBasis basis = new MinimumCycleBasis(graph);
                int totalLength = Arrays.stream(basis.paths())
                        .mapToInt(path -> path.length - 1)
                        .sum();
                List<int[]> relevant = normalizeCycles(
                        new RelevantCycles(graph).paths(),
                        edges
                );

                System.out.printf(
                        "%s\t%d\t%d\t%s%n",
                        source,
                        basis.size(),
                        totalLength,
                        encodeCycles(relevant)
                );
            }
        }
    }

    private static int[][] parseEdges(String source) {
        if (source.isEmpty()) {
            return new int[0][2];
        }
        String[] encoded = source.split(";");
        int[][] edges = new int[encoded.length][2];
        for (int index = 0; index < encoded.length; index++) {
            String[] endpoints = encoded[index].split(",");
            edges[index][0] = Integer.parseInt(endpoints[0]);
            edges[index][1] = Integer.parseInt(endpoints[1]);
        }
        return edges;
    }

    private static int[][] adjacencyList(int nodeCount, int[][] edges) {
        int[] degrees = new int[nodeCount];
        for (int[] edge : edges) {
            degrees[edge[0]]++;
            degrees[edge[1]]++;
        }
        int[][] graph = new int[nodeCount][];
        for (int node = 0; node < nodeCount; node++) {
            graph[node] = new int[degrees[node]];
        }
        int[] positions = new int[nodeCount];
        for (int[] edge : edges) {
            graph[edge[0]][positions[edge[0]]++] = edge[1];
            graph[edge[1]][positions[edge[1]]++] = edge[0];
        }
        return graph;
    }

    private static List<int[]> normalizeCycles(int[][] paths, int[][] edges) {
        Map<Long, Integer> edgeIds = new HashMap<>();
        for (int edge = 0; edge < edges.length; edge++) {
            edgeIds.put(edgeKey(edges[edge][0], edges[edge][1]), edge);
        }

        List<int[]> cycles = new ArrayList<>(paths.length);
        for (int[] path : paths) {
            int[] cycle = new int[path.length - 1];
            for (int index = 0; index + 1 < path.length; index++) {
                Integer edge = edgeIds.get(edgeKey(path[index], path[index + 1]));
                if (edge == null) {
                    throw new IllegalStateException("CDK returned a non-edge");
                }
                cycle[index] = edge;
            }
            Arrays.sort(cycle);
            cycles.add(cycle);
        }
        cycles.sort(CYCLE_ORDER);
        return cycles;
    }

    private static long edgeKey(int first, int second) {
        int low = Math.min(first, second);
        int high = Math.max(first, second);
        return ((long) low << 32) | (high & 0xffffffffL);
    }

    private static String encodeCycles(List<int[]> cycles) {
        if (cycles.isEmpty()) {
            return "-";
        }
        List<String> encoded = new ArrayList<>(cycles.size());
        for (int[] cycle : cycles) {
            encoded.add(encodeIntegers(cycle));
        }
        return String.join(";", encoded);
    }

    private static String encodeIntegers(int[] values) {
        StringBuilder encoded = new StringBuilder();
        for (int index = 0; index < values.length; index++) {
            if (index > 0) {
                encoded.append(',');
            }
            encoded.append(values[index]);
        }
        return encoded.toString();
    }
}

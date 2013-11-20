#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

const int MAX_VERTICES = 1000;

typedef uint32_t Graph[MAX_VERTICES][MAX_VERTICES];

// Read graph from input stream and return number of vertices.
int readGraph(std::istream& in, Graph G) {
    // Read number of vertices.
    int N; in >> N;

    // Read vertex coordinates.
    double x[MAX_VERTICES];
    double y[MAX_VERTICES];
    for (int i = 0; i < N; ++i) {
        in >> x[i] >> y[i];
    }

    // Calculate distances.
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            G[i][j] = std::round(std::sqrt(std::pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2)));
        }
    }

    return N;
}

// Calculate and print MST using Prim's algorithm.
void primMST(Graph G, int N) {
    bool inMST[MAX_VERTICES];
    uint32_t cost[MAX_VERTICES];
    int parent[MAX_VERTICES];

    std::fill(inMST, inMST + N, false);
    std::fill(parent, parent + N, -1);
    std::fill(cost, cost + N, std::numeric_limits<uint32_t>::max());
    cost[0] = 0;

    for (int added = 0; added < N - 1; ++added) {
        // Find lowest cost vertex not in MST and add it to MST.
        int u = -1;
        uint32_t c = std::numeric_limits<int>::max();
        for (int v = 0; v < N; ++v) {
            if (cost[v] < c && !inMST[v]) {
                u = v;
                c = cost[v];
            }
        }
        inMST[u] = true;

        // For each neighbor v of u not already in MST.
        for (int v = 0; v < N; ++v) {
            if (0 < G[u][v] && G[u][v] < cost[v] && !inMST[v]) {
                // Distance from u to v is lower than current cost of v, so mark u
                // as parent of v and set cost of v to the distance from u to v.
                parent[v] = u;
                cost[v] = G[u][v];
            }
        }
    }

    // Print MST.
    for (int v = 1; v < N; ++v) {
        std::cout << v << " " << parent[v] << std::endl;
    }
}

int main(void) {

    // Read graph.
    Graph G;
    int N = readGraph(std::cin, G);

    // Calculate and print MST.
    primMST(G, N);

    return 0;
}

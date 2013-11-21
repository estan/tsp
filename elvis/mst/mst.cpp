#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <cassert>

const int MAX_VERTICES = 1000;
const int NO_VERTEX = -1;

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

// Calculate MST rooted at 0 using Prim's algorithm.
void primMST(Graph G, int parent[MAX_VERTICES], int N) {
    bool inMST[MAX_VERTICES];
    uint32_t cost[MAX_VERTICES];

    std::fill(inMST, inMST + N, false);
    std::fill(parent, parent + N, NO_VERTEX);
    std::fill(cost, cost + N, std::numeric_limits<uint32_t>::max());
    cost[0] = 0;

    for (int added = 0; added < N - 1; ++added) {
        // Find lowest cost vertex not in MST and add it to MST.
        int u = NO_VERTEX;
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
}

// Test function for the primMST function.
void testPrimMST() {
    // Run the function on a small test case.
    int N = 5;
    Graph G = {{0, 2, 0, 6, 0},
               {2, 0, 3, 8, 5},
               {0, 3, 0, 0, 7},
               {6, 8, 0, 0, 9},
               {0, 5, 7, 9, 0}};
    int parent[MAX_VERTICES];
    std::fill(parent, parent + N, 42); // Just in case.

    primMST(G, parent, N);

    // Check results.
    assert(parent[0] == NO_VERTEX);
    assert(parent[1] == 0);
    assert(parent[2] == 1);
    assert(parent[3] == 0);
    assert(parent[4] == 1);
}

int main(void) {
    // Run a small test case first.
    testPrimMST();

    // Read graph from standard input.
    Graph G;
    int N = readGraph(std::cin, G);

    // Calculate and print MST.
    int parent[MAX_VERTICES];
    primMST(G, parent, N);

    // Print MST.
    for (int v = 1; v < N; ++v) {
        std::cout << v << " " << parent[v] << std::endl;
    }

    return 0;
}

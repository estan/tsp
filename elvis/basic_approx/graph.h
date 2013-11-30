#ifndef GRAPH_H
#define GRAPH_H

#include <cstdint>

#define NO_VERTEX -1

typedef uint32_t** Graph;

/**
 * Create a new empty graph with N vertices and return it.
 */
Graph newGraph(std::size_t N) {
    Graph G = new uint32_t*[N];
    for (std::size_t i = 0; i < N; ++i) {
        G[i] = new uint32_t[N];
        std::fill(G[i], G[i] + N, 0);
    }
    return G;
}

/**
 * Delete the graph G containing N vertices.
 */
void deleteGraph(Graph G, std::size_t N) {
    if (G != 0) {
        for (std::size_t i = 0; i < N; ++i) {
            delete G[i];
        }
        delete G;
        G = 0;
    }
}

/**
 * Read graph from an input stream and return it.
 *
 * @param in Input stream.
 * @param N The number of vertices in the graph (output).
 * @return G The read graph.
 */
Graph readGraph(std::istream& in, std::size_t *N) {
    // Read number of vertices.
    in >> *N;

    // Read vertex coordinates.
    std::vector<double> x(*N);
    std::vector<double> y(*N);
    for (std::size_t i = 0; i < *N; ++i) {
        in >> x[i] >> y[i];
    }

    // Calculate distances.
    Graph G = newGraph(*N);
    for (std::size_t i = 0; i < *N; ++i) {
        for (std::size_t j = 0; j < *N; ++j) {
            G[i][j] = std::round(std::sqrt(std::pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2)));
        }
    }

    return G;
}

// Prints graph G to standard output.
void printGraph(Graph G, std::size_t N) {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            std::cout << std::left << std::setw(3) << G[i][j];
        }
        std::cout << std::endl;
    }
}

#endif // GRAPH_H

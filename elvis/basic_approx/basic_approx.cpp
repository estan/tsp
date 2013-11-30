#include <iostream>
#include <string>
#include <iomanip>
#include <stack>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "graph.h"

/**
 * Calculate minimum spanning tree using Prim's algorithm.
 *
 * @param G Input graph.
 * @param M Minimum spanning tree of G.
 * @param N Number of vertices in G.
 */
void primMST(Graph G, Graph M, std::size_t N) {
    std::vector<bool> inMST(N);
    std::vector<uint32_t> cost(N);
    std::vector<uint32_t> parent(N);

    std::fill(inMST.begin(), inMST.end(), false);
    std::fill(parent.begin(), parent.end(), NO_VERTEX);
    std::fill(cost.begin(), cost.end(), std::numeric_limits<uint32_t>::max());
    cost[0] = 0;

    for (std::size_t added = 0; added < N - 1; ++added) {
        // Find lowest cost vertex not in MST and add it to MST.
        std::size_t u = NO_VERTEX;
        uint32_t c = std::numeric_limits<int>::max();
        for (std::size_t v = 0; v < N; ++v) {
            if (cost[v] < c && !inMST[v]) {
                u = v;
                c = cost[v];
            }
        }
        inMST[u] = true;

        // For each neighbor v of u not already in MST.
        for (std::size_t v = 0; v < N; ++v) {
            if (0 < G[u][v] && G[u][v] < cost[v] && !inMST[v]) {
                // Distance from u to v is lower than current cost of v, so mark u
                // as parent of v and set cost of v to the distance from u to v.
                parent[v] = u;
                cost[v] = G[u][v];
            }
        }
    }

    // Build MST graph.
    for (std::size_t v = 1; v < N; ++v) {
        M[parent[v]][v] = G[parent[v]][v];
    }
}

/**
 * Calculates the 2-approximation of TSP from an MST.
 *
 * @param M Input MST.
 * @param N Number of vertices in M.
 * @return A 2-approximation of the TSP.
 */
std::vector<uint32_t> twoApprox(Graph M, std::size_t N) {
    std::vector<uint32_t> tour;
    std::vector<bool> visited(N);
    std::stack<int> stack;
    stack.push(0);

    while (tour.size() < N) {
        std::size_t u = stack.top();
        stack.pop();
        if (!visited[u]) {
            tour.push_back(u);
            visited[u] = true;
            for (std::size_t v = 0; v < N; ++v) {
                if (M[u][v] != 0) {
                    stack.push(v);
                }
            }
        }
    }

    return tour;
}

/**
 * Calculate nearest neighbors graph from a given graph.
 *
 * @param G Input graph.
 * @param N Number of vertices in G.
 * @return A matrix where row i contains the the neighbors of i
 *         in ascending order of proximity.
 */
Graph createNeighborsGraph(Graph G, std::size_t N) {
    Graph nG = newGraph(N);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            nG[i][j] = j;
        }
        std::sort(nG[i], nG[i] + N, [&](std::size_t j1, std::size_t j2) {
            return G[i][j1] != 0 && G[i][j1] < G[i][j2];
        });
    }
    return nG;
}

/**
 * Reverse a segment of a tour.
 *
 * This functions reverses the segment [start, end] of the given
 * tour and updates the pos vector accordingly.
 *
 * @param tour The input tour.
 * @param start Start index of segment to reverse.
 * @param end End index of segment to reverse.
 * @param pos Vector containing positions of cities in the tour, will
 *            be updated to reflect the reversal.
 */
template<typename T>
void reverse(std::vector<T> &tour, std::size_t start, std::size_t end, std::vector<std::size_t>& pos) {
    const std::size_t N = tour.size();
    bool wrapped = start > end;
    while ((!wrapped && start < end) || (wrapped && end < start)) {
        std::swap(tour[start], tour[end]);
        pos[tour[start]] = start;
        pos[tour[end]] = end;
        ++start;
        --end;
        if (start == N) {
            wrapped = !wrapped;
            start = 0;
        }
        if (end == -1) {
            wrapped = !wrapped;
            end = N - 1;
        }
    }
}

/**
 * Perform a 2-opt pass on the given tour.
 *
 * @param tour The tour to optimize.
 * @param G The input graph.
 * @param nG Nearest neighbor graph.
 * @param N Number of vertices in G.
 */
void twoOpt(std::vector<uint32_t>& tour, Graph G, Graph nG, std::size_t N) {
    // The smallest possible edge in a tour (minimum city distance in G).
    uint32_t min = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            min = std::min(min, G[i][j]);
        }
    }

    // Maximum distance in current tour.
    uint32_t max = std::numeric_limits<std::size_t>::min();
    for (std::size_t i = 0, j = 1; j < N; ++i, ++j) {
        max = std::max(max, G[tour[i]][tour[j]]);
    }
    max = std::max(max, G[tour[N - 1]][tour[0]]);

    // Tour position vector (pos[i] = j means that i is the j:th city in tour)
    std::vector<std::size_t> pos(N);
    for (std::size_t n = 0; n < N; ++n) {
        pos[tour[n]] = n;
    }

    for (std::size_t n1 = 0; n1 < N; ++n1) {
        std::size_t m1 = (n1 - 1 + N) % N;
        for (std::size_t j2 = 0; j2 < N; ++j2) {
            std::size_t c2 = nG[n1][j2]; // c2 is the j2:th closest city to n1.
            std::size_t n2 = pos[c2];    // n2 is the position of c2 in the tour.
            std::size_t m2 = (n2 - 1 + N) % N;
            if (G[tour[n1]][tour[n2]] + min > G[tour[m1]][tour[n1]] + max) {
                break; // Go to next n1.
            }
            if (m1 == m2 || n1 == n2) {
                continue; // Don't consider connecting a city with itself.
            }
            if (G[tour[m2]][tour[m1]] + G[tour[n2]][tour[n1]] <
                G[tour[m1]][tour[n1]] + G[tour[m2]][tour[n2]]) {
                /* We have the following situation:
                 *
                 * --m1 m2--      <-m1-m2--
                 *     X     ===>
                 * <-n2 n1->      --n2-n1->
                 *
                 * That is, the pair m2->m1 and n2->n1 are shorter than the current
                 * m1->n1 and m2->n2. So simply reverse the order of the sub-tour
                 * n2...m1.
                 */
                reverse(tour, n2, m1, pos);

                // Update max link.
                max = std::max(max, std::max(G[tour[m2]][tour[m1]], G[tour[n2]][tour[n1]]));

                break; // Go to next n1.
            }
        }
    }
}

/**
 * Returns the length of the given tour of G.
 */
uint64_t length(const std::vector<uint32_t>& tour, Graph G) {
    uint64_t length = 0;
    for (std::size_t i = 0; i < tour.size() - 1; ++i) {
        length += G[tour[i]][tour[i + 1]];
    }
    length += G[tour[tour.size() - 1]][tour[0]];
    return length;
}

int main(int argc, char *argv[]) {
    // Read graph from standard input.
    std::size_t N;
    Graph G = readGraph(std::cin, &N);

    // Calculate MST.
    Graph M = newGraph(N);
    primMST(G, M, N);

    // Calculate 2-approximation.
    std::vector<uint32_t> tour = twoApprox(M, N);

    // Optimize using 2-opt.
    if (argc == 2 && std::string(argv[1]).compare("-2") == 0) {
        Graph nG = createNeighborsGraph(G, N);
        uint64_t oldLength = std::numeric_limits<uint64_t>::max();
        uint64_t newLength = length(tour, G);
        do {
            twoOpt(tour, G, nG, N);
            oldLength = newLength;
            newLength = length(tour, G);
        } while (newLength < oldLength);
        deleteGraph(nG, N);
    }

    // Print tour.
    for (auto city : tour) {
        std::cout << city << std::endl;
    }

    deleteGraph(G, N);
    deleteGraph(M, N);

    return 0;
}


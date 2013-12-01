#include <iostream>
#include <string>
#include <iomanip>
#include <stack>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "matrix.h"

using namespace std;

/**
 * Calculate minimum spanning tree using Prim's algorithm.
 *
 * @param d Distance matrix.
 * @param MST Minimum spanning tree of G.
 * @param N Number of cities.
 */
void primMST(uint32_t **d, uint32_t **MST, size_t N) {
    vector<bool> inMST(N);
    vector<uint32_t> cost(N);
    vector<uint32_t> parent(N);
    fill(inMST.begin(), inMST.end(), false);
    fill(cost.begin(), cost.end(), numeric_limits<uint32_t>::max());
    cost[0] = 0;

    for (size_t added = 0; added < N - 1; ++added) {
        // Find lowest cost city not in MST and add it to MST.
        size_t u = 0;
        uint32_t c = numeric_limits<int>::max();
        for (size_t v = 0; v < N; ++v) {
            if (cost[v] < c && !inMST[v]) {
                u = v;
                c = cost[v];
            }
        }
        inMST[u] = true;

        // For each neighbor v of u not already in MST.
        for (size_t v = 0; v < N; ++v) {
            if (0 < d[u][v] && d[u][v] < cost[v] && !inMST[v]) {
                // Distance from u to v is lower than current cost of v, so mark u
                // as parent of v and set cost of v to the distance from u to v.
                parent[v] = u;
                cost[v] = d[u][v];
            }
        }
    }

    // Build MST matrix.
    for (size_t v = 1; v < N; ++v) {
        MST[parent[v]][v] = d[parent[v]][v];
    }
}

/**
 * Calculates a 2-approximation TSP tour from an MST.
 *
 * @param MST Input MST.
 * @param N Number of cities.
 * @return 2-approximation of the TSP.
 */
vector<size_t> twoApprox(uint32_t **MST, size_t N) {
    vector<size_t> tour;
    vector<bool> visited(N);
    stack<int> stack;
    stack.push(0);

    while (tour.size() < N) {
        size_t u = stack.top();
        stack.pop();
        if (!visited[u]) {
            tour.push_back(u);
            visited[u] = true;
            for (size_t v = 0; v < N; ++v) {
                if (MST[u][v] != 0) {
                    stack.push(v);
                }
            }
        }
    }

    return tour;
}

/**
 * Calculate nearest neighbors matrix from a distance matrix.
 *
 * @param d Distance matrix.
 * @param N Number of cities.
 * @return A matrix where row i contains the neighbors of city i
 *         in ascending order of proximity.
 */
uint32_t **createNeighborsMatrix(uint32_t **d, size_t N) {
    uint32_t **neighbor = newMatrix(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            neighbor[i][j] = j;
        }
        sort(neighbor[i], neighbor[i] + N, [&](size_t j1, size_t j2) {
            return d[i][j1] != 0 && d[i][j1] < d[i][j2];
        });
    }
    return neighbor;
}

/**
 * Reverse a segment of a tour.
 *
 * This functions reverses the segment [start, end] of the given
 * tour and updates the position vector accordingly.
 *
 * @param tour The input tour.
 * @param start Start index of segment to reverse.
 * @param end End index of segment to reverse.
 * @param position Vector containing positions of cities in the tour, will
 *                 be updated to reflect the reversal.
 */
template<typename T>
void reverse(vector<T> &tour, size_t start, size_t end, vector<size_t>& position) {
    const size_t N = tour.size();
    bool wrapped = start > end;
    while ((!wrapped && start < end) || (wrapped && end < start)) {
        swap(tour[start], tour[end]);
        position[tour[start]] = start;
        position[tour[end]] = end;
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
 * Reverse a segment of a tour.
 *
 * This functions reverses the segment [start, end] of the given tour.
 *
 * @param tour The input tour.
 * @param start Start index of segment to reverse.
 * @param end End index of segment to reverse.
 */
void reverse(vector<size_t> &tour, size_t start, size_t end) {
    const size_t N = tour.size();
    while (start < end) {
        swap(tour[start % N], tour[end % N]);
        ++start;
        --end;
    }
}

/**
 * Peforms a 2-opt pass on the given tour.
 *
 * This functions uses the slow and stupid O(N^2) approach.
 *
 * @param tour The tour to optimize.
 * @param d Distance matrix.
 */
void twoOptSlow(vector<size_t>& tour, uint32_t **d) {
    const size_t N = tour.size();
    size_t u, v, w, z;
    for (size_t i = 0; i < N; ++i) {
        u = tour[i];
        v = tour[(i + 1) % N];
        for (size_t j = i + 2; (j + 1) % N != i; ++j) {
            w = tour[j % N];
            z = tour[(j + 1) % N];
            if (d[u][w] + d[z][v] < d[u][v] + d[w][z]) {
                //   --u w--        --u-w->
                //      X     ===>
                //   <-z v->        <-z-v--
                reverse(tour, i + 1, j);
            }
        }
    }
}

/**
 * Perform a 2-opt pass on the given tour.
 *
 * This function uses the fast approach described on page 12 of "Large-Step
 * Markov Chains for the Traveling Salesman Problem" (Martin/Otto/Felten, 1991)
 *
 * @param tour The tour to optimize.
 * @param d Distance matrix.
 * @param neighbor Nearest neighbors matrix.
 * @param N Number of cities.
 */
void twoOpt(vector<size_t>& tour, uint32_t **d, uint32_t **neighbor, size_t N) {
    // The shortest possible link in a tour (shortest distance in G).
    uint32_t min = numeric_limits<size_t>::max();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i != j)
                min = std::min(min, d[i][j]);
        }
    }

    // Maximum distance in current tour.
    uint32_t max = numeric_limits<size_t>::min();
    for (size_t i = 0, j = 1; j < N; ++i, ++j) {
        max = std::max(max, d[tour[i]][tour[j]]);
    }
    max = std::max(max, d[tour[N - 1]][tour[0]]);

    // Tour position vector (position[i] = j means i is the j:th city in tour)
    vector<size_t> position(N);
    for (size_t n = 0; n < N; ++n) {
        position[tour[n]] = n;
    }

    for (size_t n1 = 0; n1 < N; ++n1) {
        size_t m1 = (n1 - 1 + N) % N;
        for (size_t j2 = 0; j2 < N - 1; ++j2) {
            // n2 is the tour position of the j2:th neighbor of the n1:th city in the tour.
            size_t n2 = position[neighbor[n1][j2]];
            size_t m2 = (n2 - 1 + N) % N;
            if (d[tour[n1]][tour[n2]] + min > d[tour[m1]][tour[n1]] + max) {
                break; // Go to next n1.
            }
            if (m1 == m2 || n1 == n2) {
                continue; // Don't consider connecting a city with itself.
            }
            if (d[tour[m2]][tour[m1]] + d[tour[n2]][tour[n1]] <
                d[tour[m1]][tour[n1]] + d[tour[m2]][tour[n2]]) {
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
                reverse(tour, n2, m1, position);

                // Update max link.
                max = std::max(max, std::max(d[tour[m2]][tour[m1]], d[tour[n2]][tour[n1]]));

                break; // Go to next n1.
            }
        }
    }
}

/**
 * Returns the total length of a tour.
 *
 * @param tour The input tour.
 * @param d Distance matrix.
 */
uint64_t length(const vector<size_t>& tour, uint32_t **d) {
    uint64_t length = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        length += d[tour[i]][tour[i + 1]];
    }
    length += d[tour[tour.size() - 1]][tour[0]];
    return length;
}

int main(int argc, char *argv[]) {
    // Create distance matrix from standard input.
    size_t N;
    uint32_t **d = createDistanceMatrix(cin, &N);

    // Calculate MST and 2-approximation.
    uint32_t **MST = newMatrix(N);
    primMST(d, MST, N);
    vector<size_t> tour = twoApprox(MST, N);
    deleteMatrix(MST, N);

    // Optimize using 2-opt.
    uint32_t **neighbor = createNeighborsMatrix(d, N);
    uint64_t oldLength = numeric_limits<uint64_t>::max();
    uint64_t newLength = length(tour, d);
    do {
        //twoOpt(tour, d, neighbor, N);
        twoOptSlow(tour, d);
        oldLength = newLength;
        newLength = length(tour, d);
    } while (newLength < oldLength);
    deleteMatrix(neighbor, N);

    // Print tour.
    for (auto city : tour) {
        cout << city << endl;
    }

    deleteMatrix(d, N);

    return 0;
}


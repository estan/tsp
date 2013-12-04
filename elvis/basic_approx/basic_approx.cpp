#include <iostream>
#include <string>
#include <iomanip>
#include <stack>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cassert>

#include "matrix.h"

#ifdef DEBUG
    #define LOG(x) std::cerr << x
#else
    #define LOG(x)
#endif

using namespace std;

// Largest neighborhood to explore during 2-opt.
const size_t TWOOPT_MAX_NEIGHBORS = 100;

// Output stream operator (convenient for printing tours).
ostream& operator<<(ostream& os, const vector<size_t>& v) {
    for (auto e : v)
        os << e << ' ';
    os << endl;
    return os;
}

/**
 * Create a distance matrix from an input stream and return it.
 *
 * @param in Input stream.
 * @return The read distance matrix.
 */
Matrix<uint32_t> createDistanceMatrix(istream& in) {
    // Read vertex coordinates.
    size_t N;
    in >> N;
    vector<double> x(N);
    vector<double> y(N);
    for (size_t i = 0; i < N; ++i) {
        in >> x[i] >> y[i];
    }

    // Calculate distance matrix.
    Matrix<uint32_t> d(N, N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            d[i][j] = d[j][i] = round(sqrt(pow(x[i]-x[j], 2) + pow(y[i]-y[j], 2)));
        }
    }

    return d;
}

/**
 * Calculate minimum spanning tree using Prim's algorithm.
 *
 * @param d Distance matrix.
 * @return The minimum spanning tree.
 */
Matrix<uint32_t> primMST(const Matrix<uint32_t>& d) {
    size_t N = d.rows();

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
                // d[u][v] is lower than current cost of v, so mark u
                // as parent of v and set cost of v to d[u][v].
                parent[v] = u;
                cost[v] = d[u][v];
            }
        }
    }

    // Build MST matrix.
    Matrix<uint32_t> MST(N, N);
    for (size_t v = 1; v < N; ++v) {
        MST[parent[v]][v] = d[parent[v]][v];
    }

    return MST;
}

/**
 * Calculates a 2-approximation TSP tour from a minimum spanning tree.
 *
 * @param MST Input minimum spanning tree.
 * @return 2-approximation TSP tour.
 */
vector<uint16_t> twoApprox(const Matrix<uint32_t>& MST) {
    size_t N = MST.rows();

    vector<uint16_t> tour;
    vector<bool> visited(N);
    stack<uint16_t> stack;
    stack.push(0);

    while (tour.size() < N) {
        size_t u = stack.top();
        stack.pop();
        if (!visited[u]) {
            tour.push_back(u);
            visited[u] = true;
            for (uint16_t v = 0; v < N; ++v) {
                if (MST[u][v] != 0) {
                    stack.push(v);
                }
            }
        }
    }

    return tour;
}

/**
 * Calculates a greedy TSP tour.
 *
 * This is the naive algorithm given in the Kattis problem description.
 *
 * @param d Distance matrix.
 * @return Greedy TSP tour.
 */
vector<uint16_t> greedy(const Matrix<uint32_t>& d) {
    size_t N = d.rows();
    vector<uint16_t> tour(N);
    vector<bool> used(N, false);
    used[0] = true;
    for (size_t i = 1; i < N; ++i) {
        // Find k, the closest city to the (i - 1):th city in tour.
        int32_t k = -1;
        for (uint16_t j = 0; j < N; ++j) {
            if (!used[j] && (k == -1 || d[tour[i-1]][j] < d[tour[i-1]][k])) {
                k = j;
            }
        }
        tour[i] = k;
        used[k] = true;
    }
    return tour;
}

/**
 * Calculate K-nearest neighbors matrix from a distance matrix.
 *
 * @param d Distance matrix.
 * @return d.rows() x (d.cols - 1) matrix where element i,j is
 *         the j:th nearest neighbor of city i.
 */
Matrix<uint16_t> createNeighborsMatrix(const Matrix<uint32_t>& d, size_t K) {
    size_t N = d.rows();
    size_t M = d.cols() - 1;
    K = min(M, K);
    Matrix<uint16_t> neighbor(N, K);
    std::vector<uint16_t> row(M); // For sorting.

    for (size_t i = 0; i < N; ++i) {
        // Fill row with 0, 1, ..., i - 1, i + 1, ..., M - 1.
        uint16_t k = 0;
        for (size_t j = 0; j < M; ++j, ++k) {
            row[j] = (i == j) ? ++k : k;
        }
        // Sort K nearest row elements by distance to i.
        partial_sort(row.begin(), row.begin() + K, row.end(),
            [&](uint16_t j, uint16_t k) {
                return d[i][j] < d[i][k];
            }
        );
        // Copy first K elements (now sorted) to neighbor matrix.
        copy(row.begin(), row.begin() + K, neighbor[i]);
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
void reverse(vector<uint16_t> &tour, size_t start, size_t end,
        vector<uint16_t>& position) {
    size_t N = tour.size();
    size_t numSwaps = (((start <= end ? end - start : (end + N) - start) + 1)/2);
    uint16_t i = start;
    uint16_t j = end;
    for (size_t n = 0; n < numSwaps; ++n) {
        swap(tour[i], tour[j]);
        position[tour[i]] = i;
        position[tour[j]] = j;
        i = (i + 1) % N;
        j = ((j + N) - 1) % N;
    }
}

/**
 * Returns the shortest distance d[i][j], i != j in the given distance matrix.
 *
 * @param d Distance matrix.
 * @return Minimum distance in d.
 */
uint32_t minDistance(const Matrix<uint32_t>& d) {
    size_t N = d.rows();
    uint32_t min = numeric_limits<size_t>::max();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i != j)
                min = std::min(min, d[i][j]);
        }
    }
    return min;
}

/**
 * Returns the longest inter-city distance in the given tour.
 *
 * @param tour Input tour.
 * @param d Distance matrix.
 * @return Maximum inter-city distance in tour.
 */
uint32_t maxDistance(const vector<uint16_t>& tour, const Matrix<uint32_t>& d) {
    size_t N = tour.size();
    uint32_t max = 0;
    for (size_t i = 0, j = 1; i < N; ++i, ++j) {
        max = std::max(max, d[tour[i]][tour[j % N]]);
    }
    return max;
}

/**
 * Returns the total length of a tour.
 *
 * @param tour The input tour.
 * @param d Distance matrix.
 * @return The total length of the tour.
 */
uint64_t length(const vector<uint32_t>& tour, const Matrix<uint32_t>& d) {
    size_t N = tour.size();
    uint64_t length = 0;
    for (size_t i = 0, j = 1; i < N; ++i, ++j) {
        length += d[tour[i]][tour[j % N]];
    }
    return length;
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
 * @param min Shortest possible inter-city distance.
 * @return true if the tour was improved, otherwise false.
 */
bool twoOpt(vector<uint16_t>& tour, const Matrix<uint32_t>& d,
        const Matrix<uint16_t>& neighbor, uint32_t min) {
    size_t N = d.rows();
    bool didImprove = false;

    // Initialize city tour position vector.
    vector<uint16_t> position(N);
    for (uint16_t i = 0; i < N; ++i) {
        position[tour[i]] = i; // tour[i] is the i:th city in tour.
    }

    uint16_t u, v, w, z;
    size_t u_i, v_i, w_i, z_i;
    uint32_t max = maxDistance(tour, d); // Longest link in current tour.

    for (u_i = 0, v_i = 1; u_i < N; ++u_i, ++v_i) {
        // For each edge (u, v).
        u = tour[u_i];
        v = tour[v_i % N];
        for (size_t n = 0; n < neighbor.cols(); ++n) {
            // Visit nearby edges (w, z).
            w_i = position[neighbor[u][n]];
            z_i = w_i + 1;
            w = tour[w_i]; // w is the n:th closest neighbor of u.
            z = tour[z_i % N];

            if (v == w || z == u) {
                continue; // Skip adjacent edges.
            }

            // d[u][w] + min is a lower bound on new length.
            // d[u][v] + max is an upper bound on old length.
            if (d[u][w] + min > d[u][v] + max) {
                break; // Go to next edge (u, v).
            }

            if (d[u][w] + d[v][z] < d[u][v] + d[w][z]) {
                //   --u w--        --u-w->
                //      X     ===>
                //   <-z v->        <-z-v--
                reverse(tour, v_i % N, w_i, position);
                
                // FIXME: Is this enough to update max?
                max = std::max(max, std::max(d[u][w], d[v][z]));

                didImprove = true;
                break;
            }
        }
    }

    return didImprove;
}

int main(int argc, char *argv[]) {
    // Create distance matrix from standard input.
    Matrix<uint32_t> d = createDistanceMatrix(cin);

    // Calculate 2-approximation tour from minimum spanning tree.
    //Matrix<uint32_t> MST = primMST(d);
    //vector<uint16_t> tour = twoApprox(MST);
    
    // Calculate a greedy tour.
    vector<uint16_t> tour = greedy(d);

    // Calculate nearest neighbors and smallest possible city distance.
    Matrix<uint16_t> neighbor = createNeighborsMatrix(d, TWOOPT_MAX_NEIGHBORS);
    uint32_t min = minDistance(d);

    // Optimize tour using 2-opt.
    bool didImprove;
    do {
        didImprove = twoOpt(tour, d, neighbor, min);
    } while (didImprove);

    // Print tour.
    for (auto city : tour) {
        cout << city << endl;
    }

    return 0;
}


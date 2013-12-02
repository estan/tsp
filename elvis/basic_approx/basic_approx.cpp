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
    Matrix<uint32_t> d(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            d[i][j] = round(sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2)));
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
    size_t N = d.size();

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
    Matrix<uint32_t> MST(N);
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
vector<size_t> twoApprox(const Matrix<uint32_t>& MST) {
    size_t N = MST.size();

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
 * @return A matrix where row i contains the neighbors of city i
 *         in ascending order of proximity.
 */
Matrix<uint32_t> createNeighborsMatrix(const Matrix<uint32_t>& d) {
    size_t N = d.size();
    Matrix<uint32_t> neighbor(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            neighbor[i][j] = j;
        }
        sort(neighbor[i], neighbor[i] + N, [&](size_t j1, size_t j2) {
            return i != j1 && (i == j2 || d[i][j1] < d[i][j2]);
        });
    }
    return neighbor;
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
void reverse(vector<size_t> &tour, size_t start, size_t end, vector<size_t>& position) {
    size_t N = tour.size();
    size_t numSwaps = (((start <= end ? end - start : (end + N) - start) + 1)/2);
    int i = start;
    int j = end;
    for (size_t n = 0; n < numSwaps; ++n) {
        swap(tour[i], tour[j]);
        position[tour[i]] = i;
        position[tour[j]] = j;
        i = (i + 1) % N;
        j = (j - 1 + N) % N;
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
#if 0
void twoOptSlow(vector<size_t>& tour, const Matrix<uint32_t>& d) {
    size_t N = tour.size();
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
#endif

/**
 * Returns the minimum distance d[i][j], i != j in the given distance matrix.
 *
 * @param d Distance matrix.
 * @return Minimum distance in d.
 */
uint32_t minDistance(const Matrix<uint32_t>& d) {
    size_t N = d.size();
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
 * Returns the maximum inter-city distance in the given tour.
 *
 * @param tour Input tour.
 * @param d Distance matrix.
 * @return Maximum inter-city distance in tour.
 */
uint32_t maxDistance(const vector<size_t>& tour, const Matrix<uint32_t>& d) {
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
 */
uint64_t length(const vector<size_t>& tour, const Matrix<uint32_t>& d) {
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
 */
bool twoOpt(vector<size_t>& tour, const Matrix<uint32_t>& d, const Matrix<uint32_t>& neighbor, uint32_t min) {
    size_t N = d.size();
    bool isTwoOpt = true;

    // Initialize city tour position vector.
    vector<size_t> position(N);
    for (size_t i = 0; i < N; ++i) {
        position[tour[i]] = i; // tour[i] is the i:th city in tour.
    }

    size_t u, v, w, z;                   // We consider exchanging (u, v) and (w, z).
    size_t u_i, v_i, w_i, z_i;           // City positions in tour.
    uint32_t max = maxDistance(tour, d); // Longest link in current tour.

    for (u_i = 0, v_i = 1; u_i < N; ++u_i, ++v_i) {
        u = tour[u_i];
        v = tour[v_i % N];
        for (size_t n = 0; n < N - 1; ++n) {
            w_i = position[neighbor[u][n]];
            z_i = w_i + 1;
            w = tour[w_i]; // w is the n:th closest neighbor of u.
            z = tour[z_i % N];

            if (v == w || z == u) {
                continue;
            }

            // d[u][w] + min is a lower bound on new length.
            // d[u][v] + max is an upper bound on old length.
            if (d[u][w] + min > d[u][v] + max) {
                break; // Go to next edge (u, v).
            }

            if (d[u][w] + d[v][z] < d[u][v] + d[w][z]) { // Or < ?
                //   --u w--        --u-w->
                //      X     ===>
                //   <-z v->        <-z-v--
                reverse(tour, v_i % N, w_i, position);
                
                // FIXME: Update max incrementally.
                //max = std::max(max, std::max(d[u][w], d[v][z]));
                max = maxDistance(tour, d);

                isTwoOpt = false;
                break;
            }
        }
    }

    return isTwoOpt;
}

int main(int argc, char *argv[]) {
    // Create distance matrix from standard input.
    Matrix<uint32_t> d = createDistanceMatrix(cin);

    // Calculate MST and 2-approximation.
    Matrix<uint32_t> MST = primMST(d);
    vector<size_t> tour = twoApprox(MST);

    // Optimize using 2-opt.
    Matrix<uint32_t> neighbor = createNeighborsMatrix(d);
    uint32_t min = minDistance(d);
    bool isTwoOpt = false;
    do {
        isTwoOpt = twoOpt(tour, d, neighbor, min);
    } while (!isTwoOpt);

    // Print tour.
    for (auto city : tour) {
        cout << city << endl;
    }

    return 0;
}


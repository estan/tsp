#include <iostream>
#include <string>
#include <iomanip>
#include <stack>
#include <limits>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>
#include <cstdint>
#include <cassert>

#include "matrix.h"

#ifdef DEBUG
    #define LOG(x) std::cerr << x << std::endl
#else
    #define LOG(x)
#endif

using namespace std;

// Random number generator.
random_device rd;
default_random_engine rng(rd());

// Size of nearest neighbors matrix.
const static size_t MAX_K = 20;

/*
 * Some small helper functions.
 */

// Return the current time.
static inline std::chrono::time_point<std::chrono::high_resolution_clock> now() {
    return std::chrono::high_resolution_clock::now();
}

// Return a duration as number of milliseconds.
template<typename T, typename E>
static inline T milliseconds(const chrono::duration<T, E>& d) {
    return chrono::duration_cast<std::chrono::milliseconds>(d).count();
}

// 3-argument maximum function.
template<typename T>
static inline T maximum(const T& a, const T& b, const T& c) {
    return std::max(a, std::max(b, c));
}

// 4-argument maximum function.
template<typename T>
static inline T maximum(const T& a, const T& b, const T& c, const T& d) {
    return std::max(a, std::max(b, std::max(c, d)));
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
 * Calculates a greedy TSP tour starting at city u.
 *
 * This is the naive algorithm given in the Kattis problem description.
 *
 * @param d Distance matrix.
 * @param u Start city of tour.
 * @return Greedy TSP tour.
 */
inline vector<uint16_t> greedy(const Matrix<uint32_t>& d, size_t u) {
    size_t N = d.rows();
    vector<uint16_t> tour(N);
    vector<bool> used(N, false);
    tour[0] = u;
    used[u] = true;
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
 * @return d.rows() x K matrix where element i,j is the j:th nearest
 *         neighbor of city i.
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
 * @param position Position of each city in the input tour. Will be updated.
 */
inline void reverse(vector<uint16_t> &tour, size_t start, size_t end,
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
 * Returns the total length of a tour.
 *
 * @param tour The input tour.
 * @param d Distance matrix.
 * @return The total length of the tour.
 */
inline uint64_t length(const vector<uint16_t>& tour, const Matrix<uint32_t>& d) {
    size_t N = tour.size();
    uint64_t length = 0;
    for (size_t i = 0, j = 1; i < N; ++i, ++j) {
        length += d[tour[i]][tour[j % N]];
    }
    return length;
}

/**
 * Order three edges by tour position.
 *
 * This function sets edges AB, CD, EF to GH, IJ, KL ordered by their tour
 * positions modulo the tour length (given in G_i, H_i et.c). This is a
 * helper function for the inner loop of threeOpt(...). GH, IJ and KL are
 * assumed to be disjoint.
 *
 * E.g if GH, IJ and KL have the order ..->GH->..->IJ->..->KL->.., then
 * AB = GH, CD = IJ, EF = KL, else AB = IJ, CD = GH, EF = KL.
 */
inline void ordered(
        uint16_t& A, size_t& A_i, uint16_t& B, size_t& B_i,
        uint16_t& C, size_t& C_i, uint16_t& D, size_t& D_i,
        uint16_t& E, size_t& E_i, uint16_t& F, size_t& F_i,

        uint16_t G, size_t G_i, uint16_t H, size_t H_i,
        uint16_t I, size_t I_i, uint16_t J, size_t J_i,
        uint16_t K, size_t K_i, uint16_t L, size_t L_i) {
    E = K; E_i = K_i;
    F = L; F_i = L_i;
    if ((I_i < G_i && G_i < K_i) ||
        (K_i < I_i && I_i < G_i) ||
        (G_i < K_i && K_i < I_i)) {
        A = I; A_i = I_i;
        B = J; B_i = J_i;
        C = G; C_i = G_i;
        D = H; D_i = H_i;
    } else {
        A = G; A_i = G_i;
        B = H; B_i = H_i;
        C = I; C_i = I_i;
        D = J; D_i = J_i;
    }
}

/**
 * Optimizes the given tour using 2-opt.
 *
 * This function uses the fast approach described on page 12-13 of "Large-Step
 * Markov Chains for the Traveling Salesman Problem" (Martin/Otto/Felten, 1991)
 *
 * @param tour The tour to optimize.
 * @param d Distance matrix.
 * @param neighbor Nearest neighbors matrix.
 * @param position Position of each city in the input tour. Will be updated.
 * @param max Longest inter-city distance in input tour. Will be updated.
 * @param min Shortest possible inter-city distance.
 */
inline void twoOpt(vector<uint16_t>& tour, const Matrix<uint32_t>& d,
        const Matrix<uint16_t>& neighbor, vector<uint16_t> &position,
        uint32_t& max, uint32_t min) {
    size_t N = d.rows(); // Number of cities.

    // Candidate edges uv, wz and their positions in tour.
    uint16_t u, v, w, z;
    size_t u_i, v_i, w_i, z_i;

    bool locallyOptimal = false;
    while (!locallyOptimal) {
        locallyOptimal = true;

        // For each edge uv.
        for (u_i = 0, v_i = 1; u_i < N; ++u_i, ++v_i) {
            u = tour[u_i];
            v = tour[v_i % N];

            // For each edge wz (w k:th closest neighbor of u).
            for (size_t k = 0; k < neighbor.cols(); ++k) {
                // Visit nearby edges (w, z).
                w_i = position[neighbor[u][k]];
                z_i = w_i + 1;
                w = tour[w_i];
                z = tour[z_i % N];

                if (v == w || z == u) {
                    continue; // Skip adjacent edges.
                }

                // d[u][w] + min is a lower bound on new length.
                // d[u][v] + max is an upper bound on old length.
                if (d[u][w] + min > d[u][v] + max) {
                    break; // Go to next edge uv.
                }

                if (d[u][w] + d[v][z] < d[u][v] + d[w][z]) {
                    //   --u w--        --u-w->
                    //      X     ===>
                    //   <-z v->        <-z-v--
                    reverse(tour, v_i % N, w_i, position);
                    max = maximum(max, d[u][w], d[v][z]);
                    locallyOptimal = false;
                    break;
                }
            }
        }
    }
}

/**
 * Optimizes the given tour using 3-opt.
 *
 * This function uses the fast approach described on page 12-15 of "Large-Step
 * Markov Chains for the Traveling Salesman Problem" (Martin/Otto/Felten, 1991)
 *
 * The algorithm will only consider "real" 3-exchanges involving all three
 * edges. So for best results, the tour should be preprocessed with the 2-opt
 * algorithm first.
 *
 * @param tour The tour to optimize.
 * @param d Distance matrix.
 * @param neighbor Nearest neighbors matrix.
 * @param position Position of each city in the input tour. Will be updated.
 * @param max Longest inter-city distance in input tour. Will be updated.
 * @param min Shortest possible inter-city distance.
 */
inline void threeOpt(vector<uint16_t>& tour, const Matrix<uint32_t>& d,
        const Matrix<uint16_t>& neighbor, vector<uint16_t>& position,
        uint32_t& max, uint32_t min,
        const chrono::time_point<chrono::high_resolution_clock>& deadline) {
    const size_t N = d.rows();   // Number of cities.

    // Candidate edges PQ, RS, TU and their positions in tour.
    uint16_t P, Q, R, S, T, U;
    size_t P_i, Q_i, R_i, S_i, T_i, U_i;

    // AB, CD, EF is PQ, RS, TU, but in tour order.
    uint16_t A, B, C, D, E, F;
    size_t A_i, B_i, C_i, D_i, E_i, F_i;

    bool locallyOptimal = false;
    while (!locallyOptimal) {
        locallyOptimal = true;

        // For each edge PQ.
        for (size_t i = 0; i < N; ++i) {
            P_i = i;
            Q_i = (P_i + 1) % N;
            P = tour[P_i];
            Q = tour[Q_i];

            if (chrono::high_resolution_clock::now() > deadline)
                return; // Deadline has passed, return early.

            // For each edge RS (S j:th nearest neighbor of P).
            for (size_t j = 0; j < neighbor.cols(); ++j) {
                S_i = position[neighbor[P][j]];
                R_i = (S_i + N - 1) % N;
                R = tour[R_i];
                S = tour[S_i];

                if (P == R || R == Q) // RS same as, or follows, PQ.
                    continue; // Go to next edge RS.

                if (d[P][S] + 2 * min > d[P][Q] + 2 * max)
                    break;

                if (d[P][S] + 2 * min > d[P][Q] + d[R][S] + max)
                    continue;

                // For each edge TU (U k:th nearest neighbor of P).
                for (size_t k = 0; k < neighbor.cols(); ++k) {
                    U_i = position[neighbor[P][k]];
                    T_i = (U_i + N - 1) % N;
                    T = tour[T_i];
                    U = tour[U_i];

                    if (U == S || // TU same as RS.
                        T == S || // TU follows RS.
                        U == R || // TU preceeds RS.
                        T == P || // TU same as PQ.
                        T == Q)   // TU follows PQ.
                        continue; // Go to next edge TU.

                    if (d[P][S] + d[Q][U] + min > d[P][Q] + d[R][S] + max)
                        break;

                    // Let AB, CD, EF be the edges PQ, RS, TU in tour order.
                    ordered(A, A_i, B, B_i, C, C_i, D, D_i, E, E_i, F, F_i,
                            P, P_i, Q, Q_i, R, R_i, S, S_i, T, T_i, U, U_i);

                    // Try exchanging AB, CD and EF.
                    bool changed = false;
                    uint32_t d_AB_CD_EF = d[A][B] + d[C][D] + d[E][F];
                    if (d[D][A] + d[F][B] + d[C][E] < d_AB_CD_EF) {
                        // Change AB, CD, EF to DA, FB, CE.
                        reverse(tour, F_i, A_i, position);
                        reverse(tour, D_i, E_i, position);
                        max = maximum(max, d[D][A], d[F][B], d[C][E]);
                        changed = true;
                    } else if (d[B][D] + d[E][A] + d[F][C] < d_AB_CD_EF) {
                        // Change AB, CD, EF to BD, EA, FC.
                        reverse(tour, F_i, A_i, position);
                        reverse(tour, B_i, C_i, position);
                        max = maximum(max, d[B][D], d[E][A], d[F][C]);
                        changed = true;
                    } else if (d[A][C] + d[B][E] + d[D][F] < d_AB_CD_EF) {
                        // Change AB, CD, EF to AC, BE, DF.
                        reverse(tour, B_i, C_i, position);
                        reverse(tour, D_i, E_i, position);
                        max = maximum(max, d[A][C], d[B][E], d[D][F]);
                        changed = true;
                    } else if (d[B][E] + d[D][A] + d[F][C] < d_AB_CD_EF) {
                        // Change AB, CD, EF to BE, DA, FC.
                        reverse(tour, A_i, F_i, position);
                        reverse(tour, B_i, C_i, position);
                        reverse(tour, D_i, E_i, position);
                        max = maximum(max, d[B][E], d[D][A], d[F][C]);
                        changed = true;
                    }

                    if (changed) {
                        locallyOptimal = false;
                        goto next_AB; // Go to next edge AB.
                    }
                }
            }
            next_AB: continue;
        }
    }
}

/**
 * Perform a random 4-opt ("double bridge") move on a tour.
 *
 * E.g.
 *
 *    A--B             A  B
 *   /    \           /|  |\
 *  H      C         H------C
 *  |      |   -->     |  |
 *  G      D         G------D
 *   \    /           \|  |/
 *    F--E             F  E
 *
 * Where edges AB, CD, EF and GH are chosen randomly.
 *
 * @param tour Input tour (must have at least 8 cities).
 * @param d Distance matrix.
 * @return The new tour.
 */
inline vector<uint16_t> doubleBridge(
        const vector<uint16_t>& tour,
        const Matrix<uint32_t>& d) {
    const size_t N = tour.size();
    vector<uint16_t> newTour;
    newTour.reserve(N);
    uniform_int_distribution<size_t> dist(1, N / 4);
    size_t A = dist(rng);
    size_t B = A + dist(rng);
    size_t C = B + dist(rng);
    copy(tour.begin(), tour.begin() + A, back_inserter(newTour));
    copy(tour.begin() + C, tour.end(), back_inserter(newTour));
    copy(tour.begin() + B, tour.begin() + C, back_inserter(newTour));
    copy(tour.begin() + A, tour.begin() + B, back_inserter(newTour));
    return newTour;
}

/**
 * Approximates optimal TSP tour through graph read from the given input stream.
 *
 * The function will try to return before the given deadline, but expect some
 * variations.
 *
 * @param in Input stream.
 * @param deadline Time available in milliseconds.
 * @return An approximation of the optimal TSP tour.
 */
std::vector<uint16_t> approximate(istream &in,
        const chrono::time_point<chrono::high_resolution_clock>& deadline) {
    /*
     * Initialization.
     */

    // Calculate distance / K-nearest neighbors matrix.
    const Matrix<uint32_t> d = createDistanceMatrix(in);
    const Matrix<uint16_t> neighbor = createNeighborsMatrix(d, MAX_K);
    const uint32_t min = minDistance(d); // Shortest distance.
    const size_t N = d.rows();           // Number of cities.

    // Generate initial greedy tour.
    vector<uint16_t> tour = greedy(d, 0);

    // Create max / position for initial 2-opt.
    vector<uint16_t> position(N);
    uint32_t max = 0;
    for (uint16_t i = 0; i < N; ++i) {
        max = std::max(max, d[i][(i + 1) % N]); // Maximum distance in tour.
        position[tour[i]] = i;                  // tour[i] is j:th city in tour.
    }

    // Optimize tour with 2-opt.
    twoOpt(tour, d, neighbor, position, max, min);

    // Some main loop statistics.
    size_t i = 0;                          // Number of iterations of main loop.
    chrono::duration<long, milli> totT(0); // Total time spent in main loop.
    chrono::duration<long, milli> avgT(0); // Average main loop iteration time.

    // Deadline for 3-opt inside main loop is 50 ms before hard deadline.
    chrono::duration<long, milli> fifty_ms(50);
    auto threeOptDeadline = deadline - fifty_ms;

    vector<uint16_t> shortestTour = tour;          // Best tour found.
    uint64_t shortestTourLength = length(tour, d); // Length of best tour found.

    /*
     * Main loop.
     *
     * We repeatedly
     *
     *   1) "Kick" the tour with a random 4-exchange.
     *   2) Optimize the tour with 2-opt + 3-opt.
     *
     * until only max(50, 2 * average iteration time) milliseconds remains
     * before deadline, and then pick the shortest tour we found.
     */
    for (i = 0; (now() + std::max(fifty_ms, 2 * avgT)) < deadline; ++i) {
        auto start = now();

        if (N >= 8) {
            // Perform random 4-opt "double bridge" move.
            tour = doubleBridge(tour, d);
        } else {
            // Tiny tour, so just shuffle it instead.
            shuffle(tour.begin(), tour.end(), rng);
        }

        // Update max / position needed by fast 2/3-opt.
        max = 0;
        for (uint16_t j = 0; j < N; ++j) {
            max = std::max(max, d[tour[j]][tour[(j + 1) % N]]);
            position[tour[j]] = j;
        }

        // Optimize tour with 2-opt + 3-opt.
        twoOpt(tour, d, neighbor, position, max, min);
        threeOpt(tour, d, neighbor, position, max, min, threeOptDeadline);

        uint64_t tourLength = length(tour, d);
        if (tourLength < shortestTourLength) {
            // Shorter tour found.
            shortestTour = tour;
            shortestTourLength = tourLength;
        }

        // Collect statistics.
        auto end = now();
        totT += chrono::duration_cast<std::chrono::milliseconds>(end - start);
        avgT = totT / (i + 1);
    }

    LOG("i: " << i);
    LOG("totT: " << milliseconds(totT) << " ms");
    LOG("avgT: " << milliseconds(avgT) << " ms");

    return shortestTour;
}

int main(int argc, char *argv[]) {

    // Establish a deadline 1950 ms into the future.
    auto deadline = chrono::high_resolution_clock::now();
    deadline += chrono::duration<int, milli>(1950);

    // Approximate a TSP tour.
    vector<uint16_t> tour = approximate(cin, deadline);

    // Print the tour.
    for (auto city : tour) {
        cout << city << endl;
    }

    return 0;
}


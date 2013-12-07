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
    #define LOG(x) std::cerr << x
#else
    #define LOG(x)
#endif

using namespace std;

// Random number generator.
random_device rd;
default_random_engine rng(rd());

// Largest neighborhood to explore during 2-opt/3-opt optimization.
const size_t OPTIMIZATION_MAX_NEIGHBORS = 100;

// Output stream operator (convenient for printing tours).
ostream& operator<<(ostream& os, const vector<uint16_t>& v) {
    for (auto e : v)
        os << e << ' ';
    os << endl;
    return os;
}

// 3-argument maximum function.
template<typename T>
inline T maximum(const T& a, const T& b, const T& c) {
    return std::max(a, std::max(b, c));
}

// 4-argument maximum function.
template<typename T>
inline T maximum(const T& a, const T& b, const T& c, const T& d) {
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
 * Calculates a greedy TSP tour starting at city u.
 *
 * This is the naive algorithm given in the Kattis problem description.
 *
 * @param d Distance matrix.
 * @param u Start city of tour.
 * @param length Tour length (output).
 * @return Greedy TSP tour.
 */
inline vector<uint16_t> greedy(const Matrix<uint32_t>& d, size_t u, uint64_t& length) {
    size_t N = d.rows();
    vector<uint16_t> tour(N);
    vector<bool> used(N, false);
    tour[0] = u;
    used[u] = true;
    length = 0;
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
        length += d[tour[i - 1]][k];
    }
    length += d[tour[N - 1]][tour[0]];
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
uint64_t length(const vector<uint16_t>& tour, const Matrix<uint32_t>& d) {
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
 * @param length Length of the input tour. Will be updated.
 */
inline void twoOpt(vector<uint16_t>& tour, const Matrix<uint32_t>& d,
        const Matrix<uint16_t>& neighbor, vector<uint16_t> &position,
        uint32_t& max, uint32_t min, uint64_t& length) {
    size_t N = d.rows();         // Number of cities.
    uint64_t oldLength = length; // Length of old tour.

    // Candidate edges uv, wz and their positions in tour.
    uint16_t u, v, w, z;
    size_t u_i, v_i, w_i, z_i;

    do {
        oldLength = length;

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

                int32_t delta = (d[u][w] + d[v][z]) - (d[u][v] + d[w][z]);
                if (delta < 0) {
                    //   --u w--        --u-w->
                    //      X     ===>
                    //   <-z v->        <-z-v--
                    reverse(tour, v_i % N, w_i, position);
                    max = maximum(max, d[u][w], d[v][z]);
                    length += delta;
                    break;
                }
            }
        }
    } while (length < oldLength);
}

/**
 * Optimizes the given tour using 3-opt.
 *
 * This function uses the fast approach described on page 12-15 of "Large-Step
 * Markov Chains for the Traveling Salesman Problem" (Martin/Otto/Felten, 1991)
 *
 * @param tour The tour to optimize.
 * @param d Distance matrix.
 * @param neighbor Nearest neighbors matrix.
 * @param position Position of each city in the input tour. Will be updated.
 * @param max Longest inter-city distance in input tour. Will be updated.
 * @param min Shortest possible inter-city distance.
 * @param length Length of the input tour. Will be updated.
 */
inline void threeOpt(vector<uint16_t>& tour, const Matrix<uint32_t>& d,
        const Matrix<uint16_t>& neighbor, vector<uint16_t>& position,
        uint32_t& max, uint32_t min, uint64_t& length,
        const chrono::time_point<chrono::high_resolution_clock>& deadline) {
    const size_t N = d.rows();   // Number of cities.
    uint64_t oldLength = length; // Length of old tour.

    // Candidate edges PQ, RS, TU and their positions in tour.
    uint16_t P, Q, R, S, T, U;
    size_t P_i, Q_i, R_i, S_i, T_i, U_i;

    // AB, CD, EF is PQ, RS, TU, but in tour order.
    uint16_t A, B, C, D, E, F;
    size_t A_i, B_i, C_i, D_i, E_i, F_i;

    do {
        oldLength = length;

        // For each edge PQ.
        for (size_t i = 0; i < N; ++i) {
            P_i = i;
            Q_i = (P_i + 1) % N;
            P = tour[P_i];
            Q = tour[Q_i];

            if (chrono::high_resolution_clock::now() > deadline) {
                return;
            }

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

                    int32_t delta = 0;

                    // 2-edge exchanges.
                    if (d[A][C] + d[B][D] < d[A][B] + d[C][D]) {
                        // Change AB, BC to AC, BD, but keep EF.
                        reverse(tour, B_i, C_i, position);
                        max = maximum(max, d[A][C], d[B][D]);
                        delta = (d[A][C] + d[B][D]) - (d[A][B] + d[C][D]);
                    } else if (d[C][E] + d[D][F] < d[C][D] + d[E][F]) {
                        // Change CD, EF to CE, DF, but keep AB.
                        reverse(tour, D_i, E_i, position);
                        max = maximum(max, d[C][E], d[D][F]);
                        delta = (d[C][E] + d[D][F]) - (d[C][D] + d[E][F]);
                    } else if (d[E][A] + d[F][B] < d[A][B] + d[E][F]) {
                        // Change AB, EF to EA, FB, but keep CD.
                        reverse(tour, F_i, A_i, position);
                        max = maximum(max, d[E][A], d[F][B]);
                        delta = (d[E][A] + d[F][B]) - (d[A][B] + d[E][F]);
                    } else {
                        // 3-edge exchanges.
                        uint32_t d_AB_CD_EF = d[A][B] + d[C][D] + d[E][F];
                        if (d[D][A] + d[F][B] + d[C][E] < d_AB_CD_EF) {
                            // Change AB, CD, EF to DA, FB, CE.
                            reverse(tour, F_i, A_i, position);
                            reverse(tour, D_i, E_i, position);
                            max = maximum(max, d[D][A], d[F][B], d[C][E]);
                            delta = (d[D][A] + d[F][B] + d[C][E]) - d_AB_CD_EF;
                        } else if (d[B][D] + d[E][A] + d[F][C] < d_AB_CD_EF) {
                            // Change AB, CD, EF to BD, EA, FC.
                            reverse(tour, F_i, A_i, position);
                            reverse(tour, B_i, C_i, position);
                            max = maximum(max, d[B][D], d[E][A], d[F][C]);
                            delta = (d[B][D] + d[E][A] + d[F][C]) - d_AB_CD_EF;
                        } else if (d[A][C] + d[B][E] + d[D][F] < d_AB_CD_EF) {
                            // Change AB, CD, EF to AC, BE, DF.
                            reverse(tour, B_i, C_i, position);
                            reverse(tour, D_i, E_i, position);
                            max = maximum(max, d[A][C], d[B][E], d[D][F]);
                            delta = (d[A][C] + d[B][E] + d[D][F]) - d_AB_CD_EF;
                        } else if (d[B][E] + d[D][A] + d[F][C] < d_AB_CD_EF) {
                            // Change AB, CD, EF to BE, DA, FC.
                            reverse(tour, A_i, F_i, position);
                            reverse(tour, B_i, C_i, position);
                            reverse(tour, D_i, E_i, position);
                            max = maximum(max, d[B][E], d[D][A], d[F][C]);
                            delta = (d[B][E] + d[D][A] + d[F][C]) - d_AB_CD_EF;
                        }
                    }

                    length += delta; // Update tour length.

                    if (delta != 0)
                        goto next_AB; // Go to next edge AB.
                }
            }
            next_AB: continue;
        }
    } while (length < oldLength);
}

/**
 * Approximates optimal TSP tour through graph read from the given input stream.
 *
 * The function will try to return before the given deadline, but expect some
 * variations.
 *
 * @param d Distance matrix.
 * @param deadline Time available in milliseconds.
 * @return An approximation of the optimal TSP tour.
 */
std::vector<uint16_t> approximate(
        istream &in,
        const chrono::time_point<chrono::high_resolution_clock>& deadline) {

    // 3-opt takes some time to finish.
    // Deadline before which last optimization step must finish.
    //auto optimizeDeadline = deadline - chrono::duration<int, milli>(500);
    auto realDeadline = deadline - chrono::duration<int, milli>(200);

    // Create distance / nearest neighbors matrix.
    Matrix<uint32_t> d = createDistanceMatrix(in);
    Matrix<uint16_t> neighbor = createNeighborsMatrix(d, OPTIMIZATION_MAX_NEIGHBORS);

    const size_t N = d.rows();           // Number of cities.
    const uint32_t min = minDistance(d); // Shortest inter-city distance.

    // Create vector of random possible start cities for greedy heuristic.
    vector<uint16_t> randomCity(N);
    for (uint16_t i = 0; i < N; ++i) {
        randomCity[i] = i;
    }
    shuffle(randomCity.begin(), randomCity.end(), rng);

    // Index and length of shortest found tour.
    size_t minTour = 0;
    uint64_t minLength = numeric_limits<uint64_t>::max();
    vector<vector<uint16_t> > tours; // Found tours.
    vector<uint16_t> position(N);    // Maps cities to their position in tour.

    chrono::duration<int, milli> totHeurDuration(0);
    chrono::duration<int, milli> avgHeurDuration(0);
    chrono::duration<int, milli> totOptDuration(0);
    chrono::duration<int, milli> avgOptDuration(0);
    
    // The main loop of the algorithm.
    int numIters = 0;
    for (size_t i = 0; i < N; ++i) {
        uint64_t tourLength; // Length of current tour.

        // Create greedy tour starting in a random city.
        tours.emplace_back(N);

        auto before = chrono::high_resolution_clock::now();
        tours[i] = greedy(d, randomCity[i], tourLength);
        auto duration = chrono::high_resolution_clock::now() - before;
        totHeurDuration += chrono::duration_cast<std::chrono::milliseconds>(duration);
        avgHeurDuration = totHeurDuration / (i + 1);

        // Maximum inter-city distance in tour.
        uint32_t max = maxDistance(tours[i], d);

        // Initialize city tour position vector.
        for (uint16_t j = 0; j < N; ++j) {
            position[tours[i][j]] = j; // tours[i][j] is the j:th city in tours[i].
        }

        // Optimize tour using 2- or 3-opt.
        before = chrono::high_resolution_clock::now();
        twoOpt(tours[i], d, neighbor, position, max, min, tourLength);
        //threeOpt(tours[i], d, neighbor, position, max, min, tourLength, realDeadline);
        duration = chrono::high_resolution_clock::now() - before;
        totOptDuration += chrono::duration_cast<std::chrono::milliseconds>(duration);
        avgOptDuration = totOptDuration / (i + 1);

        // Check if we got a shorter tour.
        if (tourLength < minLength) {
            minTour = i;
            minLength = tourLength;
        }

        numIters++;
        // Check if deadline has passed.
        if ((chrono::high_resolution_clock::now() + 2 * avgOptDuration) > deadline) {
            break;
        }
    }

    /*
    cerr << "numIters: " << numIters << endl;
    cerr << "totHeurDuration: " << chrono::duration_cast<std::chrono::milliseconds>(totHeurDuration).count() << " ms" << endl;
    cerr << "avgHeurDuration: " << chrono::duration_cast<std::chrono::milliseconds>(avgHeurDuration).count() << " ms" << endl;
    cerr << "totOptDuration: " << chrono::duration_cast<std::chrono::milliseconds>(totOptDuration).count() << " ms" << endl;
    cerr << "avgOptDuration: " << chrono::duration_cast<std::chrono::milliseconds>(avgOptDuration).count() << " ms" << endl;*/
    return tours[minTour]; // Return shortest found tour.
}

int main(int argc, char *argv[]) {
    // Establish a deadline 1900 ms into the future.
    auto deadline = chrono::high_resolution_clock::now();
    deadline += chrono::duration<int, milli>(1900);

    // Approximate a TSP tour.
    vector<uint16_t> tour = approximate(cin, deadline);

    // Print the tour.
    for (auto city : tour) {
        cout << city << endl;
    }

    return 0;
}


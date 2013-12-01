#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cmath>
#include <cstdint>

/**
 * Create a new N x N matrix and return it.
 */
uint32_t **newMatrix(std::size_t N) {
    uint32_t **M = new uint32_t*[N];
    for (std::size_t i = 0; i < N; ++i) {
        M[i] = new uint32_t[N];
        std::fill(M[i], M[i] + N, 0);
    }
    return M;
}

/**
 * Delete the N x N matrix M.
 */
void deleteMatrix(uint32_t **M, std::size_t N) {
    for (std::size_t i = 0; i < N; ++i) {
        delete M[i];
    }
    delete[] M;
}

/**
 * Create a distance matrix from an input stream and return it.
 *
 * @param in Input stream.
 * @param N The number of cities (output).
 * @return The read distance matrix.
 */
uint32_t **createDistanceMatrix(std::istream& in, std::size_t *N) {
    // Read number of vertices.
    in >> *N;

    // Read vertex coordinates.
    std::vector<double> x(*N);
    std::vector<double> y(*N);
    for (std::size_t i = 0; i < *N; ++i) {
        in >> x[i] >> y[i];
    }

    // Calculate distance matrix.
    uint32_t **d = newMatrix(*N);
    for (std::size_t i = 0; i < *N; ++i) {
        for (std::size_t j = 0; j < *N; ++j) {
            d[i][j] = std::round(std::sqrt(std::pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2)));
        }
    }

    return d;
}

/**
 * Prints the N x N matrix M to standard output.
 *
 * @param M The matrix to print.
 * @param N Number of vertices in matrix.
 */
void printMatrix(uint32_t **M, std::size_t N) {
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            std::cout << std::left << std::setw(3) << M[i][j];
        }
        std::cout << std::endl;
    }
}

#endif // MATRIX_H

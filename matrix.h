#ifndef MATRIX_H
#define MATRIX_H

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * A very simple N x M matrix class.
 */
template<typename T>
class Matrix {
public:
    Matrix(std::size_t N, std::size_t M) :
        m_data(N * M, T()),
        m_rows(N),
        m_cols(M) {}

    inline T* operator[](int i) {
        return &m_data[i * m_cols];
    }

    inline T const* operator[](int i) const {
        return &m_data[i * m_cols];
    }

    inline std::size_t rows() const {
        return m_rows;
    }

    inline std::size_t cols() const {
        return m_cols;
    }

private:
    std::vector<T> m_data;
    std::size_t m_rows;
    std::size_t m_cols;
};

// Stream output operator.
template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        for (std::size_t j = 0; j < matrix.cols(); ++j) {
            os << std::left << std::setw(3) << matrix[i][j];
        }
        os << std::endl;
    }
    return os;
}

#endif // MATRIX_H

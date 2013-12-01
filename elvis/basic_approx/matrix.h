#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstdint>

/**
 * A very simple N x N matrix class.
 */
template<typename T>
class Matrix {
public:
    Matrix(std::size_t N) : m_data(N * N, T()), m_size(N) {}

    inline T* operator[](int i) {
        return &m_data[i * m_size];
    }

    inline T const* operator[](int i) const {
        return &m_data[i * m_size];
    }

    inline std::size_t size() const {
        return m_size;
    }

private:
    std::vector<T> m_data;
    std::size_t m_size;
};

// Stream output operator.
template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        for (std::size_t j = 0; j < matrix.size(); ++j) {
            os << std::left << std::setw(3) << matrix[i][j];
        }
        os << std::endl;
    }
    return os;
}

#endif // MATRIX_H

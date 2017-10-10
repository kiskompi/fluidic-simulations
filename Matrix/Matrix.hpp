#ifndef MATRIX_H
#define MATRIX_H

#include "../kokkos/containers/src/Kokkos_Vector.hpp"
#include <array>

template<typename value_type>
class Matrix{

private:
    Kokkos::vector<std::Kokkos::vector<value_type>> m_matrix;

public:
    Matrix();
    ~Matrix(){}

    size_t get_sizeX() const {return sizeX;}
    size_t get_sizeY() const {return sizeY;}

    Matrix   operator=(const value_type &rval);
    std::array<value_type, sizeY>& operator[](const size_t i) {
        return m_matrix[i];
    }
    const std::array<value_type, sizeY>& operator[](const size_t i) const {
        return m_matrix[i];
    }

};

template<typename T, size_t X, size_t Y>
inline Matrix<T, X, Y>::Matrix():
    m_matrix(std::array<std::array<T, Y>, X>())
{
    for (auto arr: m_matrix){
        arr.fill(0);
    }
}
#endif
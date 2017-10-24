#ifndef MATRIX_H
#define MATRIX_H

#include <Kokkos_Vector.hpp>

template<typename value_type>
class Matrix{

private:
    Kokkos::vector<Kokkos::vector<value_type>> m_matrix;

public:
    Matrix();
    ~Matrix(){}

    size_t get_sizeX() const {return m_matrix.size();}
    size_t get_sizeY() const {return m_matrix[0].size();}

    Kokkos::vector<value_type>& operator=(const value_type &rval);
    Kokkos::vector<value_type>& operator[](const size_t i) {
        return m_matrix[i];
    }
    const Kokkos::vector<value_type>& operator[](const size_t i) const {
        return m_matrix[i];
    }

};

template<typename T>
inline Matrix<T>::Matrix():
    m_matrix (Kokkos::vector<Kokkos::vector<T>>())
{
    for (Kokkos::vector<T> arr: m_matrix){
        for (T elem: arr){
            elem = 0;
        }
    }
}
#endif
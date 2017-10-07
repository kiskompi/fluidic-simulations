#ifndef MATRIX_H
#define MATRIX_H

#include <array>

template<typename value_type, size_t sizeX, size_t sizeY>
class Matrix{
template<typename T, size_t size>
using STDArray = std::array<T, size>;

private:
    STDArray<STDArray<value_type, sizeY>, sizeX> m_matrix;

public:
    Matrix();
    ~Matrix(){}

    size_t get_sizeX() const {return sizeX;}
    size_t get_sizeY() const {return sizeY;}

    Matrix   operator=(const value_type &rval);
    STDArray<value_type, sizeX> operator[](const size_t i);

};

template<typename T, size_t X, size_t Y>
inline Matrix<T, X, Y>::Matrix():
    m_matrix()
{
    for (auto arr: m_matrix){
        arr.fill(0);
    }
}
#endif
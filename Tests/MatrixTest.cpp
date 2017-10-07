#include "../Matrix.hpp"
template <typename T, size_t X, size_t Y>
void MatrixTest(){
    Matrix<T, X, Y> a = Matrix<T, X, Y>();
}

int main(){
    MatrixTest<int, 100, 100>();
}

#include "../Matrix.hpp"
template <typename T, size_t X, size_t Y>
void MatrixTest(){
    Matrix<T, X, Y> a = Matrix<T, X, Y>();
}
template <typename T, size_t X, size_t Y>
void ViewTest(){
   Kokkos::View<T[X][Y], Kokkos::OpenMP> b;
}
int main(){
    MatrixTest<int, 100, 100>();
    ViewTest<double, 100, 100>();
}

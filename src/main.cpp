#include "matrix_t.h"

const double eps0 = pow(10, -3);

int main() {
    int n = 3;
    Matrix_t A(n, n);
    Matrix_t b(n, 1);
    A.Enter(); b.Enter();
    SeidelMethod(A, b, eps0).first.Print();
}
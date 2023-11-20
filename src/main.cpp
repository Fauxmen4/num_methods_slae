#include "matrix_t.h"

std::pair<Matrix_t, int> SimpleIterMethod(Matrix_t A, Matrix_t b, double eps);
/*std::pair<Matrix_t, int>*/ void SeidelMethod(Matrix_t A, Matrix_t b, double eps);

const double eps0 = pow(10, -3);

int main() {
    int n = 3;
    Matrix_t A(n, n);
    Matrix_t b(n, 1);
    A.Enter(); b.Enter();
    SeidelMethod(A, b, eps0);
}

//! Добавить проверку на норму и соответсвующие действия если проверка не проходится
/*std::pair<Matrix_t, int>*/ void SeidelMethod(Matrix_t A, Matrix_t b, double eps) {
    int n = b.lines;
    Matrix_t d(n, 1);
    for (int i = 0; i < n; i++) {
        d.matrix[i][0] = b.matrix[i][0]/A.matrix[i][i];
    }
    Matrix_t c(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            } else {
                c.matrix[i][j] = (-1)*A.matrix[i][j]/A.matrix[i][i];
            }
        }
    }
    
}

//! Добавить проверку на норму и соответсвующие действия если проверка не проходится
std::pair<Matrix_t, int> SimpleIterMethod(Matrix_t A, Matrix_t b, double eps) {
    int n = b.lines;
    double mu = 1/A.InfiniteNorm();
    Matrix_t B(n); B = B - mu*A;
    Matrix_t c(n, 1); c = mu*b;
    Matrix_t x0(c);
    Matrix_t x1(n, 1); x1 = B*x0+c;
    int count = 1;
    while (B.InfiniteNorm()/(1-B.InfiniteNorm())*(x1-x0).InfiniteNorm() > eps) {
        x0 = x1;
        x1 = B*x0+c;
        count += 1;
    }
    return std::pair<Matrix_t, int>{x1, count};
}
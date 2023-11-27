#include "matrix_t.h"

//! Добавить проверку на норму и соответсвующие действия если проверка не проходится
std::pair<Matrix_t, int> SeidelMethod(Matrix_t A, Matrix_t b, double eps) {
    int n = b.lines;
    Matrix_t d(n, 1);
    for (int i = 0; i < n; i++) {
        d.matrix[i][0] = b.matrix[i][0] / A.matrix[i][i];
    }
    Matrix_t C(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                C.matrix[i][j] = 0;
            else 
                C.matrix[i][j] = (-1)*A.matrix[i][j]/A.matrix[i][i];
        }
    }
    Matrix_t x0(d);
    Matrix_t x1(n, 1);
    int count = 0;
    while ((A*x1-b).InfNorm() > eps) {
        if (count != 0) {
            x0 = x1;
        }
        int zeroIndex = 0;
        x1.Clear();
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if (j < zeroIndex)
                    x1.matrix[i][0] += C.matrix[i][j]*x1.matrix[j][0];
                else if (j == zeroIndex)
                    continue;
                else 
                    x1.matrix[i][0] += C.matrix[i][j]*x0.matrix[j][0];
            }
            x1.matrix[i][0] += d.matrix[i][0];
            zeroIndex++;
        }
        count++;
    } 
    return std::pair<Matrix_t, int>{x1, count};
}


//! Добавить проверку на норму и соответсвующие действия если проверка не проходится
std::pair<Matrix_t, int> SimpleIterMethod(Matrix_t A, Matrix_t b, double eps) {
    int n = b.lines;
    double mu = 1/A.InfNorm();
    Matrix_t B(n); B = B - mu*A;
    Matrix_t c(n, 1); c = mu*b;
    Matrix_t x0(c);
    Matrix_t x1(n, 1);
    int count = 0;
    while (B.InfNorm()/(1-B.InfNorm())*(x1-x0).InfNorm() > eps) {
        if (count != 0) {
            x0 = x1;
        }
        x1 = B*x0+c;
        count += 1;
    }
    return std::pair<Matrix_t, int>{x1, count};
}
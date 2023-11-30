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
    double mu = 1/A.FirstNorm();
    Matrix_t B(n); B = B - mu*A;
    Matrix_t c(n, 1); c = mu*b;
    Matrix_t x0(c);
    Matrix_t x1(n, 1);
    int count = 0;
    while (B.FirstNorm()/(1-B.FirstNorm())*(x1-x0).FirstNorm() > eps) {
        if (count != 0) {
            x0 = x1;
        }
        x1 = B*x0+c;
        count += 1;
    }
    return std::pair<Matrix_t, int>{x1, count};
}

Matrix_t QR(Matrix_t A, Matrix_t b) {
    int n = A.lines;
    Matrix_t Q(n);
    Matrix_t R(A);
    
    for (int i = 0; i < n-1; i++) {
        Matrix_t z(Ort(n-i));
        Matrix_t y(R.Cut(i,n-1,i,i));
        double a = y.SecondVectorNorm();
        Matrix_t w = (y-a*z)*(1/(y-a*z).SecondVectorNorm());
        Matrix_t Q1(n);
        Q1.Insert(i,n-1,i,n-1,E(n-i)-2*w*w.Transpose());
        R.Insert(i,n-1,i,n-1,(Q1.Cut(i,n-1,i,n-1)*R.Cut(i,n-1,i,n-1)));
        Q = Q*Q1;
    }

    Matrix_t y(Q.Transpose()*b);
    Matrix_t x(n, 1);
    x.matrix[n-1][0] = y.matrix[n-1][0]/R.matrix[n-1][n-1];
    for (int i = n-1; i >= 0; i--) {
        x.matrix[i][0] = y.matrix[i][0];
        for (int j = n-1; j > i; j--) {
            x.matrix[i][0] -= R.matrix[i][j]*x.matrix[j][0];
        }
        x.matrix[i][0] /= R.matrix[i][i];
    }

    return x;
}

Matrix_t LU(Matrix_t A, Matrix_t b) {
    int n = A.lines;
    Matrix_t P(n);
    Matrix_t M(A);
    for (int i = 0; i < n-1; i++)
    {
        double max_el = 0;
        int lineToSwap = -1;
        for (int j = i; j < n; j++) {
            if (abs(M.matrix[j][i]) > abs(max_el)) {
                max_el = std::max(M.matrix[j][i], max_el);
                lineToSwap = j;
            }
        }
        M.SwapLines(i, lineToSwap);
        P.SwapLines(i, lineToSwap);
        //? Теперь преобразуем M:
            //? Сначала пройдемся по i-тому столбцу
            //? Потом по элементам справа от элементов под позицией [i][i]
        for(int j = i+1; j < n; j++) {
            M.matrix[j][i] = M.matrix[j][i]/M.matrix[i][i];
        }    
        for(int j = i+1; j < n; j++) {
            for (int k = i+1; k < n; k++) {
                M.matrix[j][k] = M.matrix[j][k] - M.matrix[j][i]*M.matrix[i][k];
            }
        }
    }
    Matrix_t U(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            //std::cout << M.matrix[i][j] << ' ';
            U.matrix[i][j] = M.matrix[i][j];
        }
    }
    Matrix_t L(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            L.matrix[i][j] = M.matrix[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        L.matrix[i][i] = 1;
    }
    //? L и U готовы
    //? Решение уравнения Ly = Pb
    b = P*b;
    Matrix_t y(n,1);
    y.matrix[0][0] = b.matrix[0][0];
    for (int i = 1; i < n; i++) {
        y.matrix[i][0] = b.matrix[i][0];
        for (int j = 0; j < i; j++) {
            y.matrix[i][0] -= L.matrix[i][j]*y.matrix[j][0];
        }
    }
    //? Решение уравнения Ux = y
    Matrix_t x(n, 1);
    x.matrix[n-1][0] = y.matrix[n-1][0]/U.matrix[n-1][n-1];
    for (int i = n-1; i >= 0; i--) {
        x.matrix[i][0] = y.matrix[i][0];
        for (int j = n-1; j > i; j--) {
            x.matrix[i][0] -= U.matrix[i][j]*x.matrix[j][0];
        }
        x.matrix[i][0] /= U.matrix[i][i];
    }
    return x;
}
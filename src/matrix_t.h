#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <tuple>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

class Matrix_t 
{
    friend Matrix_t operator * (double s, const Matrix_t& other);

public:
    int lines;
    int columns;
    double** matrix;

public:
    Matrix_t(int lines, int columns);
    Matrix_t(const Matrix_t &other);
    Matrix_t(int n);
    ~Matrix_t();

    void Enter();
    void Print();

    bool operator == (const Matrix_t& other);
    bool operator != (const Matrix_t& other);

    Matrix_t& operator = (const Matrix_t& other);
    
    Matrix_t operator + (const Matrix_t& other);
    Matrix_t operator - (const Matrix_t& other);
    Matrix_t operator * (const Matrix_t& other);
    Matrix_t operator * (double s);

    Matrix_t Transpose();

    int Det();      //! Not implemented yet

    double FirstNorm();
    double SecondVectorNorm();   
    double InfNorm();

    bool NotNull();

    void Clear();

    Matrix_t Cut(int b1, int b2, int a1, int a2);
    void Insert(int b1, int b2, int a1, int a2, const Matrix_t& other);
    void SwapLines(int l1, int l2);
};

Matrix_t Ort(int n);
Matrix_t E(int n);

std::pair<Matrix_t, int> SimpleIterMethod(Matrix_t A, Matrix_t b, double eps);
std::pair<Matrix_t, int> SeidelMethod(Matrix_t A, Matrix_t b, double eps);

Matrix_t Gauss(Matrix_t A, Matrix_t b);
Matrix_t QR(Matrix_t A, Matrix_t b);
Matrix_t LU(Matrix_t A, Matrix_t b);
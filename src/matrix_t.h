#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <tuple>

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
    double SecondNorm();   
    double InfNorm();

    bool NotNull();

    void Clear();
};

std::pair<Matrix_t, int> SimpleIterMethod(Matrix_t A, Matrix_t b, double eps);
std::pair<Matrix_t, int> SeidelMethod(Matrix_t A, Matrix_t b, double eps);
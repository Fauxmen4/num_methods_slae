#include "matrix_t.h"

Matrix_t::Matrix_t(int lines, int columns) {
    matrix = new double*[lines];
    for (int i = 0; i < lines; i++) {
        matrix[i] = new double[columns]; 
    }
    this->lines = lines;
    this->columns = columns;
}

Matrix_t::Matrix_t(const Matrix_t &other) {
    this->lines = other.lines;
    this->columns = other.columns;
    this->matrix = new double*[lines];
    for(int i = 0; i < this->lines; i++) {
        matrix[i] = new double[columns];
    }
    for (int i = 0; i < this->lines; i++) {
        for (int j = 0; j < this->columns; j++) {
            this->matrix[i][j] = other.matrix[i][j];
        }
    }
}

Matrix_t::Matrix_t(int n) {
    this->lines = n;
    this->columns = n;
    this->matrix = new double*[lines];
    for(int i = 0; i < this->lines; i++) {
        matrix[i] = new double[columns];
    }
    for (int i = 0; i < lines; i++) {
        this->matrix[i][i] = 1;
    }
}

Matrix_t::~Matrix_t() {
    for (int i = 0; i < this->lines; i++) {
        delete [] this->matrix[i];
    }
    delete [] this->matrix;
}

void Matrix_t::Print() {
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < columns; j++) {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

void Matrix_t::Enter() {
    this->matrix = new double*[lines];
    for (int i = 0; i < this->lines; i++) {
        this->matrix[i] = new double[columns];
    }
    for (int i = 0; i < this->lines; i++) {
        for (int j = 0; j < this->columns; j++) {
            int tmp; std::cin >> tmp; 
            this->matrix[i][j] = tmp;
        }
    }
}

bool Matrix_t::operator == (const Matrix_t& other) {
    if (this->lines != other.lines || this->columns != other.columns) {
        return false;
    }
    for (int i = 0; i < this->lines; i++) {
        for (int j = 0; j < this->columns; j++) {
            if (this->matrix[i][j] != other.matrix[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix_t::operator != (const Matrix_t& other) {
    if (this->lines != other.lines || this->columns != other.columns) {
        return true;
    }
    for (int i = 0; i < this->lines; i++) {
        for (int j = 0; j < this->columns; j++) {
            if (this->matrix[i][j] != other.matrix[i][j]) {
                return true;
            }
        }
    }
    return false;
}

Matrix_t& Matrix_t::operator = (const Matrix_t& other) {
    for (int i = 0; i < this->lines; i++) {
        delete [] this->matrix[i];
    }
    delete [] this->matrix;
    this->lines = other.lines;
    this->columns = other.columns;
    this->matrix = new double*[other.lines];
    for(int i = 0; i < other.lines; i++) {
        this->matrix[i] = new double[other.columns];
    }
    for(int i = 0; i < this->lines; i++) {
        for (int j = 0; j < other.lines; j++) {
            this->matrix[i][j] = other.matrix[i][j];
        }
    }
    return *this;
}

Matrix_t Matrix_t::operator + (const Matrix_t& other) {
    Matrix_t tmp(this->lines, this->columns);
    for (int i = 0; i < tmp.lines; i++) {
        for (int j = 0; j < tmp.columns; j++) {
            tmp.matrix[i][j] = other.matrix[i][j] + this->matrix[i][j];
        }
    }
    return tmp;
}

Matrix_t Matrix_t::operator - (const Matrix_t& other) {
    Matrix_t tmp(this->lines, this->columns);
    for (int i = 0; i < tmp.lines; i++) {
        for (int j = 0; j < tmp.columns; j++) {
            tmp.matrix[i][j] = this->matrix[i][j] - other.matrix[i][j];
        }
    }
    return tmp;
}

Matrix_t Matrix_t::operator * (double s) {
    Matrix_t tmp(this->lines, this->columns);
    for (int i = 0; i < this->lines; i++) {
        for (int j = 0; j < this->columns; j++) {
            tmp.matrix[i][j] = this->matrix[i][j]*s;
        }
    }
    return tmp;
}

Matrix_t operator * (double s, const Matrix_t& other) {
    Matrix_t tmp(other.lines, other.columns);
    for (int i = 0; i < other.lines; i++) {
        for (int j = 0; j < other.columns; j++) {
            tmp.matrix[i][j] = other.matrix[i][j]*s;
        }
    }
    return tmp;
}

Matrix_t Matrix_t::operator * (const Matrix_t& other) {
    Matrix_t tmp(this->lines, other.columns);
    for (int i = 0; i < tmp.lines; i++) {
        for (int j = 0; j < tmp.columns; j++) {
            for (int k = 0; k < this->columns; k++) {
                tmp.matrix[i][j] += this->matrix[i][k]*other.matrix[k][j];
            }
        }
    }
    return tmp;
}

Matrix_t Matrix_t::Transpose() {
    Matrix_t tmp(this->columns, this->lines);
    for (int i = 0; i < this->lines; i++) {
        for (int j = 0; j < this->columns; j++) {
            tmp.matrix[j][i] = this->matrix[i][j];
        }
    }
    return tmp;
}

//! Not implemented yet
int Matrix_t::Det() {
    return 0;   
}

double Matrix_t::FirstNorm() {
    double res = -100000;
    for (int i = 0; i < this->columns; i++) {
        double tmp = 0;
        for (int j = 0; j < this->lines; j++) {
            tmp += abs(this->matrix[j][i]);
        }
        res =std::max(tmp, res);
    }
    return res;   
}

double Matrix_t::SecondNorm() {
    return 0;
}

double Matrix_t::InfiniteNorm() {
    double res = -100000;
    for (int i = 0; i < this->lines; i++) {
        double tmp = 0;
        for (int j = 0; j < this->columns; j++) {
            tmp += abs(this->matrix[i][j]);
        }
        res =std::max(tmp, res);
    }
    return res;
}

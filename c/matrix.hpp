#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <tuple>
using namespace std;

class Matrix {
    public:
    unsigned int rows;
    unsigned int cols;
    vector<vector<double>> M;

    Matrix();
    Matrix(vector<vector<double>> A);
    Matrix(unsigned int r, unsigned int c);
    Matrix(unsigned int n);

    vector<double>& operator [](unsigned int i);
    Matrix operator *(const Matrix& B);
    vector<double> operator *(const vector<double>& B);
    vector<double> get_row(unsigned int i);
    vector<double> get_col(unsigned int j);

    Matrix copy();
    Matrix t();
    Matrix i();
    double det();
    bool isSPD();

    private:
    char spd_flag = 0;
};


Matrix transpose(Matrix A);
Matrix prod(Matrix A, Matrix B);
double scal_prod(vector<double> x, vector<double> y);
double determinant(Matrix A);

tuple<Matrix, Matrix, Matrix> lu(Matrix A);
Matrix inv_lu(Matrix A);

Matrix cholesky(Matrix S);
Matrix inv_chol(Matrix A);

Matrix inverse(Matrix A);

vector<double> solve_upper(Matrix U, vector<double> b);
vector<double> solve_lower(Matrix L, vector<double> b);

#endif
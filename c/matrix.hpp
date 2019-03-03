#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <tuple>
using namespace std;

class Matrix {
    public:
    unsigned int rows;
    unsigned int cols;
	double *M;

	Matrix();
	Matrix(const Matrix& A);
	Matrix(vector<vector<double>> A);
	Matrix(unsigned int r, unsigned int c);
	Matrix(unsigned int n);
	~Matrix();

	double& operator ()(unsigned int i, unsigned int j);
	Matrix operator +(Matrix B);
	Matrix operator -(Matrix B);
	Matrix operator *(Matrix& B);
	Matrix operator *(Matrix B);

    Matrix& operator =(Matrix B);

    vector<double> get_row(unsigned int i);
    vector<double> get_col(unsigned int j);
    double& at(unsigned int i, unsigned int j);

    Matrix copy();
    Matrix t();
    Matrix i();
    double det();
    bool isSPD();

    private:
    char spd_flag = 0;
};

Matrix to_row(vector<double> x);
Matrix to_col(vector<double> x);
Matrix transpose(Matrix& A);
Matrix prod(Matrix& A, Matrix& B);
double determinant(Matrix& A);

tuple<Matrix*, Matrix*, Matrix*> lu(Matrix& A);
Matrix inv_lu(Matrix& A);

Matrix cholesky(Matrix& S);
Matrix inv_chol(Matrix& A);

Matrix inverse(Matrix& A);

vector<double> solve_upper(Matrix *U, vector<double> b);
vector<double> solve_lower(Matrix *L, vector<double> b);

#endif
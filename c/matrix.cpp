#include <math.h>
#include "matrix.hpp"
#define TRUE 'T'
#define FALSE 'F'

/*
    Matrix class implementation
    Provides useful methods to manage 2D arrays and access their data.
*/
Matrix::Matrix() {}
Matrix::Matrix(const Matrix& A) {
	rows = A.rows;
	cols = A.cols;
	M = new double[rows*cols];
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
			M[i*cols + j] = A.M[i*cols + j];
}
Matrix::Matrix(vector<vector<double>> A) {
	rows = A.size();
	cols = A[0].size();

	M = new double[rows*cols];
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
			M[i*cols + j] = A[i][j];
}
Matrix::Matrix(unsigned int r, unsigned int c) {
	rows = r;
	cols = c;
	M = new double[rows*cols];
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
			M[i*cols + j] = 0;
}
Matrix::Matrix(unsigned int n) {
	rows = n;
	cols = n;
	M = new double[n*n];
	for (unsigned int i = 0; i < n; i++)
		for (unsigned int j = 0; j < n; j++)
			M[i*cols + j] = 0;
}
Matrix::~Matrix() {
	delete[] M;
}

double& Matrix::operator ()(unsigned int i, unsigned int j) {
	return M[i*cols + j];
}

////////
Matrix Matrix::operator +(Matrix B) {
	Matrix C = Matrix(rows,cols);
	if (B.rows != rows || B.cols != cols)	return C;

	for (unsigned long i = 0; i < rows*cols; i++)
		C.M[i] = M[i] + B.M[i];
	return C;
}
Matrix Matrix::operator -(Matrix B) {
	Matrix C = Matrix(rows,cols);
	if (B.rows != rows || B.cols != cols)	return C;

	for (unsigned long i = 0; i < rows*cols; i++)
		C.M[i] = M[i] - B.M[i];
	return C;
}

Matrix& Matrix::operator =(Matrix B) {
	if (rows != B.rows || cols != B.cols) {
		delete[] M;
		rows = B.rows;
		cols = B.cols;
		M = new double[rows*cols];
	}

	for (unsigned int i = 0; i < rows; i++)
        for (unsigned int j = 0; j < cols; j++)
    		M[i*cols+j] = B.M[i*cols+j];
	return *(this);
}
////////

Matrix Matrix::operator *(Matrix& B) {
    return prod(*(this), B);
}
Matrix Matrix::operator *(Matrix B) {
    return prod(*(this), B);
}

vector<double> Matrix::get_row(unsigned int i) {
	vector<double> r(cols);
	for (unsigned int j = 0; j < cols; j++)
		r[j] = M[i*cols + j];
	return r;
}

vector<double> Matrix::get_col(unsigned int j) {
	vector<double> c(rows);
	for (unsigned int i = 0; i < rows; i++)
		c[i] = M[i*cols + j];
	return c;
}

double& Matrix::at(unsigned int i, unsigned int j) {
	return M[i*cols + j];
}

Matrix Matrix::copy() {
	Matrix copy = Matrix();
	copy.rows = rows;
	copy.cols = cols;
	copy.M = new double[rows*cols];
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
			copy(i, j) = M[i*cols + j];
	return copy;
}
Matrix Matrix::t() {
    return transpose(*(this));
}
Matrix Matrix::i() {
    return inverse(*(this));
}
double Matrix::det() {
    return determinant(*(this));
}
bool Matrix::isSPD() {
    if (spd_flag == TRUE)  return true;
    if (spd_flag == FALSE)  return false;

    // Symmetric
    for (unsigned int i = 0; i < rows; i++)
        for (unsigned int j = 0; j < i; j++)
            if (M[i*cols + j] != M[j*cols + i]) {
                spd_flag = FALSE;
                return false;
            }
    
    // Positive definite (x'Mx > 0)
    vector<double> x = {1, 1, 1};
	Matrix Xr = to_row(x);
	Matrix Xc = to_col(x);
	Matrix T = prod(*(this), Xc);
	Matrix sp = prod(Xr, T);

    if (sp.M[0] > 0) {
        spd_flag = TRUE;
        return true;
    } else {
        spd_flag = FALSE;
        return false;
    }
}


/*
    Algorithms and operations with matrices
    Transpose:  O = A'
    Product:    O = A*B
    Sc. prod:   o = x'*y
*/
Matrix transpose(Matrix& A) {
	Matrix tA = Matrix(A.cols, A.rows);
	for (unsigned int i = 0; i < tA.rows; i++)
		for (unsigned int j = 0; j < tA.cols; j++)
			tA(i, j) = A(j, i);
	return tA;
}

Matrix prod(Matrix& A, Matrix& B) {
	Matrix C = Matrix(A.rows, B.cols);
	for (unsigned int i = 0; i < A.rows; i++)
		for (unsigned int j = 0; j < B.cols; j++)
			for (unsigned int k = 0; k < A.cols; k++)
				C(i, j) += A(i, k) * B(k, j);
	return C;
}

Matrix to_row(vector<double> x) {
	return Matrix({ x });
}

Matrix to_col(vector<double> x) {
	unsigned int n = x.size();
	Matrix X = Matrix(n, 1);
	for (unsigned int i = 0; i < n; i++)
		X(i, 0) = x[i];
	return X;
}


/*
    Back Substitution method to obtain a upper triangular matrix
    The determinant is just the product of the elements on the diagonal
*/
double determinant(Matrix& A) {
	unsigned int i, j, k, n = A.rows;
	double l, max, det = 1;
	Matrix U = A.copy();

	for (k = 0; k < n - 1; k++) {
		max = 0;
		for (i = k; i < n; i++) {
			if (abs(U(i, k)) > max) {
				j = i;
				max = abs(U(i, k));
			}
		}

		if (j != k) {
			for (i = k; i < n; i++) {
				swap(U(j, i), U(k, i));
				det *= -1;
			}
		}

		for (i = k + 1; i < n; i++) {
			l = U(i, k) / U(k, k);
			for (j = k; j < n; j++)
				U(i, j) -= l * U(k, j);
		}
	}

	for (i = 0; i < n; i++)
		det = det * U(i, i);
	return det;
}

/*
    LU factorization (with partial pivoting): decomposes the input matrix A such that
        LU = PA     where L and U are respectively lower and upper triangular matrices.
    Computational cost:  O(2/3 n^3)
*/
tuple<Matrix*, Matrix*, Matrix*> lu(Matrix& A) {
	unsigned int n = A.rows;
	unsigned int i, j, k;
	double max = 0;

	Matrix *U = new Matrix(A);
	Matrix *L = new Matrix(n);
	Matrix *P = new Matrix(n);

	for (i = 0; i < n; i++) {
		L->at(i, i) = 1;
		P->at(i, i) = 1;
	}

	for (k = 0; k < n - 1; k++) {
		max = 0;
		for (i = k; i < n; i++) {
			if (abs(U->at(i, k)) > max) {
				j = i;
				max = abs(U->at(i, k));
			}
		}

		for (i = 0; i < n; i++) {
			if (i >= k)  swap(U->at(j, i), U->at(k, i));
			swap(P->at(k, i), P->at(j, i));
		}

		for (i = k + 1; i < n; i++) {
			L->at(i, k) = U->at(i, k) / U->at(k, k);
			for (j = k; j < n; j++)
				U->at(i, j) -= L->at(i, k) * U->at(k, j);
		}
	}

	return make_tuple(L, U, P);
}

/*
    Matrix inversion using the result of the LU factorization.
    Passing a matrix with range = 0 will produce an indefinite output.
    Computational cost:  O(8/3 n^3)
*/
Matrix inv_lu(Matrix& A) {
	unsigned int n = A.rows;

	Matrix *L, *U, *P;
	Matrix iA = Matrix(n);
	tie(L, U, P) = lu(A);
	vector<double> x(n, 0);

	for (unsigned int j = 0; j < n; j++) {
		x = solve_lower(L, P->get_col(j));
		x = solve_upper(U, x);
		for (unsigned int i = 0; i < n; i++)
			iA(i, j) = x[i];
	}

	delete L;
	delete U;
	delete P;
	return iA;
}

/*
    Similar to the LU: for SPD (Symmetric Positive Definite) matrices we have
        LL' = A     where L is a non-singular lower triangular matrix
*/
Matrix cholesky(Matrix& S) {
    unsigned int n = S.rows;
    Matrix C = Matrix(n);
    
    double sum;
    for (unsigned int i = 0; i < n; i++) { 
        for (unsigned int j = 0; j <= i; j++) { 
            sum = 0;
            for (unsigned int k = 0; k < j; k++) 
                sum += C(i,k) * C(j,k);
            
            if (i == j)   C(i,i) = sqrt(S(i,i) - sum); 
            else   C(i,j) = (S(i,j) - sum) / C(j,j); 
        }
    }

    return C;
}

/*
    This technique speeds up the inverse computation using the matrix symmetry.
    Computational cost:  O(1/3 n^3)
*/
Matrix inv_chol(Matrix& A) {
    unsigned int n = A.rows;
    Matrix cA = Matrix(n);
    vector<double> x(n);

    Matrix *C = new Matrix(cholesky(A));
	Matrix *tC = new Matrix(C->t());
    
    for (unsigned int j = 0; j < n; j++) {
        x = vector<double>(n, 0.);
        x[j] = 1;
        x = solve_lower(C, x);
        x = solve_upper(tC, x);
        for (unsigned int i = 0; i < n; i++)
            cA(i,j) = x[i];
    }
	
	delete C;
	delete tC;
    return cA;
}

Matrix inverse(Matrix& A) {
    if (A.isSPD())  return inv_chol(A);
    else  return inv_lu(A);
}


/*
    Linear systems solver for triangular matrices.
    Used only for matrix inversion, better approaches can be found in linear.cpp
*/
vector<double> solve_upper(Matrix *U, vector<double> b) {
	unsigned int i, j, n = U->rows;
	vector<double> x(n, 0);
	double sp;

	x[n-1] = b[n-1] / U->at(n-1, n-1);
	for (i = n - 2; i < n; i--) {
		sp = 0;
		for (j = i + 1; j < n; j++)
			sp += U->at(i,j)*x[j];

		x[i] = (b[i] - sp) / U->at(i,i);
	}
	return x;
}

vector<double> solve_lower(Matrix *L, vector<double> b) {
	unsigned int i, j, n = L->rows;
	vector<double> x(n, 0);
	double sp;

	x[0] = b[0] / L->at(0,0);
	for (i = 1; i < n; i++) {
		sp = 0;
		for (j = 0; j < i; j++)
			sp += L->at(i,j)*x[j];

		x[i] = (b[i] - sp) / L->at(i,i);
	}
	return x;
}

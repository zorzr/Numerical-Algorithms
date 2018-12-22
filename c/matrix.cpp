#include <math.h>
#include "matrix.hpp"
#define SWAP(a,b)  do { typeof(a) tmp = a;   a = b;   b = tmp; } while(0);
#define TRUE 'T'
#define FALSE 'F'

/*
    Matrix class implementation
    Provides useful methods to manage 2D arrays and access their data.
*/
Matrix::Matrix(){}
Matrix::Matrix(vector<vector<double>> A) {
    rows = A.size();
    cols = A[0].size();
    M = A;
}
Matrix::Matrix(unsigned int r, unsigned int c) {
    rows = r;
    cols = c;
    M = vector<vector<double>>(r, vector<double>(c, 0.));
}
Matrix::Matrix(unsigned int n) {
    rows = n;
    cols = n;
    M = vector<vector<double>>(n, vector<double>(n, 0.));
}

vector<double>& Matrix::operator [](unsigned int i) {
    return M[i];
}

Matrix Matrix::operator *(const Matrix& B) {
    return prod(*(this), B);
}
vector<double> Matrix::operator *(const vector<double>& b) {
    vector<double> x(rows);
    for (unsigned int i = 0; i < rows; i++)
        x[i] = scal_prod(M[i], b);
    return x;
}

vector<double> Matrix::get_row(unsigned int i) {
    return M[i];
}

vector<double> Matrix::get_col(unsigned int j) {
    vector<double> c(rows);
    for (unsigned int i = 0; i < rows; i++)
        c[i] = M[i][j];
    return c;
}

Matrix Matrix::copy() {
    return Matrix(M);
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
            if (M[i][j] != M[j][i]) {
                spd_flag = FALSE;
                return false;
            }
    
    // Positive definite (x'Mx > 0)
    vector<double> x = {1, 1, 1};
    if (scal_prod(x, (*(this))*x) > 0) {
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
Matrix transpose(Matrix A) {
    Matrix tA = Matrix(A.cols, A.rows);
    for (unsigned int i = 0; i < tA.rows; i++)
        for (unsigned int j = 0; j < tA.cols; j++)
            tA[i][j] = A[j][i];
    return tA;
}

Matrix prod(Matrix A, Matrix B) {
    Matrix C = Matrix(A.rows, B.cols);
    for (unsigned int i = 0; i < A.rows; i++)
        for (unsigned int j = 0; j < B.cols; j++)
            for (unsigned int k = 0; k < A.cols; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

double scal_prod(vector<double> x, vector<double> y) {
    if (x.size() != y.size())  return 0;
    
    double sp = 0;
    for (unsigned int i = 0; i < x.size(); i++)
        sp += x[i] * y[i];
    return sp;
}

/*
    Back Substitution method to obtain a upper triangular matrix
    The determinant is just the product of the elements on the diagonal
*/
double determinant(Matrix A) {
    unsigned int i, j, k, n = A.rows;
    double l, max, det = 1;
    Matrix U = A.copy();

    for (k = 0; k < n-1; k++) {
        max = 0;
        for (i = k; i < n; i++) {
            if (abs(U[i][k]) > max) {
                j = i;
                max = abs(U[i][k]);
            }
        }
        
        for (i = k; i < n; i++)
            SWAP(U[j][i], U[k][i])
        
        for (i = k+1; i < n; i++) {
            l = U[i][k] / U[k][k];
            for (j = k; j < n; j++)
                U[i][j] -= l * U[k][j];
        }
    }

    for (i = 0; i < n; i++)
        det = det * U[i][i];
    return det;
}

/*
    LU factorization (with partial pivoting): decomposes the input matrix A such that
        LU = PA     where L and U are respectively lower and upper triangular matrices.
    Computational cost:  O(2/3 n^3)
*/
tuple<Matrix, Matrix, Matrix> lu(Matrix A) {
    if (A.rows != A.cols)  return make_tuple(NULL,NULL,NULL);
    unsigned int n = A.rows;
    unsigned int i, j, k;
    double max = 0;

    Matrix U = A.copy();
    Matrix L = Matrix(n);
    Matrix P = Matrix(n);

    for (i = 0; i < n; i++) {
        L[i][i] = 1;
        P[i][i] = 1;
    }

    for (k = 0; k < n-1; k++) {
        max = 0;
        for (i = k; i < n; i++) {
            if (abs(U[i][k]) > max) {
                j = i;
                max = abs(U[i][k]);
            }
        }
        
        for (i = 0; i < n; i++) {
            if (i >= k)  SWAP(U[j][i], U[k][i])
            SWAP(P[k][i], P[j][i])
        }
        
        for (i = k+1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (j = k; j < n; j++)
                U[i][j] -= L[i][k] * U[k][j];
        }
    }

    return make_tuple(L, U, P);
}

/*
    Matrix inversion using the result of the LU factorization.
    Passing a matrix with range = 0 will produce an indefinite output.
    Computational cost:  O(8/3 n^3)
*/
Matrix inv_lu(Matrix A) {
    unsigned int n = A.rows;

    Matrix L, U, P;
    Matrix iA = Matrix(n);
    tie(L, U, P) = lu(A);
    vector<double> x(n);
    
    for (unsigned int j = 0; j < n; j++) {
        x = solve_lower(L, P.get_col(j));
        x = solve_upper(U, x);
        for (unsigned int i = 0; i < n; i++)
            iA[i][j] = x[i];
    }

    return iA;
}

/*
    Similar to the LU: for SPD (Symmetric Positive Definite) matrices we have
        LL' = A     where L is a non-singular lower triangular matrix
*/
Matrix cholesky(Matrix S) {
    unsigned int n = S.rows;
    Matrix C = Matrix(n);
    
    double sum;
    for (unsigned int i = 0; i < n; i++) { 
        for (unsigned int j = 0; j <= i; j++) { 
            sum = 0;
            for (unsigned int k = 0; k < j; k++) 
                sum += C[i][k] * C[j][k];
            
            if (i == j)   C[i][i] = sqrt(S[i][i] - sum); 
            else   C[i][j] = (S[i][j] - sum) / C[j][j]; 
        }
    }

    return C;
}

/*
    This technique speeds up the inverse computation using the matrix symmetry.
    Computational cost:  O(1/3 n^3)
*/
Matrix inv_chol(Matrix A) {
    unsigned int n = A.rows;

    Matrix C = cholesky(A);
    Matrix cA = Matrix(n);
    vector<double> x(n);
    
    for (unsigned int j = 0; j < n; j++) {
        x = vector<double>(n, 0.);
        x[j] = 1;
        x = solve_lower(C, x);
        x = solve_upper(C.t(), x);
        for (unsigned int i = 0; i < n; i++)
            cA[i][j] = x[i];
    }

    return cA;
}

Matrix inverse(Matrix A) {
    if (A.isSPD())  return inv_chol(A);
    else  return inv_lu(A);
}


/*
    Linear systems solver for triangular matrices.
    Used only for matrix inversion, better approaches can be found in linear.cpp
*/
vector<double> solve_upper(Matrix U, vector<double> b) {
    unsigned int i, j, n = U.rows;
    vector<double> x(n, 0);
    double sp;

    x[n-1] = b[n-1] / U[n-1][n-1];
    for (i = n-2; i < n; i--) {
        sp = 0;
        for (j = i+1; j < n; j++)
            sp += U[i][j]*x[j];
        
        x[i] = (b[i] - sp) / U[i][i];
    }
    return x;
}

vector<double> solve_lower(Matrix L, vector<double> b) {
    unsigned int i, j, n = L.rows;
    vector<double> x(n, 0);
    double sp;

    x[0] = b[0] / L[0][0];
    for (i = 1; i < n; i++) {
        sp = 0;
        for (j = 0; j < i; j++)
            sp += L[i][j]*x[j];

        x[i] = (b[i] - sp) / L[i][i];
    }
    return x;
}

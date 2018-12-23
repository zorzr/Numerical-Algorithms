#include <math.h>
#include "linear.hpp"
#define SWAP(a,b)  do { typeof(a) tmp = a;   a = b;   b = tmp; } while(0);

/*
    Gaussian Elimination: direct method for solving linear systems.
    Partial pivoting is used. Similar to LU factorization.
    Computational cost:  O(2/3 n^3)
*/
vector<double> gaussian_elimination(Matrix A, vector<double> b) {
    unsigned int i, j, k, n = b.size();
    vector<double> x = b;
    Matrix U = A.copy();
    double l, max;

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
        SWAP(x[j], x[k])
        
        for (i = k+1; i < n; i++) {
            l = U[i][k] / U[k][k];
            x[i] -= l * x[k];
            for (j = k; j < n; j++)
                U[i][j] -= l * U[k][j];
        }
    }

    return solve_upper(U, x);
}

/*
    Jacobi method: iterative algotithm for solving linear systems.
    Splitting methods represent the matrix A as:     A = M - N
     then at each iteration they solve    Mx = b + Nx
    Here the matrix M only contains the diagonal of A.
    Computational cost:  O(2 n^2) for each iteration
*/
vector<double> jacobi(Matrix A, vector<double> b, double tol, int maxits, vector<double> x0) {
    vector<double> x = x0;
    vector<double> r = sub(b, A*x0);
    double res0 = norm(r);

    Matrix itM = Matrix(A.rows);
    for (unsigned int i = 0; i < A.rows; i++)
        itM[i][i] = A[i][i];
    itM = inverse(itM);

    for (int k = 0; k < maxits; k++) {
        x = add(x, itM*r);
        r = sub(b, A*x);

        if (norm(r)/res0 < tol)
            break;
    }

    return x;
}

/*
    Gauss-Seidel method: iterative algotithm for solving linear systems.
    Splitting method similar to Jacobi's, but here M is a matrix with the lower
     triangular values of A.
    Computational cost:  O(2 n^2) for each iteration
*/
vector<double> gauss_seidel(Matrix A, vector<double> b, double tol, int maxits, vector<double> x0) {
    vector<double> x = x0;
    vector<double> r = sub(b, A*x0);
    double res0 = norm(r);

    Matrix itM = Matrix(A.rows);
    for (unsigned int i = 0; i < A.rows; i++)
        for (unsigned int j = 0; j <= i; j++)
            itM[i][j] = A[i][j];
    itM = inverse(itM);

    for (int k = 0; k < maxits; k++) {
        x = add(x, itM*r);
        r = sub(b, A*x);

        if (norm(r)/res0 < tol)
            break;
    }

    return x;
}

// Needed to make things easier
vector<double> add(vector<double> x, vector<double> y) {
    vector<double> z(x.size());
    for (unsigned int i = 0; i < x.size(); i++)
        z[i] = x[i] + y[i];
    return z;
}
vector<double> sub(vector<double> x, vector<double> y) {
    vector<double> z(x.size());
    for (unsigned int i = 0; i < x.size(); i++)
        z[i] = x[i] - y[i];
    return z;
}
double norm(vector<double> x) {
    double ssv = 0;
    for (unsigned int i = 0; i < x.size(); i++)
        ssv += pow(x[i], 2);
    return sqrt(ssv);
}
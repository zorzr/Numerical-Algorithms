#include <math.h>
#include "linear.hpp"


/*
    Gaussian Elimination: direct method for solving linear systems.
    Partial pivoting is used. Similar to LU factorization.
    Computational cost:  O(2/3 n^3)
*/
vector<double> gaussian_elimination(Matrix& A, vector<double> b) {
    unsigned int i, j, k, n = b.size();
    vector<double> x = b;
    Matrix U = A.copy();
    double l, max;

    for (k = 0; k < n-1; k++) {
		max = 0;
		for (i = k; i < n; i++) {
			if (abs(U(i,k)) > max) {
				j = i;
				max = abs(U(i,k));
			}
		}
        
        for (i = k; i < n; i++)
            swap(U(j,i), U(k,i));
        swap(x[j], x[k]);
        
        for (i = k+1; i < n; i++) {
            l = U(i,k) / U(k,k);
            x[i] -= l * x[k];
            for (j = k; j < n; j++)
                U(i,j) -= l * U(k,j);
        }
    }

    return solve_upper(&U, x);
}

/*
    Jacobi method: iterative algotithm for solving linear systems.
    Splitting methods represent the matrix A as:     A = M - N
     then at each iteration they solve    Mx = b + Nx
    Here the matrix M only contains the diagonal of A.
    Computational cost:  O(2 n^2) for each iteration
*/
vector<double> jacobi(Matrix &A, vector<double> b, double tol, int maxits, vector<double> x0) {
    Matrix B = to_col(b);
    Matrix X = to_col(x0);
    Matrix R = B - prod(A,X);
    double res0 = norm(R);

    Matrix itM = Matrix(A.rows);
    for (unsigned int i = 0; i < A.rows; i++)
        itM(i,i) = 1 / A(i,i);
    
    int k;
    for (k = 0; k < maxits; k++) {
        X = X + prod(itM, R);
        R = B - prod(A,X);

        if (norm(R)/res0 < tol)
            break;
    }

    printf("%d iterations\n", k);
    return X.get_col(0);
}

/*
    Gauss-Seidel method: iterative algotithm for solving linear systems.
    Splitting method similar to Jacobi's, but here M is a matrix with the lower
     triangular values of A. Faster convergence.
    Computational cost:  O(2 n^2) for each iteration
*/
vector<double> gauss_seidel(Matrix &A, vector<double> b, double tol, int maxits, vector<double> x0) {
    Matrix B = to_col(b);
    Matrix X = to_col(x0);
    Matrix R = B - prod(A,X);
    double res0 = norm(R);

    Matrix itM = Matrix(A.rows);
    for (unsigned int i = 0; i < A.rows; i++)
        for (unsigned int j = 0; j <= i; j++)
            itM(i,j) = A(i,j);
    itM = inverse(itM);

    int k;
    for (k = 0; k < maxits; k++) {
        X = X + prod(itM, R);
        R = B - prod(A,X);

        if (norm(R)/res0 < tol)
            break;
    }

    printf("%d iterations\n", k);
    return X.get_col(0);
}

double norm(Matrix vect) {
    double ssv = 0;
    for (unsigned long i = 0; i < vect.rows*vect.cols; i++)
        ssv += vect.M[i]*vect.M[i];
    return sqrt(ssv);
}
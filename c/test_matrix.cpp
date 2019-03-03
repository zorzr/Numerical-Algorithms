#include <stdio.h>
#include <vector>
#include "matrix.hpp"
#include "tools.h"

// Tests
void test_prod() {
    vector<vector<double> > V1 {
        { 1,  0, 2, -1, 1 },
        {-2,  1, 0,  3, 0 },
        { 0, -1, 1,  0, 1 }
    };
    vector<vector<double> > V2 {
        {  1,  1 },
        { -1,  0 },
        {  0, -2 },
        { -1,  1 },
        {  0,  1 }
    };

    Matrix A = Matrix(V1);
    Matrix B = Matrix(V2);

    printf("[PRODUCT]\n");
    print_matrix(prod(A,B));
}

void test_lu() {
    vector<vector<double> > V {
        { 1,  0,  2 },
        { 2, -1,  3 },
        { 4,  1,  8 }
    };

	Matrix A = Matrix(V);
	Matrix *L, *U, *P;
	tie(L, U, P) = lu(A);

	printf("[LU FACTORIZATION]\nL:\n");
	print_matrix(*L);
	printf("U:\n");
	print_matrix(*U);
	printf("P:\n");
	print_matrix(*P);
}

void test_inverse() {
    vector<vector<double> > V {
        { 1,  0,  2 },
        { 2, -1,  3 },
        { 4,  1,  8 }
    };

    Matrix A = Matrix(V);

    printf("[INVERSE (LU FACT.)]\n");
    print_matrix(inv_lu(A));
}

void test_chol() {
    /*
    vector<vector<double> > V {
        { 25, 15, -5 },
        { 15, 18,  0 },
        { -5,  0, 11 }
    };
    */
    vector<vector<double> > V {
        {  3, -2,  0 },
        { -2,  2,  0 },
        {  0,  0,  1 }
    };

    Matrix S = Matrix(V);
    
    printf("[INVERSE (CHOLESKY ALG.)]\n");
    print_matrix(inv_chol(S));
}


int main(int argc, char **argv) {
    test_prod();
    test_lu();
    test_inverse();
    test_chol();
    return 0;
}

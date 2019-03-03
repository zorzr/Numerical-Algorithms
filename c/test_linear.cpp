#include <stdio.h>
#include "linear.hpp"
#include "tools.h"

// Tests
void test_gem() {
    vector<vector<double> > V {
        {  2,  2,  0 },
        {  1,  1, -1 },
        {  3, -2,  4 }
    };
    vector<double> b = { 4, 1, 5 };
    Matrix A = Matrix(V);

    printf("[GAUSSIAN ELIMINATION]\n");
    print_vect(gaussian_elimination(A, b));
}

void test_jacobi() {
    vector<vector<double> > V {
        {  7,  6,  3 },
        {  2,  5, -4 },
        { -4, -3,  8 }
    };
    vector<double> b = { 16, 3, 1 };
    Matrix A = Matrix(V);

    printf("[JACOBI METHOD]\n");
    print_vect(jacobi(A,b,0.000001,100,{0,0,0}));
}

void test_gauss_seidel() {
    vector<vector<double> > V {
        {  7,  6,  3 },
        {  2,  5, -4 },
        { -4, -3,  8 }
    };
    vector<double> b = { 16, 3, 1 };
    Matrix A = Matrix(V);

    printf("[GAUSS-SEIDEL METHOD]\n");
    print_vect(gauss_seidel(A,b,0.000001,100,{0,0,0}));
}


int main(int argc, char **argv) {
    test_gem();
    test_jacobi();
    test_gauss_seidel();
    return 0;
}
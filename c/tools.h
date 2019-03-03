#ifndef TOOLS_H
#define TOOLS_H

#include "matrix.hpp"

void print_matrix(Matrix& A) {
    for (unsigned int i = 0; i < A.rows; i++) {
        for(unsigned int j = 0; j < A.cols; j++)
            printf("%lf\t", A(i,j));
        printf("\n");
    }
    printf("\n");
}

void print_vect(vector<double> v) {
    for (unsigned int i = 0; i < v.size(); i++)
        printf("%lf\t", v[i]);
    printf("\n\n");
}

#endif
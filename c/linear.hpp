#ifndef LINEAR_HPP
#define LINEAR_HPP

#include <vector>
#include "matrix.hpp"
using namespace std;

vector<double> gaussian_elimination(Matrix A, vector<double> b);
vector<double> jacobi(Matrix A, vector<double> b, double tol, int maxits, vector<double> x0);
vector<double> gauss_seidel(Matrix A, vector<double> b, double tol, int maxits, vector<double> x0);

vector<double> add(vector<double> x, vector<double> y);
vector<double> sub(vector<double> x, vector<double> y);
double norm(vector<double> x);

#endif
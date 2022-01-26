#include "cblas.h"
#include "Matrix_sum.h"

Matrix_sum::Matrix_sum(const int N){
    n = N;
}

double Matrix_sum::Find_BLAS_Sum(double *matrix_value, const int matrix_size){

    // take in a reference to first value of matrix
    int ix = 1;
    int iy = 0;
    double y = 1.0;

    return cblas_ddot(matrix_size, matrix_value, ix, &y, iy);
}

double Matrix_sum::Find_BLAS_Sum2(std::vector<double> &matrix){

    std::vector<double> x = matrix;
    int ix = 1;
    int iy = 0;
    double y = 1.0;

    return cblas_ddot(n, &x[0], ix, &y, iy);
}

double Matrix_sum::Find_Manual_Sum(double (&matrix)[], int matrix_size){
    return 0.0;
}
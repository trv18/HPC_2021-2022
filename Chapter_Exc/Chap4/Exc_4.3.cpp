#include <iostream>
#include <vector>
#include <array>
#include "cblas.h"

#include "Exc_4.3_funcs.hpp"

void Vec2Array(std::vector<std::vector<double> > &vals, int N)
// double** Vec2Array(std::vector<std::vector<double> > &vals, int N)

{
        // NOTE: this doesnt work because pointers do not reference continous 
        // piece of memory for different rows
        
//    double** temp = new double*[N];
//    for(unsigned i=0; (i < N); i++)
//    { 
//       temp[i] = new double[N];
//       for(unsigned j=0; (j < N); j++)
//       {
//           temp[i][j] = vals[i][j];
//       } 
//    }
//    return temp;
}

int main(){
    const int N = 3;
    std::vector<double> Matr(N*N);
    std::vector<double> X(N,1);

    symmetric(Matr, N); 

    //  Print out matrix
    print_matrix(Matr, N, N);

    
    CBLAS_LAYOUT Layout = CblasRowMajor;
    CBLAS_TRANSPOSE transA = CblasTrans;
    CBLAS_TRANSPOSE transB = CblasNoTrans;

    const int M = N;
    const int K = N;
    const double alpha = 1.0;
    const double lda = N;
    double C[N][N] = {0.0}; // "A" matrix in Algorithm
    double y[N] = {0.0};    // "b" matrix in Algorithm
    const double icx = 1.0;
    const double icy = 1.0;

    
    // Inputs: ... , ... , ... , Arow, Bcols, Acols, 
    cblas_dgemm(Layout, transA, transB, M, N, K, alpha, &Matr[0], lda,  &Matr[0], lda, 1.0, &C[0][0], lda);
    cblas_dgemv(Layout, transA, M, N, alpha, &Matr[0], lda,  &X[0], icx, 1.0, &y[0], icy);



}
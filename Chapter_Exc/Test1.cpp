#include <iostream>
#include <vector>
#include "cblas.h"
#include <cmath>
#include <cstdlib>

typedef std::vector<double>dvec; 

void generate_y(const int N, dvec &y){
    for(int i=0; i<N; i++){
        y[i] = i;
    }
}

void print_matrix(std::vector<double> &vec, const int N, const int M){
     for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            // std::cout << M_ref[i][j] << std::endl;
            std::cout << vec[i*N+j] <<" ";
        }   
        std::cout << "\n";
    }
    std::cout << "\n";

}

/**
 * @brief Takes a reference to matrix of zeros A and fils it in using triangular 
 * 
 * @param N 
 * @param A 
 */
void generate_A(const int N, dvec &A){
    for(int i = 0; i < N ; ++i){
        for(int j = i; j < N; ++j){
            if(j==i){
                A.push_back(1.0);
            }
            else if(j==i+1){
                A.push_back(0.5);
            }
            else{
                A.push_back(0);
            }
        }

    }
}

int main(){

    const int N = 20;  // matrix dimension
    const double alpha = 0.9;  // alpha from algorithm;
    const int iter = 1000;  // number of iterations 



    dvec b(N,0); // initialize beta vector
    dvec y(N,0); // initialize y vector
    dvec A; // matrix of coef
    dvec X_k(N,0); // vector for X_{k} iteration
    dvec X_kp1(N,0); // vector for X_{k+1} iteration
    dvec temp_vec(N,0); // vector to be overwritten when calling BLAS
    dvec zero_vec(N,0); // vector to reinitialize temp_vec to 0


    CBLAS_LAYOUT Layout = CblasRowMajor;

    const int M = N;
    const int K = N;
    const int a = 1;
    const double lda = N;
    const double icx = 1.0;
    const double icr = 1.0;

    generate_A(N, A); // fill in A
    generate_y(N, y); // fill in y

     // Initialise X0 guess
     srand(time(NULL));
    for (int i = 0; i < N; i++){
        X_k[i] = ((double) rand() / (RAND_MAX) - 0.5 )*2;  
    }


    // calculate b vector
    cblas_dspmv(Layout,CblasUpper, N, a, &A[0],  &y[0], icx, 1.0, &b[0], icr); //dgemv overwrites b <- Ay + b



    for(int i = 0; i < iter; i++){
        
        // allow temp_vec to be overwritten with b-Ax_{k}
        temp_vec = b;
        cblas_dspmv(Layout, CblasUpper, N, -1.0, &A[0],  &X_k[0], icx, 1.0, &temp_vec[0], icr); //dgemv overwrites b <- Ay + b
        // calculate X_kp1 
        cblas_daxpy(N, alpha, &temp_vec[0], icx, &X_k[0], icr); // overwrties b2 ==> p0 = b - Ax

    }

    // calculate Error
    temp_vec = zero_vec;
    cblas_dspmv(Layout, CblasUpper, N, alpha, &A[0],  &X_k[0], icx, 1.0, &temp_vec[0], icr); //dgemv overwrites tempvec <- Ax
    
    
    cblas_daxpy(N, -1*alpha, &b[0], icx, &temp_vec[0], icr); // overwrties tempvec < b- Ax
    
    double error = cblas_dnrm2(N, &temp_vec[0],icx);
    std::cout << "error: " << error << std::endl;
}
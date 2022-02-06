#include <vector>
#include <iterator> // includes begin and end
#include <iostream> 
#include <stdlib.h>
#include <time.h>
#include "cblas.h"


void symmetric(std::vector<double> &M,  const int N){
    srand(time(NULL));

    for (int i = 0; i < N; i++){

        for(int j = i; j < N; j++){
            M[i*N+j] = M[j*N+i] = ((double) rand() / (RAND_MAX) - 0.5 )*2;  
        }
    }

    // const int alpha = 1;
    // int ix = 1;
    // int iy = 0;
    // double y = 1.0;

    // // double Mt = cblas_ddot(N*N, &M[0][0], ix, &y, iy);
    // // std::cout << "ref value: " << Mt << "\n";

    // // cblas_caxpy(N*N, &alpha ,&Mt, ix, &M, ix);    

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



  // @param [A]: matrix of coeff
  // @param [b]: vector for right hand side
  // @param [X_k]: Initial guess for solution X
  // @returns X solution overwritten in X_k
std::vector<double> Conj_Grad_Method(std::vector<double> &A, std::vector<double> &b, std::vector<double> &X_k)
  

{

    size_t N = b.size();
    int k = 0   ;

    CBLAS_LAYOUT Layout = CblasRowMajor;
    CBLAS_TRANSPOSE transA = CblasNoTrans;
    CBLAS_TRANSPOSE transB = CblasTrans;

    const int M = N;
    const int K = N;
    const double alpha = 1.0;
    const double lda = N;
    const double icx = 1.0;
    const double icr = 1.0;

    std::vector<double> temp_vec(N,0); // vector to be overwritten when calling BLAS
    std::vector<double> zero_vec(N,0); // vector to reinitialize temp_vec to 0
    std::vector<double> r_k;           // "r_{k}" vector in Algorithm
    std::vector<double> r_kp1;         // "r_{k+1}" vector in Algorithm
    std::vector<double> p_k;           // "p_{k}" vector in Algorithm
    std::vector<double> b2 = b;        // copy of b matrix
    std::vector<double> all_r;         // Vector all error norms.

    double a_k;               // alpha_{k} in Algorithm
    double b_k;               // beta{k} in Algorithm

    // Initialise r0, p0
    cblas_dgemv(Layout, transA, M, N, alpha, &A[0], lda,  &X_k[0], icx, 1.0, &temp_vec[0], icr); //dgemv overwrites r 
    cblas_daxpy(N, -1*alpha, &temp_vec[0], icx, &b2[0], icr); // overwrties b2 ==> p0 = b - Ax
    r_k = r_kp1 = p_k = b2;
    temp_vec = zero_vec;


    while(cblas_dnrm2(N, &r_kp1[0],icx) > 1e-3 && k<10000) {
        r_k = r_kp1;
        a_k = cblas_ddot(N, &r_k[0], icx, &r_k[0], icx); // outputs double precision

        cblas_dgemv(Layout, transA, M, N, alpha, &A[0], lda,  &p_k[0], icx, 1.0, &temp_vec[0], icr); //dgemv overwrites temp_vec

        a_k /= cblas_ddot(N, &p_k[0], icx, &temp_vec[0], icx); // calcualte alpha_{k}
        temp_vec = zero_vec;

        cblas_daxpy(N, a_k, &p_k[0], icx, &X_k[0], icx); // Overwrite X_{k} with X_{k+1} using alpha_{k}*p_{k} + X_{k}
        
        cblas_dgemv(Layout, transA, M, N, alpha, &A[0], lda,  &p_k[0], icx, 1.0, &temp_vec[0], icr); //dgemv overwrites temp_vec with A*p_k
        cblas_daxpy(N, -1*a_k, &temp_vec[0], icx, &r_kp1[0], icx); // calculate r_{k+1}
        
        b_k =  cblas_ddot(N, &r_kp1[0], icx, &r_kp1[0], icx) / cblas_ddot(N, &r_k[0], icx, &r_k[0], icx); //calculate b_{k}
        
        cblas_dscal(N,b_k, &p_k[0],icx);
        cblas_daxpy(N, 1, &r_kp1[0] , icx, &p_k[0], icx); // Overwrite r_{k} with X_{k+1} using beta_{k}*p_{k} + r_{k}
        all_r.push_back(cblas_dnrm2(N, &r_kp1[0],icx));
        
        k++; // add 1 to iteration
    }

    std::cout << "Error: " <<cblas_dnrm2(N, &r_kp1[0] ,icx) <<"\n";
    std::cout << "Iteration: " << k <<"\n";

    return all_r;

}
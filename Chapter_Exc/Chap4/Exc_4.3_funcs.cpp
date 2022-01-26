#include <vector>
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


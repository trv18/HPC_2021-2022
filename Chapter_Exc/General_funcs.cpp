#include <vector>
#include <iostream>
#include <iomanip>      // std::setw


// void print_matrix(std::vector<double> &vec, const int N, const int M){}

void print_matrix_banded(std::vector<double> &vec, const int N, const int M, const int ku, const int kl){
    for(int i=0; i<N; i++){

        for(int j=0; j<M; j++){

            if(std::max(0, i-ku-1) <= j && j <= std::min(N-1, i+kl-1)){
                std::cout << std::setw(10) << vec[j* + i*(ku+j-i)] << " ";
            } 
            else{
                std::cout << 0 << " ";
            }
        }
        std::cout << "\n";
    }
}

// 
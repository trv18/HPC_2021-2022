#include <iostream>
#include "cblas.h"
#include <vector>
#include <cmath>

typedef std::vector<double> dvec;

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

double f_func(double x){

    return sin(M_PI*x);
}

void generate_f(dvec &f, const int N, const double h){
    for(int i=0; i < N; i++){
        f[i] =  f_func(i*h);
    }
}

/* generates D matrix using triangular sroage Conj_Grad_Method
 * @param D: reference of vector of 0s D to be filled in.
 * @param N: matrix size NxN 
 */ 
void generate_D(dvec &D, const int N, const double h){
    dvec options = {2.0, -5.0, 4.0, -1.0};
    for(int i=0; i < N; i++){
        int k = 0;
        for(int j=i; j < N; j++){
            D.push_back(1/pow(h,2)*options[k]);
            k++;
            
        }
    }
}

void generate_b(dvec &b, const int N, const double h){  
    const double L = (N-1)*h;
    b[N-3] = -f_func(L+h);
    b[N-2] = 4*f_func(L+h)  -  f_func(L+2*h);
    b[N-1] = -5*f_func(L+h) + 4*f_func(L+2*h) - f_func(L+3*h);
}

int main(){
    const int N = 5;
    const double h = 0.001;

    dvec D;
    dvec b(N, 0);
    dvec f(N);

    generate_D(D, N, h);
    generate_b(b, N, h);
    generate_f(f, N, h);

    print_matrix(D, 1, (N*N+N)/2);
}
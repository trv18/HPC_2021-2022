#include <iostream>
#include <vector>
#include <array>
#include <stdlib.h>
#include <time.h>
#include "cblas.h"
#include <cmath>

#include "Exc_4.4_funcs.hpp"
#include "Exc_4.3_funcs.hpp"
#include "boost/program_options.hpp"
namespace po = boost::program_options;

#define lambda 1

// g++ Exc_4.3.cpp Exc_4.3_funcs.cpp -lblas -I/usr/include/python3.8 
// -I/usr/lib/python3/dist-packages/numpy/core/include/numpy/ -lpython3.8 -o Exc_4.3

void Vec2Array(std::vector<std::vector<double> > &vals, int N)
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


double func_f(double x){ std::cout << x; return -(lambda + pow(M_PI,2))*sin(M_PI*x);  }
// double func_u(double x){ return (cos();  }


int main(int argc, char* argv[]){

    
 // ------------------------------ Command Line Inputs ----------------------------------------
    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
    ("version,v", "print version string")
    ("help", "produce help message")
    ;

    // Declare a group of options that will be 
    // allowed both on command line and in 
    // config file
    po::options_description config("Configuration");
    config.add_options()
    ("SysDim,N", po::value<int>()->default_value(5),
    "Square Matrix Dimension.")
    ("Print_Matrix,P", po::value<bool>()->default_value(false), 
    "Output Matrices to be solved")
    ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);


    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);

    notify(vm);

    if (vm.count("help")) {

        std::cout << "Allows Options: " << std::endl;
        std::cout << generic << std::endl;
        std::cout << config << std::endl;

        return 0;
    }

    const int n = vm["SysDim"].as<int>();
    bool pm = vm["Print_Matrix"].as<bool>();
    std::cout << "Chosen Dimension: "<< n << std::endl;


    // ------------------------------- Blas Functions -------------------------------------------
    srand(time(NULL));
    
    const int a = 0, b = 2;
    const double dx = (b-a)/(double)(n-1);
    // const double lambda = 1.0;
    const double alpha = -2*pow(dx,-2) - lambda;
    const double beta  = pow(dx,-2);
    std::vector<double> Matr; // "A" matrix in Algorithm
    std::vector<double> f_vec; // "b" vector in Algorithm
    std::vector<double> X(n-2); // initial guess. using vector.push_back so no initial size


    // auto func_f = [&](double x){return (-lambda + pow(M_PI,2))*sin(M_PI*x); };

    SymmetricStored(Matr, n, alpha, beta);      
    create_f(f_vec, n, dx, a, b, *func_f); // u0 & un-1 = 0
    

    // Initialise X0 guess
    for (int i = 0; i < n-2; i++){
        X[i] = ((double) rand() / (RAND_MAX) - 0.5 )*2;  
    }
    
    
    std::vector<double> Errors = Conj_Grad_Method_v2(Matr, f_vec, X);

    if(pm){

        print_matrix(Matr, 1, Matr.size());
        print_matrix(f_vec, 1, n-2);

        print_matrix(X, 1, X.size());
        print_matrix(Errors, 1, Errors.size());
        std::cout << sin((a+dx)*M_PI) << " " << a+dx << "\n";
            }
    




}
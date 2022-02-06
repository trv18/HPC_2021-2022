#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <cblas.h>

#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include "Exc_4.3_funcs.hpp"

#define F77NAME(x) x##_
extern "C" {
void F77NAME(dgesv)(const int &n, const int &nrhs, double *a,
                    const int &lda, int *ipvt, double* B,
                    const int &ldb, int& info );

void F77NAME(dgemv)(const char &   trans, const int & m, 
                    const int &    n, const double &  alpha,
                    const double * A, const int &     lda,
                    const double * x, const int &     incx,
                    const double & beta, double *     y, 
                    const int &    incy);

}

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
    ("Matrix_Dimension,N", po::value<int>()->default_value(5),
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

    const int N = vm["Matrix_Dimension"].as<int>();
    bool pm = vm["Print_Matrix"].as<bool>();
    std::cout << "Chosen Dimension: "<< N << std::endl;


    // ------------------------------- Blas Functions -------------------------------------------

    srand(time(NULL));

    std::vector<double> Matr(N*N); // "A" matrix in Algorithm
    std::vector<double> X(N);

    symmetric(Matr, N); 

    // Initialise actual solution X
    for (int i = 0; i < N; i++){
        X[i] = ((double) rand() / (RAND_MAX) - 0.5 )*2;  
    }


    //  Print out matrix
    if (pm){
        print_matrix(Matr, N, N);
        print_matrix(X,1,N);
    }

    const int M = N;
    const int K = N;
    const double alpha = 1.0;
    const double lda = N;
    std::vector<double> C(N*N,0) ;         // "A" matrix in Algorithm
    std::vector<double> y(N,0);    // "b" matrix in Algorithm
    const double icx = 1.0;
    const double icy = 1.0;

    CBLAS_LAYOUT Layout = CblasRowMajor;
    CBLAS_TRANSPOSE transA = CblasTrans;
    CBLAS_TRANSPOSE transB = CblasNoTrans;
    CBLAS_SIDE sideA = CblasLeft;

    
    // Inputs: ... , ... , ... , Ar ow, Bcols, Acols, 
    // cblas_dgemm(Layout, transA, transB, M, N, K, alpha, &Matr[0], lda,  &Matr[0], lda, 1.0, &C[0], lda); // overwrites C
    // cblas_dsymm (Layout, sideA, CblasUpper, M, N, alpha, &Matr[0], lda,  &Matr[0], lda, 1.0, &C[0], lda); // overwrites C
    F77NAME(dgemv)('N', N, M, alpha, &Matr[0], lda, &X[0], icx, 0.0, &y[0], icy); // overwrites y

    // Initialise X0 guess
    for (int i = 0; i < N; i++){
        X[i] = ((double) rand() / (RAND_MAX) - 0.5 )*2;  
    }

    std::vector<double> Errors = Conj_Grad_Method(Matr, y, X);

    if(pm){
        print_matrix(X,1,N);
    }

    if(pm){
        print_matrix(Errors, 1, Errors.size());
    }

}